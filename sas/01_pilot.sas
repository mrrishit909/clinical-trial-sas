/*  CDISC pilot study CDISCPILOT01: xanomeline skin patch for mild-to-moderate Alzheimer's disease, 254 patients,
    placebo vs low dose (54 mg) vs high dose (81 mg), 24 weeks. Re-analysed from its public ADaM datasets the way a
    submission's tables are made:
      1. who was in the study (demographics, intent-to-treat population)
      2. who finished it (disposition; who stopped because of side effects)
      3. did it work: ADAS-Cog(11) change from baseline to week 24, LOCF, ANCOVA (the protocol's primary endpoint)
      4. was it safe: time to the first skin reaction (Kaplan-Meier, log-rank, Cox), adverse events, glucose
    Ran in SAS Enterprise Guide (SAS 9.4) on USF's virtual desktop; the log is sas/01_pilot.log. Lines starting CSV|
    in the log are the results; the repository's check.py recomputes them in Python and compares them with the
    R Consortium's FDA submission pilot outputs for the same study.                                              */

options nodate nonumber validvarname=v7;
%put NOTE: SAS &sysvlong on &sysscp &sysscpl;

/* 1. Read the ADaM transport files from the PhUSE repository, pinned to one commit */
%let src = https://raw.githubusercontent.com/phuse-org/phuse-scripts/fbd239d9450a99621f0199323d7abf3ae22036db/data/adam/cdiscpilot01;
%macro get(ds);                              /* filerefs and librefs are 8 characters at most: reuse f and x */
  filename f temp;
  proc http url="&src/&ds..xpt" out=f; run;
  libname x xport "%sysfunc(pathname(f))";
  data &ds; set x.&ds; run;
  libname x clear; filename f clear;
%mend;
%get(adsl) %get(adqsadas) %get(adtte) %get(adae) %get(adlbc)

/* Results go to the log as CSV|tag|value|value... lines (explicit separators, so blank values keep their place) */
%macro csv(ds, tag, vars, where=1);
  data _null_;
    set &ds;
    where &where;
    length line $500;
    line = cats('CSV|', "&tag"
    %let i = 1;
    %do %while(%length(%scan(&vars, &i, %str( ))));
      , '|', %scan(&vars, &i, %str( ))
      %let i = %eval(&i + 1);
    %end;
    );
    putlog line;
  run;
%mend;

proc format;
  value trt 0 = 'Placebo' 54 = 'Xanomeline low dose' 81 = 'Xanomeline high dose';
run;

/* 2. Who was in the study: intent-to-treat population */
ods output summary=demo;
proc means data=adsl n mean std median min max stackodsoutput;
  where ittfl = 'Y';
  class trt01pn;
  var age heightbl weightbl bmibl mmsetot;
run;
%csv(demo, demo, trt01pn variable n mean stddev median min max)
proc freq data=adsl noprint;
  where ittfl = 'Y';
  tables trt01pn * agegr1 / out=agegrp;
  tables trt01pn * race / out=race;
run;
%csv(agegrp, agegr1, trt01pn agegr1 count)
%csv(race, race, trt01pn race count)

/* 3. Who finished: reason for leaving the study, and a test of stopping for side effects by arm */
proc freq data=adsl noprint;
  tables trt01pn * dcdecod / out=disp;
run;
%csv(disp, disp, trt01pn dcdecod count)
data adsl2;
  set adsl;
  ae_stop = (dsraefl = 'Y');
run;
ods output FishersExact=fisher;
proc freq data=adsl2;
  tables trt01pn * ae_stop / fisher;
run;
%csv(fisher, fisher, name1 nvalue1)

/* 4. Did it work: ADAS-Cog(11) change from baseline to week 24 (LOCF), efficacy population */
data adas;
  set adqsadas;
  where efffl = 'Y' and ittfl = 'Y' and paramcd = 'ACTOT' and anl01fl = 'Y' and avisitn = 24;
run;
ods output summary=adas_desc;
proc means data=adas n mean std median min max stackodsoutput;
  class trtpn;
  var base aval chg;
run;
%csv(adas_desc, adas_desc, trtpn variable n mean stddev median min max)

/* dose response: treatment as a continuous dose (0, 54, 81), site group as a factor, baseline as a covariate */
ods output ParameterEstimates=dose;
proc glm data=adas;
  class sitegr1;
  model chg = trtpn sitegr1 base / solution;
run; quit;
%csv(dose, dose, parameter estimate stderr probt, where=upcase(parameter) = 'TRTPN')

/* pairwise: treatment as a factor */
ods output LSMeans=lsm LSMeanCL=lsmcl Estimates=est;
proc glm data=adas;
  class trtpn sitegr1;
  model chg = trtpn sitegr1 base / clparm;
  lsmeans trtpn / stderr cl;
  estimate 'Low - Placebo'  trtpn -1 1 0;
  estimate 'High - Placebo' trtpn -1 0 1;
  estimate 'High - Low'     trtpn 0 -1 1;
run; quit;
%csv(lsm, adas_lsm, trtpn lsmean stderr)
%csv(est, adas_diff, parameter estimate stderr probt lowercl uppercl)

/* mean change at each visit, for the chart */
data adas_t;
  set adqsadas;
  where efffl = 'Y' and ittfl = 'Y' and paramcd = 'ACTOT' and anl01fl = 'Y' and avisitn in (8, 16, 24);
run;
proc means data=adas_t noprint nway;
  class trtpn avisitn;
  var chg;
  output out=adas_time n=n mean=mean stderr=se;
run;
%csv(adas_time, adas_time, trtpn avisitn n mean se)

/* 5. Was it safe: time to the first dermatologic event, safety population */
data tte;
  set adtte;
  where saffl = 'Y' and paramcd = 'TTDE';
run;
ods output HomTests=logrank Quartiles=quart;
proc lifetest data=tte plots=survival(atrisk=0 to 200 by 20);
  time aval * cnsr(1);
  strata trtan;
  format trtan trt.;
run;
%csv(logrank, logrank, test chisq df probchisq)
%csv(quart, km_median, trtan estimate lowerlimit upperlimit, where=percent = 50)
data grid;
  do t = 0 to 200 by 20; output; end;
run;
proc sql;
  create table atrisk as
    select a.trtan, g.t, sum(a.aval >= g.t) as at_risk
      from grid as g, tte as a
     group by a.trtan, g.t
     order by a.trtan, g.t;
quit;
%csv(atrisk, atrisk, trtan t at_risk)
ods output ParameterEstimates=cox;
proc phreg data=tte;
  class trtan(ref='0');
  model aval * cnsr(1) = trtan / risklimits;
run;
%csv(cox, cox, classval0 estimate stderr probchisq hazardratio hrlowercl hruppercl)

/* adverse events that started on treatment: subjects with at least one, by arm */
proc sql;
  create table saf_n as select trt01an as trtan, count(*) as n from adsl where saffl = 'Y' group by trt01an;
  create table ae_any as select trtan, count(distinct usubjid) as subjects from adae
    where saffl = 'Y' and trtemfl = 'Y' group by trtan;
  create table ae_derm as select trtan, count(distinct usubjid) as subjects from adae
    where saffl = 'Y' and trtemfl = 'Y' and cq01nam ne '' group by trtan;
  create table ae_soc as select aebodsys, trtan, count(distinct usubjid) as subjects from adae
    where saffl = 'Y' and trtemfl = 'Y' group by aebodsys, trtan;
  create table ae_pt as select aedecod, trtan, count(distinct usubjid) as subjects from adae
    where saffl = 'Y' and trtemfl = 'Y' group by aedecod, trtan;
quit;
%csv(saf_n, saf_n, trtan n)
%csv(ae_any, ae_any, trtan subjects)
%csv(ae_derm, ae_derm, trtan subjects)
%csv(ae_soc, ae_soc, aebodsys trtan subjects)
%csv(ae_pt, ae_pt, aedecod trtan subjects, where=subjects >= 5)

/* glucose: change from baseline to week 20, high dose vs placebo, ANCOVA with baseline as a covariate */
data gluc;
  set adlbc;
  where paramcd = 'GLUC' and avisitn = 20 and trtpn in (0, 81) and chg ne . and base ne .;
run;
ods output LSMeanCL=g_lsmcl Estimates=g_est FitStatistics=g_fit;
proc glm data=gluc;
  class trtpn;
  model chg = trtpn base / clparm;
  lsmeans trtpn / cl;
  estimate 'High - Placebo' trtpn -1 1;
run; quit;
%csv(g_est, gluc_diff, parameter estimate stderr probt lowercl uppercl)
%csv(g_fit, gluc_fit, rootmse)
%csv(g_lsmcl, gluc_lsm, trtpn lsmean lowercl uppercl)

/* 6. Charts */
data adas_plot;
  set adas_time;
  lo = mean - 1.96 * se; hi = mean + 1.96 * se;
  format trtpn trt.;
run;
proc sgplot data=adas_plot;
  title 'ADAS-Cog(11): mean change from baseline (LOCF); higher = worse';
  series x=avisitn y=mean / group=trtpn markers name='s';
  scatter x=avisitn y=mean / group=trtpn yerrorlower=lo yerrorupper=hi;
  refline 0 / axis=y;
  xaxis values=(8 16 24) label='Week'; yaxis label='Mean change (95% CI)';
  keylegend 's' / title='';
run;
data disp_plot;
  set disp;
  length reason $24;
  reason = ifc(dcdecod = 'COMPLETED', '1 Completed', ifc(dcdecod = 'ADVERSE EVENT', '2 Stopped: adverse event', '3 Stopped: other'));
  format trt01pn trt.;
run;
proc sgplot data=disp_plot pctlevel=group;
  title 'How patients left the study';
  hbar trt01pn / response=count group=reason groupdisplay=stack stat=percent;
  xaxis label='Percent of randomised patients'; yaxis display=(nolabel);
run;
title;
