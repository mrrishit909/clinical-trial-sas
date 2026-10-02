"""The project's check: every SAS result recomputed in Python from the same ADaM files, and the published figures
compared with the R Consortium's FDA submission pilot 1 outputs for the same study
(github.com/RConsortium/submissions-pilot1, output/: tlf-demographic.out, tlf-primary.rtf, tlf-efficacy.rtf,
tlf-kmplot.pdf).

    ./venv/bin/python check.py        (needs data/*.xpt: see README; results/ from parse_log.py)

Tolerances: descriptive statistics and least-squares results to 1e-6; Cox model to 1e-4 (two optimisers).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact, t as t_dist
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import SurvfuncRight, survdiff

HERE = Path(__file__).parent
R = HERE / "results"
ARMS = {"Placebo": 0, "Xanomeline low dose": 54, "Xanomeline high dose": 81}
CODES = ["TRT01PN", "TRT01AN", "TRTPN", "TRTAN", "CNSR", "AVISITN"]

# R Consortium submission pilot 1 (placebo, low dose, high dose)
REF_DEMO = {"AGE": [(75.21, 8.59, 76, 52, 89), (75.67, 8.29, 77.5, 51, 88), (74.38, 7.89, 76, 56, 88)],
            "HEIGHTBL": [(162.57, 11.52, 162.6, 137.2, 185.4), (163.43, 10.42, 162.6, 135.9, 195.6), (165.82, 10.13, 165.1, 146.1, 190.5)],
            "WEIGHTBL": [(62.76, 12.77, 60.55, 34, 86.2), (67.28, 14.12, 64.9, 45.4, 106.1), (70, 14.65, 69.2, 41.7, 108)],
            "BMIBL": [(23.64, 3.67, 23.4, 15.1, 33.3), (25.06, 4.27, 24.3, 17.7, 40.1), (25.35, 4.16, 24.8, 13.7, 34.5)],
            "MMSETOT": [(18.05, 4.27, 19.5, 10, 23), (17.87, 4.22, 18, 10, 24), (18.51, 4.16, 20, 10, 24)]}
REF_COUNTS = {("agegr1", "<65"): [14, 8, 11], ("agegr1", "65-80"): [42, 47, 55], ("agegr1", ">80"): [30, 29, 18],
              ("race", "WHITE"): [78, 78, 74], ("race", "BLACK OR AFRICAN AMERICAN"): [8, 6, 9]}
REF_ADAS = {"BASE": [(79, 24.1, 12.19), (81, 24.4, 12.92), (74, 21.3, 11.74)],
            "AVAL": [(79, 26.7, 13.79), (81, 26.4, 13.18), (74, 22.8, 12.48)],
            "CHG": [(79, 2.5, 5.80), (81, 2.0, 5.55), (74, 1.5, 4.26)]}
REF_DOSE_P = 0.245
REF_DIFF = {"Low - Placebo": (-0.5, 0.82, 0.569, -2.1, 1.1), "High - Placebo": (-1.0, 0.84, 0.233, -2.7, 0.7),
            "High - Low": (-0.5, 0.84, 0.520, -2.2, 1.1)}
REF_ATRISK = {0: [86, 75, 65, 59, 50, 47, 45, 42, 40, 35, 0], 54: [84, 58, 31, 20, 14, 12, 8, 6, 6, 5, 0],
              81: [84, 48, 31, 14, 7, 4, 4, 4, 4, 3, 0]}
REF_GLUC = (0.07, 0.822)          # difference in LS means and p; its CI and RMSE use n and a normal quantile, see README


def rd(name):
    d = pd.read_sas(HERE / "data" / f"{name}.xpt", format="xport", encoding="latin-1")
    for c in CODES:
        if c in d:
            d[c] = d[c].round()                # transport-file zeros arrive as 5e-79
    return d


def res(tag):
    d = pd.read_csv(R / f"{tag}.csv", na_values=["."])            # SAS writes missing as "."
    for c in ["trt01pn", "trtpn", "trtan"]:
        if c in d and d[c].dtype == object:
            d[c] = d[c].map(lambda v: ARMS.get(v, v)).astype(float)
    return d


def close(a, b, tol=1e-6):
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b)))


def describe(d, by, var):
    g = d.groupby(by)[var]
    return pd.DataFrame({"n": g.count(), "mean": g.mean(), "stddev": g.std(), "median": g.median(), "min": g.min(), "max": g.max()})


def ls_means(fit, data, arm_col, arms, other_factors):
    """LS means as SAS computes them: covariates at their mean, other class effects averaged with equal weights."""
    x = fit.model.data.orig_exog
    out = {}
    for a in arms:
        rows = []
        levels = [sorted(data[f].unique()) for f in other_factors] or [[None]]
        for combo in pd.MultiIndex.from_product(levels) if other_factors else [()]:
            row = pd.Series(0.0, index=x.columns)
            row["Intercept"] = 1.0
            row["BASE"] = data["BASE"].mean()
            name = f"C({arm_col})[T.{a}]"
            if name in row.index:
                row[name] = 1.0
            for f, lv in zip(other_factors, combo):
                name = f"C({f})[T.{lv}]"
                if name in row.index:
                    row[name] = 1.0
            rows.append(row)
        L = pd.concat(rows, axis=1).mean(axis=1)
        out[a] = (float(L @ fit.params), float(np.sqrt(L @ fit.cov_params() @ L)))
    return out


def main():
    n = 0
    adsl, adas, adtte, adae, adlbc = (rd(x) for x in ["adsl", "adqsadas", "adtte", "adae", "adlbc"])
    itt = adsl[adsl.ITTFL == "Y"]

    # demographics
    demo = res("demo")
    for var in REF_DEMO:
        py = describe(itt, "TRT01PN", var)
        for i, arm in enumerate([0, 54, 81]):
            s = demo[(demo.trt01pn == arm) & (demo.variable.str.upper() == var)].iloc[0]
            for k in ["n", "mean", "stddev", "median", "min", "max"]:
                assert close(s[k], py.loc[arm, k]), (var, arm, k, s[k], py.loc[arm, k]); n += 1
            m, sd, med, lo, hi = REF_DEMO[var][i]
            assert round(s["mean"], 2) == m and round(s["stddev"], 2) == sd and close(s["median"], med) and close(s["min"], lo) and close(s["max"], hi), (var, arm, s.to_dict())
    for tag, col in [("agegr1", "AGEGR1"), ("race", "RACE")]:
        sas = res(tag).set_index(["trt01pn", tag])["count"]
        py = itt.groupby(["TRT01PN", col]).size()
        assert sas.sort_index().to_dict() == {k: v for k, v in py.sort_index().to_dict().items() if v > 0}, tag; n += len(sas)
        for (t, lv), want in REF_COUNTS.items():
            if t == tag:
                assert [int(sas.get((a, lv), 0)) for a in [0, 54, 81]] == want, (tag, lv)

    # disposition
    disp = res("disp").set_index(["trt01pn", "dcdecod"])["count"]
    assert disp.sort_index().to_dict() == adsl.groupby(["TRT01PN", "DCDECOD"]).size().sort_index().to_dict(); n += len(disp)
    fisher = res("fisher").set_index("name1")["nvalue1"]
    tab = pd.crosstab(adsl.TRT01PN, adsl.DSRAEFL.eq("Y")).to_numpy()
    assert close(fisher["P_TABLE"], fisher_exact(tab).statistic) and fisher["XP2_FISH"] < 1e-6, fisher.to_dict(); n += 1

    # primary endpoint
    a = adas[(adas.EFFFL == "Y") & (adas.ITTFL == "Y") & (adas.PARAMCD == "ACTOT") & (adas.ANL01FL == "Y") & (adas.AVISITN == 24)].copy()
    ad = res("adas_desc")
    for var in ["BASE", "AVAL", "CHG"]:
        py = describe(a, "TRTPN", var)
        for i, arm in enumerate([0, 54, 81]):
            s = ad[(ad.trtpn == arm) & (ad.variable.str.upper() == var)].iloc[0]
            for k in ["n", "mean", "stddev", "median", "min", "max"]:
                assert close(s[k], py.loc[arm, k]), (var, arm, k); n += 1
            rn, rm, rsd = REF_ADAS[var][i]
            assert s["n"] == rn and round(s["mean"], 1) == rm and round(s["stddev"], 2) == rsd, (var, arm, s.to_dict())
    m1 = smf.ols("CHG ~ TRTPN + C(SITEGR1) + BASE", a).fit()
    dose = res("dose").iloc[0]
    for k, v in [("estimate", m1.params["TRTPN"]), ("stderr", m1.bse["TRTPN"]), ("probt", m1.pvalues["TRTPN"])]:
        assert close(dose[k], v), ("dose", k); n += 1
    assert round(dose["probt"], 3) == REF_DOSE_P
    m2 = smf.ols("CHG ~ C(TRTPN) + C(SITEGR1) + BASE", a).fit()
    lsm = ls_means(m2, a, "TRTPN", [0.0, 54.0, 81.0], ["SITEGR1"])
    for _, s in res("adas_lsm").iterrows():
        est, se = lsm[s["trtpn"]]
        assert close(s["lsmean"], est) and close(s["stderr"], se), ("lsm", s.to_dict(), est, se); n += 2
    df = m2.df_resid
    for _, s in res("adas_diff").iterrows():
        hi, lo = {"Low - Placebo": (54.0, 0.0), "High - Placebo": (81.0, 0.0), "High - Low": (81.0, 54.0)}[s["parameter"]]
        L = pd.Series(0.0, index=m2.params.index)
        for arm, sign in [(hi, 1), (lo, -1)]:
            if arm:
                L[f"C(TRTPN)[T.{arm}]"] += sign
        est, se = float(L @ m2.params), float(np.sqrt(L @ m2.cov_params() @ L))
        p = 2 * t_dist.sf(abs(est / se), df)
        ci = (est - t_dist.ppf(0.975, df) * se, est + t_dist.ppf(0.975, df) * se)
        for k, v in [("estimate", est), ("stderr", se), ("probt", p), ("lowercl", ci[0]), ("uppercl", ci[1])]:
            assert close(s[k], v), ("diff", s["parameter"], k); n += 1
        r = REF_DIFF[s["parameter"]]
        assert (round(s["estimate"], 1), round(s["stderr"], 2), round(s["probt"], 3), round(s["lowercl"], 1), round(s["uppercl"], 1)) == r, (s.to_dict(), r)
    at = adas[(adas.EFFFL == "Y") & (adas.ITTFL == "Y") & (adas.PARAMCD == "ACTOT") & (adas.ANL01FL == "Y") & adas.AVISITN.isin([8, 16, 24])]
    g = at.groupby(["TRTPN", "AVISITN"]).CHG
    for _, s in res("adas_time").iterrows():
        key = (s["trtpn"], s["avisitn"])
        assert close(s["n"], g.count()[key]) and close(s["mean"], g.mean()[key]) and close(s["se"], g.sem()[key]), key; n += 3

    # time to first dermatologic event
    tte = adtte[(adtte.SAFFL == "Y") & (adtte.PARAMCD == "TTDE")].copy()
    tte["event"] = 1 - tte.CNSR
    lr = res("logrank").set_index("test").loc["Log-Rank"]
    chi, p = survdiff(tte.AVAL.to_numpy(), tte.event.to_numpy(), tte.TRTAN.to_numpy())
    assert close(lr["chisq"], chi) and close(lr["probchisq"], p, 1e-4), (lr.to_dict(), chi, p); n += 2
    for _, s in res("km_median").iterrows():
        d = tte[tte.TRTAN == s["trtan"]]
        sf = SurvfuncRight(d.AVAL.to_numpy(), d.event.to_numpy())
        tt, surv = sf.surv_times, sf.surv_prob
        below = np.where(surv < 0.5)[0]
        exact = np.where(np.isclose(surv, 0.5))[0]
        if not len(below) and not len(exact):          # never falls below 50% (placebo): no median, SAS prints "."
            assert pd.isna(s["estimate"]), ("km median", s.to_dict()); n += 1
            continue
        med = (tt[exact[0]] + tt[exact[0] + 1]) / 2 if len(exact) else tt[below[0]]     # SAS: midpoint when S(t) = 0.5 exactly
        assert close(s["estimate"], med), ("km median", s.to_dict(), med); n += 1
    ar = res("atrisk")
    for arm, want in REF_ATRISK.items():
        sas = ar[ar.trtan == arm].sort_values("t")["at_risk"].astype(int).tolist()
        py = [int((tte[tte.TRTAN == arm].AVAL >= t).sum()) for t in range(0, 201, 20)]
        assert sas == py == want, (arm, sas, py, want); n += len(sas)
    x = pd.get_dummies(tte.TRTAN, prefix="t", dtype=float)[["t_54.0", "t_81.0"]]
    cox = PHReg(tte.AVAL.to_numpy(), x.to_numpy(), status=tte.event.to_numpy(), ties="breslow").fit()
    for i, (_, s) in enumerate(res("cox").sort_values("classval0").iterrows()):
        b, se = cox.params[i], cox.bse[i]
        for k, v in [("estimate", b), ("stderr", se), ("hazardratio", np.exp(b)), ("hrlowercl", np.exp(b - 1.959964 * se)), ("hruppercl", np.exp(b + 1.959964 * se))]:
            assert close(s[k], v, 1e-4), ("cox", s["classval0"], k, s[k], v); n += 1

    # adverse events
    saf = adsl[adsl.SAFFL == "Y"]
    assert res("saf_n").set_index("trtan")["n"].to_dict() == saf.groupby("TRT01AN").size().to_dict()
    te = adae[(adae.SAFFL == "Y") & (adae.TRTEMFL == "Y")]
    checks = {"ae_any": te.groupby("TRTAN").USUBJID.nunique(),
              "ae_derm": te[te.CQ01NAM != ""].groupby("TRTAN").USUBJID.nunique()}
    for tag, py in checks.items():
        assert res(tag).set_index("trtan")["subjects"].to_dict() == py.to_dict(), tag; n += len(py)
    soc = te.groupby(["AEBODSYS", "TRTAN"]).USUBJID.nunique()
    assert res("ae_soc").set_index(["aebodsys", "trtan"])["subjects"].sort_index().to_dict() == soc.sort_index().to_dict(); n += len(soc)
    pt = te.groupby(["AEDECOD", "TRTAN"]).USUBJID.nunique()
    pt = pt[pt >= 5]
    assert res("ae_pt").set_index(["aedecod", "trtan"])["subjects"].sort_index().to_dict() == pt.sort_index().to_dict(); n += len(pt)

    # glucose, week 20, high dose vs placebo
    gl = adlbc[(adlbc.PARAMCD == "GLUC") & (adlbc.AVISITN == 20) & adlbc.TRTPN.isin([0, 81]) & adlbc.CHG.notna() & adlbc.BASE.notna()].copy()
    m3 = smf.ols("CHG ~ C(TRTPN) + BASE", gl).fit()
    k = "C(TRTPN)[T.81.0]"
    gd = res("gluc_diff").iloc[0]
    ci = m3.conf_int().loc[k]
    for key, v in [("estimate", m3.params[k]), ("stderr", m3.bse[k]), ("probt", m3.pvalues[k]), ("lowercl", ci[0]), ("uppercl", ci[1])]:
        assert close(gd[key], v), ("gluc", key); n += 1
    assert (round(gd["estimate"], 2), round(gd["probt"], 3)) == REF_GLUC
    assert close(res("gluc_fit").iloc[0]["rootmse"], np.sqrt(m3.scale)); n += 1
    glsm = ls_means(m3, gl, "TRTPN", [0.0, 81.0], [])
    for _, s in res("gluc_lsm").iterrows():
        est, se = glsm[s["trtpn"]]
        q = t_dist.ppf(0.975, m3.df_resid)
        assert close(s["lsmean"], est) and close(s["lowercl"], est - q * se) and close(s["uppercl"], est + q * se), ("gluc lsm", s.to_dict()); n += 3

    print(f"OK: {n} SAS results recomputed in Python; R Consortium pilot figures matched: demographics (5 measures x 3 arms), "
          f"age-group and race counts, ADAS-Cog n/mean/SD, dose-response p {REF_DOSE_P}, all 3 pairwise comparisons, "
          f"KM numbers at risk at 11 time points x 3 arms, glucose difference and p")


if __name__ == "__main__":
    main()
