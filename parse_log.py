"""Turn the CSV| lines of the saved SAS log into one results file per tag. Column names come from the %csv calls in
the SAS program itself, so the two can't drift apart.

    ./venv/bin/python parse_log.py sas/01_pilot.log   -> results/<tag>.csv
"""
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent


def main(log):
    sas = (HERE / "sas" / "01_pilot.sas").read_text()
    cols = {tag: vars_.split() for tag, vars_ in re.findall(r"^%csv\(\w+, (\w+), ([\w ]+?)(?:,|\))", sas, re.M)}
    rows = {tag: [] for tag in cols}
    for line in Path(log).read_text().splitlines():
        if line.startswith("CSV|"):
            parts = line.rstrip().split("|")
            rows[parts[1]].append(parts[2:])
    (HERE / "results").mkdir(exist_ok=True)
    for tag, r in rows.items():
        assert r, f"no {tag} lines in the log"
        pd.DataFrame(r, columns=[c.lower() for c in cols[tag]]).to_csv(HERE / "results" / f"{tag}.csv", index=False)
    print(", ".join(f"{t} {len(r)}" for t, r in rows.items()))


if __name__ == "__main__":
    main(sys.argv[1])
