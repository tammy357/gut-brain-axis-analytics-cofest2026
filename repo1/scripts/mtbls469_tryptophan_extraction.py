# MTBLS469 — tryptophan extraction from plasma mzXML files
#
# TMAO cannot be measured in this study: its [M+H]+ sits at m/z 76.08, and
# this study's instrument scanned 100-1000, so that mass was never recorded.
# Tryptophan is in range and is one of the metabolites named in Soto et al.
#
# The study's MAF is internally inconsistent about which m/z to expect --
# it lists mass_to_charge as 204.0898 (the neutral monoisotopic mass) while
# labelling the modification [M+H]+ (which would be observed at 205.0972).
# This searches both and reports whichever has signal. The MAF also states a
# retention time of 4.9 min, which is an independent check: a real tryptophan
# peak should appear near there.

import re
from pathlib import Path

import numpy as np
import pandas as pd
from pyteomics import mzxml

CANDIDATES = {
    "tryptophan_MH_205.0972": 205.0972,   # [M+H]+, calculated
    "tryptophan_MAF_204.0898": 204.0898,  # value listed in the study's MAF
}
TOLERANCE_PPM = 10
EXPECTED_RT_MIN = 4.9

FNAME_RE = re.compile(r"^(PR\d+)_v(\d+)_(male|female)_(arm\d)_(baseline|juice|apple)\.mzXML$")


def extract_intensity(path, target_mz, tolerance_ppm=TOLERANCE_PPM):
    delta = target_mz * tolerance_ppm / 1e6
    low, high = target_mz - delta, target_mz + delta
    best_intensity, best_rt, n_scans = 0.0, None, 0

    with mzxml.MzXML(str(path)) as reader:
        for scan in reader:
            if scan.get("msLevel") != 1:
                continue
            mz_array = scan.get("m/z array")
            intensity_array = scan.get("intensity array")
            if mz_array is None or len(mz_array) == 0:
                continue
            mask = (mz_array >= low) & (mz_array <= high)
            if not mask.any():
                continue
            n_scans += 1
            scan_max = intensity_array[mask].max()
            if scan_max > best_intensity:
                best_intensity, best_rt = scan_max, scan.get("retentionTime")

    return best_intensity, best_rt, n_scans


study_dir = Path("mtbls_data/MTBLS469")
all_mzxml = sorted(study_dir.rglob("*.mzXML"))
plasma_files = [p for p in all_mzxml if FNAME_RE.match(p.name)]
print(f"Total mzXML files: {len(all_mzxml)}")
print(f"Plasma files (PR prefix): {len(plasma_files)}")
print()

rows = []
for path in plasma_files:
    subject, visit, gender, arm, meal = FNAME_RE.match(path.name).groups()
    row = {
        "file": path.name, "subject": subject, "visit": int(visit),
        "gender": gender, "arm": arm, "meal": meal,
    }
    for label, mz in CANDIDATES.items():
        intensity, rt, n_scans = extract_intensity(path, mz)
        row[f"{label}_intensity"] = intensity
        row[f"{label}_rt"] = rt
        row[f"{label}_n_scans"] = n_scans
    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv("mtbls469_tryptophan_extraction.csv", index=False)
print(f"Saved mtbls469_tryptophan_extraction.csv ({len(df)} rows)")
print()

for label in CANDIDATES:
    col = f"{label}_intensity"
    n_signal = (df[col] > 0).sum()
    print(f"{label}: {n_signal}/{len(df)} files with signal")
    if n_signal:
        print(f"  median intensity: {df[df[col] > 0][col].median():.0f}")
        rts = df[df[col] > 0][f'{label}_rt'].dropna()
        if len(rts):
            print(f"  median RT: {rts.median():.2f} min (study states {EXPECTED_RT_MIN} min)")
    print()

print("By intervention:")
for label in CANDIDATES:
    col = f"{label}_intensity"
    if (df[col] > 0).any():
        print(f"\n{label}:")
        print(df.groupby("meal")[col].agg(["count", "median"]))
