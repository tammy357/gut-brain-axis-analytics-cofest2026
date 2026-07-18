# Cell 1 — download the mzXML files (4.46GB study; this will take a while)
!mtbls public download MTBLS469 FILES -p mtbls_data -o


# Cell 2 — check what landed
from pathlib import Path

study_dir = Path("mtbls_data/MTBLS469")
mzxml_files = sorted(study_dir.rglob("*.mzXML"))
print(f"mzXML files found: {len(mzxml_files)}")
for f in mzxml_files[:5]:
    print(" ", f.name)


# Cell 3 — targeted extraction for TMAO across the plasma files
#
# MTBLS469 filenames encode everything we need:
#   PR01_v2_male_arm1_juice.mzXML -> subject PR01, visit 2, arm1, control juice
#   PR02_v2_female_arm1_apple.mzXML -> subject PR02, visit 2, arm1, whole apple
#
# The study's own MAF only annotates apple polyphenols and tryptophan (TMAO was
# never one of its targets), so this searches the spectra directly at TMAO's
# known mass instead of relying on the study's annotation choices.

import re
import numpy as np
import pandas as pd
from pyteomics import mzxml

TMAO_MZ = 76.0757
TOLERANCE_PPM = 10


def extract_intensity(path, target_mz=TMAO_MZ, tolerance_ppm=TOLERANCE_PPM):
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


FNAME_RE = re.compile(r"^(PR\d+)_v(\d+)_(male|female)_(arm\d)_(baseline|juice|apple)\.mzXML$")

rows = []
for path in mzxml_files:
    m = FNAME_RE.match(path.name)
    if not m:
        print(f"Skipping unparsed filename: {path.name}")
        continue

    subject, visit, gender, arm, meal = m.groups()
    intensity, rt, n_scans = extract_intensity(path)
    rows.append({
        "file": path.name,
        "subject": subject,
        "visit": int(visit),
        "gender": gender,
        "arm": arm,
        "meal": meal,
        "tmao_intensity": intensity,
        "tmao_rt_min": rt,
        "tmao_n_scans": n_scans,
    })
    print(f"{path.name}: {intensity:.0f}")

df = pd.DataFrame(rows)
df.to_csv("mtbls469_tmao_extraction.csv", index=False)
print(f"\nSaved mtbls469_tmao_extraction.csv ({len(df)} rows)")
print()
print(df.groupby("meal")["tmao_intensity"].describe())
