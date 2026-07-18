"""
Targeted extraction of the trimethylamine-N-oxide (TMAO) peak from raw
mzXML files in a MetaboLights study, as an alternative to a full
untargeted metabolomics processing pipeline.

TMAO (C3H9NO) has a monoisotopic neutral mass of 75.0684 Da. In
positive-mode LC-MS, the expected ion is [M+H]+ at m/z 76.0757. This
script scans each sample's mzXML file, extracts the intensity of the
closest peak within a tolerance window around that m/z in every MS1
scan, and reports the maximum intensity found (a simple proxy for peak
height -- not a proper chromatographic peak integration, which would
need alignment across samples and a real peak-picking algorithm).

This is a bounded, single-compound extraction, not a substitute for a
full pipeline (XCMS, MZmine) if the goal is untargeted discovery -- but
it's enough to check whether TMAO is present and get a relative
intensity comparison across samples for one specific compound.

Requires: pip install pyteomics numpy pandas
"""

import numpy as np
import pandas as pd
from pyteomics import mzxml

TMAO_MZ = 76.0757
DEFAULT_TOLERANCE_PPM = 10


def mz_tolerance_window(target_mz, tolerance_ppm):
    delta = target_mz * tolerance_ppm / 1e6
    return target_mz - delta, target_mz + delta


def extract_target_intensity(mzxml_path, target_mz=TMAO_MZ,
                               tolerance_ppm=DEFAULT_TOLERANCE_PPM):
    """Return the maximum intensity observed within the tolerance window
    around target_mz, across all MS1 scans in the file, along with the
    retention time at which that maximum occurred."""
    low, high = mz_tolerance_window(target_mz, tolerance_ppm)

    best_intensity = 0.0
    best_rt = None
    n_scans_with_signal = 0

    with mzxml.MzXML(str(mzxml_path)) as reader:
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

            n_scans_with_signal += 1
            scan_max = intensity_array[mask].max()
            if scan_max > best_intensity:
                best_intensity = scan_max
                best_rt = scan.get("retentionTime")

    return {
        "max_intensity": best_intensity,
        "retention_time_min": best_rt,
        "n_scans_with_signal": n_scans_with_signal,
    }


def process_sample_list(mzxml_paths, target_mz=TMAO_MZ,
                          tolerance_ppm=DEFAULT_TOLERANCE_PPM):
    """Run extract_target_intensity over a list of mzXML file paths and
    return a summary DataFrame, one row per file."""
    rows = []
    for path in mzxml_paths:
        result = extract_target_intensity(path, target_mz, tolerance_ppm)
        result["file"] = str(path)
        rows.append(result)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python targeted_tmao_extraction.py <mzxml_file_1> [<mzxml_file_2> ...]")
        sys.exit(1)

    df = process_sample_list(sys.argv[1:])
    print(df.to_string(index=False))
    df.to_csv("tmao_extraction_results.csv", index=False)
    print("\nSaved tmao_extraction_results.csv")
