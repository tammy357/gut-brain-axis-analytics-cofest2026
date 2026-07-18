# Cell — extended targeted extraction: TMAO + 4 more paper-seed metabolites
#
# The MTBLS218 MAF file only lists 19 manually-annotated compounds, and none
# of hydroxyproline, guanidinoacetic acid, SDMA, or the carnitines are among
# them. That reflects what the original authors chose to annotate, not what
# the raw spectra contain. This searches the raw mzXML files directly for
# each compound's expected mass, independent of the study's own annotation.
#
# m/z values calculated from monoisotopic mass + proton mass for [M+H]+,
# except carnitine derivatives, which are inherently charged (trimethyl-
# ammonium cation) and observed directly as [M]+, not [M+H]+.

import numpy as np
import pandas as pd
from pyteomics import mzxml
from pathlib import Path

TARGETS = {
    "TMAO":                  76.0757,
    "Hydroxyproline":        132.0655,
    "Guanidinoacetic_acid":  118.0611,
    "SDMA":                  203.1503,
    # Carnitine cation backbone check (unesterified) -- a sanity-check mass,
    # NOT the same as C14/C5-DC/C10:2 carnitine esters, which are heavier.
    # Included here as a lower-confidence exploratory check only.
    "Carnitine_free":        162.1125,
}

DEFAULT_TOLERANCE_PPM = 10


def mz_tolerance_window(target_mz, tolerance_ppm):
    delta = target_mz * tolerance_ppm / 1e6
    return target_mz - delta, target_mz + delta


def extract_target_intensity(mzxml_path, target_mz, tolerance_ppm=DEFAULT_TOLERANCE_PPM):
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

    return best_intensity, best_rt, n_scans_with_signal


def process_all_targets(mzxml_paths, targets=TARGETS):
    rows = []
    for path in mzxml_paths:
        row = {"file": path.name}
        for compound_name, target_mz in targets.items():
            intensity, rt, n_scans = extract_target_intensity(path, target_mz)
            row[f"{compound_name}_intensity"] = intensity
            row[f"{compound_name}_rt_min"] = rt
            row[f"{compound_name}_n_scans"] = n_scans
        rows.append(row)
        print(f"Processed {path.name}")
    return pd.DataFrame(rows)


if __name__ == "__main__":
    study_dir = Path("mtbls_data/MTBLS218/FILES")
    mzxml_files = sorted(study_dir.glob("*.mzXML"))
    print(f"Found {len(mzxml_files)} mzXML files")

    df = process_all_targets(mzxml_files)
    df.to_csv("multi_metabolite_extraction_results.csv", index=False)
    print("Saved multi_metabolite_extraction_results.csv")
