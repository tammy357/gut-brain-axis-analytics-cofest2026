# Cell 1 — install everything needed
!pip install -q metabolights-utils pyteomics numpy pandas


# Cell 2 — download the raw mzXML files specifically from the FILES subfolder
# (the ISA-Tab metadata alone, which is what earlier attempts pulled, does
# NOT include the raw spectral data -- that lives under a separate FILES path)
import subprocess

result = subprocess.run(
    ["mtbls", "public", "download", "MTBLS218", "FILES", "-p", "mtbls_data", "-o"],
    capture_output=True, text=True, timeout=600,
)
print("STDOUT:")
print(result.stdout)
print("STDERR:")
print(result.stderr)
print("Return code:", result.returncode)


# Cell 3 — confirm what actually landed on disk
from pathlib import Path

study_dir = Path("mtbls_data/MTBLS218")
if study_dir.exists():
    all_files = sorted(study_dir.rglob("*"))
    mzxml_files = [f for f in all_files if f.suffix.lower() == ".mzxml"]
    print(f"Total files found under {study_dir}: {len([f for f in all_files if f.is_file()])}")
    print(f"mzXML files found: {len(mzxml_files)}")
    for f in mzxml_files[:10]:
        print(" ", f)
else:
    print(f"{study_dir} does not exist -- download did not complete")


# Cell 4 — run the targeted TMAO extraction across whatever mzXML files were found
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


if mzxml_files:
    rows = []
    for path in mzxml_files:
        result = extract_target_intensity(path)
        result["file"] = path.name
        rows.append(result)
        print(f"Processed {path.name}: max_intensity={result['max_intensity']:.0f}")

    df = pd.DataFrame(rows)
    df.to_csv("tmao_extraction_results.csv", index=False)
    print("\nFull results:")
    print(df.to_string(index=False))
    print("\nSaved tmao_extraction_results.csv")
else:
    print("No mzXML files were found -- check Cell 2's output above for download errors")
