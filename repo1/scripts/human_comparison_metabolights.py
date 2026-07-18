"""
Pull human TMAO data from MetaboLights studies MTBLS310 (red wine trial)
and MTBLS218 (bariatric surgery), and compare against the mouse TMAO
finding from the ST000885 reanalysis.

Run this in an environment with network access to ebi.ac.uk -- it
times out from a sandboxed environment with no FTP egress.

Requires: pip install metabolights-utils pandas
"""

import json
import subprocess
from pathlib import Path

import pandas as pd

STUDIES = {
    "MTBLS310": "Red wine intervention, coronary artery disease patients",
    "MTBLS218": "Bariatric surgery, metabolic status before/after",
}

DOWNLOAD_ROOT = Path("mtbls_data")


def download_study(study_id):
    """Download ISA metadata + data files for a public MetaboLights study."""
    DOWNLOAD_ROOT.mkdir(exist_ok=True)
    result = subprocess.run(
        ["mtbls", "public", "download", study_id, "-p", str(DOWNLOAD_ROOT)],
        capture_output=True, text=True, timeout=300,
    )
    print(result.stdout)
    if result.returncode != 0:
        print(f"Download failed for {study_id}: {result.stderr}")
    return result.returncode == 0


def describe_study(study_id, jsonpath=None):
    """Print or filter the study's model via the mtbls describe command."""
    cmd = ["mtbls", "public", "describe", study_id]
    if jsonpath:
        cmd.append(jsonpath)
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    print(result.stdout)
    return result.stdout


def list_all_downloaded_files(study_id):
    """List every file actually present in the downloaded study folder,
    so a naming assumption that doesn't match reality is visible
    immediately rather than silently producing zero results."""
    study_dir = DOWNLOAD_ROOT / study_id
    if not study_dir.exists():
        print(f"  {study_dir} does not exist")
        return []
    files = sorted(study_dir.rglob("*"))
    files = [f for f in files if f.is_file()]
    for f in files:
        print(f"  {f.relative_to(study_dir)}")
    return files


def find_maf_files(study_id):
    """Metabolite Annotation Format files (m_*.txt or m_*.tsv) hold the
    actual metabolite abundance tables once a study is downloaded."""
    study_dir = DOWNLOAD_ROOT / study_id
    if not study_dir.exists():
        return []
    return list(study_dir.rglob("m_*.txt")) + list(study_dir.rglob("m_*.tsv"))


def search_metabolite_in_maf(maf_path, metabolite_names):
    """Look for rows matching any of the given metabolite names in a MAF
    file. MAF files are tab-separated with a 'metabolite_identification'
    or similar column holding the compound name."""
    df = pd.read_csv(maf_path, sep="\t")
    name_col_candidates = [c for c in df.columns
                            if "identification" in c.lower() or "metabolite" in c.lower()]
    if not name_col_candidates:
        return pd.DataFrame()
    name_col = name_col_candidates[0]

    mask = df[name_col].astype(str).str.contains(
        "|".join(metabolite_names), case=False, na=False
    )
    return df[mask]


def main():
    target_metabolites = ["trimethylamin", "TMAO", "hydroxyprolin"]

    for study_id, description in STUDIES.items():
        print(f"\n=== {study_id}: {description} ===")
        download_study(study_id)
        describe_study(study_id, "$.investigation.studies[0].title")

        maf_files = find_maf_files(study_id)
        print(f"Found {len(maf_files)} MAF file(s)")

        if not maf_files:
            print("  No MAF files matched -- listing everything actually downloaded:")
            list_all_downloaded_files(study_id)

        for maf_path in maf_files:
            hits = search_metabolite_in_maf(maf_path, target_metabolites)
            if not hits.empty:
                print(f"  {maf_path.name}: {len(hits)} matching row(s)")
                print(hits.to_string(index=False))
            else:
                print(f"  {maf_path.name}: no matches for {target_metabolites}")


if __name__ == "__main__":
    main()
