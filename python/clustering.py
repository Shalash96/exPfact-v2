"""
Clustering of HDX-MS contiguous peptide regions using Mclust (via Rscript).

- Reads peptide assignments from a file
- Identifies contiguous regions between gaps in coverage
- Calls `multi.r` R script to cluster each region using Mclust

Author: Code updated by Mahmoud Shalash based on original code by E. Paci group
"""

import argparse
import subprocess
import numpy as np
from read import read_assignments
from typing import List


def find_covered_regions_between_gaps(assignments: np.ndarray) -> List[str]:
    """
    Identifies contiguous covered regions based on gaps between uncovered residues,
    skipping the first amino acid of each peptide (due to rapid back-exchange in HDX-MS).

    Parameters
    ----------
    assignments : np.ndarray
        Peptide assignment array with shape (N, 3), where columns are:
        [peptide_index, start_residue (1-based), end_residue (1-based)]

    Returns
    -------
    List[str]
        List of contiguous covered regions (as 'start-end' strings), inferred between uncovered blocks.
    """
    max_residue = int(np.max(assignments[:, 2]))
    coverage = np.zeros(max_residue, dtype=bool)

    starts = assignments[:, 1].astype(int) - 1  # convert to 0-based
    ends = assignments[:, 2].astype(int)

    for start, end in zip(starts, ends):
        if end > start + 1:
            coverage[start + 1:end] = True  # Skip first residue of peptide

    uncovered_indices = np.where(~coverage)[0]

    if len(uncovered_indices) == 0:
        return [f"1-{max_residue}"]

    uncovered_residues = uncovered_indices + 1
    regions = []

    if uncovered_residues[0] > 1:
        regions.append(f"1-{uncovered_residues[0] - 1}")

    if len(uncovered_residues) > 1:
        gaps = np.diff(uncovered_residues)
        gap_indices = np.where(gaps > 1)[0]
        for i in gap_indices:
            start = uncovered_residues[i] + 1
            end = uncovered_residues[i + 1] - 1
            regions.append(f"{start}-{end}")

    if uncovered_residues[-1] < max_residue:
        regions.append(f"{uncovered_residues[-1] + 1}-{max_residue}")

    return regions



def run_mclust_on_regions(regions: list[str], r_script: str = "../R/multi.r", sp_path: str = "all.sp") -> None:
    """
    Run Mclust clustering via Rscript for each contiguous covered region.

    Parameters
    ----------
    regions : List[str]
        List of peptides ranges in format 'start-end'.
    r_script : str
        Path to the R clustering script.
    sp_path : str
        Path to all.sp file
    """
    for region in regions:
        start, end = region.split('-')
        command = ["Rscript", r_script, start, end, sp_path]
        try:
            subprocess.run(command, check=True)
            print(f"Clustered region: {start}-{end}")
        except subprocess.CalledProcessError:
            print(f"Failed clustering for region: {start}-{end}")


def main():
    parser = argparse.ArgumentParser(description="Run Mclust clustering on contiguous regions between gaps.")
    parser.add_argument("--ass", required=True, help="Path to peptide assignment file (.list)")
    parser.add_argument("--all_sp", default="all.sp", help="Path to sp file (default: all.sp)")
    args = parser.parse_args()

    assignments = read_assignments(args.ass)
    regions = find_covered_regions_between_gaps(assignments)
    run_mclust_on_regions(regions, sp_path=args.all_sp)


if __name__ == "__main__":
    main()
