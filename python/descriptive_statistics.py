"""
Descriptive statistics tool for protection factor fitting results in exPfact.

This script:
- Selects the top X% best-scoring solutions (based on .diff files)
- Aggregates their ln(P) profiles
- Computes average, median, min/max protection factors across all residues

Author: E. Paci group (GPL-2.0), Enhanced by Mahmoud Shalash
"""

import os
import numpy as np
import pandas as pd
from read import read_pfact
from logger import log
from typing import List
import argparse

def select_top_solutions(out_file: str, top_percent: float) -> None:
    """
    Selects top scoring protection factor (.pfact) files based on average error in .diff files.

    Parameters
    ----------
    out_file : str
        Prefix of result files to consider (e.g., 'out' for out1.diff, out2.diff, etc.)
    top_percent : float
        Percentage (0-100) of best solutions to keep
    """
    diff_entries = []

    for file in os.listdir():
        if file.endswith(".diff") and out_file in file:
            try:
                diff_data = pd.read_csv(file, header=None, sep="\\s+")
                avg_diff = np.mean(diff_data[1])
                diff_entries.append((file, avg_diff))
            except Exception as e:
                log.warning(f"Skipping file {file} due to error: {e}")

    if not diff_entries:
        log.error("No matching .diff files found.")
        return

    # Sort files by ascending average cost
    sorted_entries = sorted(diff_entries, key=lambda x: x[1])
    n_top = max(1, int(len(sorted_entries) * top_percent / 100))

    with open("diff.list", "w") as f:
        for file, val in sorted_entries:
            f.write(f"{file} {val:.5f}\n")
    log.info(f"Top {n_top} solutions selected from {len(sorted_entries)} total.")

    # Collect protection factors from selected files
    selected_pfactors: List[List[float]] = []
    for file, _ in sorted_entries[:n_top]:
        pfact_file = file.replace(".diff", ".pfact")
        pfactors = read_pfact(pfact_file)
        selected_pfactors.append(pfactors)

    with open("all.sp", "w") as f:
        for row in selected_pfactors:
            f.write(" ".join(f"{val:.5f}" for val in row) + "\n")

    log.info("Top solution protection factors written to all.sp")

def run_descriptive() -> None:
    """
    Computes descriptive statistics (mean, std, median, min, max) across top solutions in all.sp
    and saves the results to respective .pfact summary files.
    """
    allsp = pd.read_csv("all.sp", header=None, sep="\\s+") # \\s+ handles any whitespace. I used it instead of delim_whitespace=True because it will be deprecated in next versions of pandas

    with open("average.pfact", "w") as f:
        for i in range(allsp.shape[1]):
            mean = np.mean(allsp[i])
            std = np.std(allsp[i])
            f.write(f"{i+1}\t{mean:.5f}\t{std:.5f}\n")
    log.info("Average and std dev saved in average.pfact")

    with open("median.pfact", "w") as f:
        for i in range(allsp.shape[1]):
            median = np.median(allsp[i])
            f.write(f"{i+1}\t{median:.5f}\n")
    log.info("Median protection factors saved in median.pfact")

    with open("minmax.pfact", "w") as f:
        for i in range(allsp.shape[1]):
            min_val = np.min(allsp[i])
            max_val = np.max(allsp[i])
            f.write(f"{i+1}\t{min_val:.5f}\t{max_val:.5f}\n")
    log.info("Min/max protection factors saved in minmax.pfact")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Select and analyze top protection factor solutions.")
    parser.add_argument("--res", required=True, help="Prefix of result files to analyze (e.g., 'out')")
    parser.add_argument("--top", type=float, default=50.0, help="Top %% of solutions to include [default: 50]")

    opts = parser.parse_args()

    log.info("Running descriptive_statistics.py")
    select_top_solutions(opts.res, opts.top)
    run_descriptive()
