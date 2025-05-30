"""
Protection Factor to Deuterium Uptake Prediction
------------------------------------------------

This script calculates the predicted deuterium uptake (Dpred) using protection factors (Pfactors)
and intrinsic exchange rates (kint). Supports multiple replicates and noise injection.

Copyright (C) 2019-2020 Emanuele Paci group.
Code upgraded by Mahmoud Shalash
Licensed under GPL-2.0
"""

import os
import argparse
from typing import Dict, List, Optional

from calc_dpred import calculate_dpred
from kint import calculate_kint_for_sequence
from read import read_assignments, read_pfact, read_seq, read_time_points
from write import write_dpred, write_combined_replicates
from logger import log


def run_pfact2dpred(
    base: str = ".",
    ass_file: str = "",
    pfact_file: str = "",
    times_file: str = "",
    temp: float = 298.15,
    pH: float = 7.0,
    seq_file: str = "",
    out: Optional[str] = None,
    nrep: int = 1,
    eps: float = 0.0
) -> None:
    """
    Run the protection factor to deuterium uptake prediction process.

    Parameters
    ----------
    base : str
        Base directory where files are located.
    ass_file : str
        Assignment file path.
    pfact_file : str
        Path to the protection factor file.
    times_file : str
        Path to file with labeling times.
    temp : float
        Temperature in Kelvin.
    pH : float
        pH of the solution.
    seq_file : str
        File with amino acid sequence.
    out : str
        Output prefix (without extension).
    nrep : int
        Number of replicates to generate.
    eps : float
        Standard deviation of Gaussian noise to add (0 disables).
    """
    log.info("Running pfact2dpred")

    # Read input files
    pfact = read_pfact(pfact_file)
    assignments = read_assignments(ass_file)
    sequence = read_seq(seq_file)
    times = read_time_points(times_file)
    res1, resn = 1, len(sequence)

    kint, prolines = calculate_kint_for_sequence(
        res1, resn, sequence, temp, pH
    )

    outfiles: List[str] = []
    for i in range(1, nrep + 1):
        outfile = f"{out}{i}" if nrep > 1 else out
        dpred = calculate_dpred(pfact, times, kint, assignments)
        write_dpred(outfile, dpred, times, eps)
        outfiles.append(f"{outfile}.Dpred")

    if nrep > 1:
        write_combined_replicates(outfiles, out=out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Predict D-uptake from protection factors and assignments."
    )

    parser.add_argument("--base", type=str, default=".")
    parser.add_argument("--ass", required=True)
    parser.add_argument("--pfact", required=True)
    parser.add_argument("--times", required=True)
    parser.add_argument("--temp", required=True, type=float)
    parser.add_argument("--pH", required=True, type=float)
    parser.add_argument("--seq", required=True)
    parser.add_argument("--out", required=False, default="dpred")
    parser.add_argument("--nrep", type=int, default=1)
    parser.add_argument("--eps", type=float, default=0.0)

    args = parser.parse_args()

    run_pfact2dpred(
        base=args.base,
        ass_file=args.ass,
        pfact_file=args.pfact,
        times_file=args.times,
        temp=args.temp,
        pH=args.pH,
        seq_file=args.seq,
        out=args.out,
        nrep=args.nrep,
        eps=args.eps
    )
