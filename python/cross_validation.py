"""
Cross-validation module for exPfact.

Performs leave-one-out cross-validation (LOOCV) across a range of harmonic penalty values (lambda)
to determine model generalizability for HDX-MS protection factor fitting.

Author: E. Paci group (GPL-2.0), Refactored for clarity and performance
"""

from exPfact import run
from read import read_assignments, read_seq, read_dexp, read_pfact, read_configuration
from calc_dpred import calculate_dpred
from kint import calculate_kint_for_sequence
from logger import log

import numpy as np
import os
import argparse
import sys
from typing import Tuple


# Define lambda range: 15 values log-spaced from 1e-15 to 1e-1
lambdas: np.ndarray = np.logspace(-15, -1, 15, endpoint=True)


def loo_dataset(dexp: np.ndarray, time_points: np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Split experimental data by leaving out time point k."""
    dexp_train = np.delete(dexp, k, axis=1)
    times_train = np.delete(time_points, k)
    dexp_test = dexp[:, k]
    times_test = np.array([time_points[k]])
    return dexp_train, times_train, dexp_test, times_test


def loo_crossval(
    dexp: np.ndarray,
    time_points: np.ndarray,
    ass: np.ndarray,
    lam: float,
    pH: float,
    temp: float,
    seq: str,
    res1: int,
    resn: int
) -> Tuple[float, float]:
    """Perform leave-one-out cross-validation for a single lambda."""
    cv_train = 0.0
    cv_test = 0.0
    kint, _ = calculate_kint_for_sequence(res1, resn, seq, temp, pH)

    for k in range(len(time_points)):
        out_file = f"CVout.rm{k}"
        dexp_train, times_train, dexp_test, times_test = loo_dataset(dexp, time_points, k)

        run(
            base_dir=os.getcwd(),
            dexp=dexp_train,
            assignments=ass,
            pfact=None,
            random_steps=None,
            time_points=times_train,
            harmonic_term=lam,
            output_file=out_file,
            tolerance=1e-10,
            weights=None,
            pH=pH,
            temperature=temp,
            seq=seq,
            res1=res1,
            resn=resn,
        )

        pfact = read_pfact(out_file + ".pfact")
        dpred_test = calculate_dpred(pfact, times_test, kint, ass)
        cost_test = np.mean([(pred - exp) ** 2 for pred, exp in zip(dpred_test, dexp_test)])

        diff = np.loadtxt(out_file + ".diff")
        train_error = diff[1] if diff.ndim == 1 else np.sum(diff[:, 1])

        cv_train += train_error
        cv_test += cost_test

        for ext in [".Dpred", ".diff", ".pfact"]:
            try:
                os.remove(out_file + ext)
            except FileNotFoundError:
                pass

    n = len(time_points)
    return cv_train / n, cv_test / n


def cross_validate(
    dexp: np.ndarray,
    time_points: np.ndarray,
    ass: np.ndarray,
    lambdas: np.ndarray,
    pH: float,
    temp: float,
    seq: str,
    res1: int,
    resn: int
) -> None:
    """Run LOOCV for each lambda value and write results."""
    with open("CVtest2.res", "w") as fout:
        for lam in lambdas:
            log.info("Calculating at lambda = %.2e" % lam)
            CVtrain, CVtest = loo_crossval(dexp, time_points, ass, lam, pH, temp, seq, res1, resn)
            fout.write(f"{lam:.5e} {CVtrain:.6f} {CVtest:.6f}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Perform LOOCV for exPfact.")
    parser.add_argument("--dexp", help="Path to experimental D-uptake file")
    parser.add_argument("--ass", help="Path to assignment (.list) file")
    parser.add_argument("--temp", help="Temperature in Kelvin", type=float)
    parser.add_argument("--pH", help="Solution pH", type=float)
    parser.add_argument("--seq", help="Path to sequence file")

    dexp: np.ndarray
    time_points: np.ndarray
    ass: np.ndarray
    temp: float
    pH: float
    seq: str
    res1: int
    resn: int

    if len(sys.argv) > 1 and sys.argv[1].endswith(".json"):
        config = read_configuration(sys.argv[1])
        dexp = config["dexp"]
        time_points = config["times"]
        ass = config["assignments"]
        temp = config["temperature"]
        pH = config["pH"]
        seq = read_seq(os.path.join(config["base"], config["seq_file"]))
        res1 = 1
        resn = len(seq)
    else:
        opts = parser.parse_args()

        # Validate all required CLI arguments
        if not all([opts.dexp, opts.ass, opts.temp, opts.pH, opts.seq]):
            log.error("Missing required arguments: --dexp, --ass, --temp, --pH, --seq")
            sys.exit(1)

        dexp, time_points = read_dexp(opts.dexp)
        ass = read_assignments(opts.ass)
        temp = float(opts.temp)
        pH = float(opts.pH)
        seq = read_seq(opts.seq)
        res1 = 1
        resn = len(seq)

    log.info("Running cross-validation on input data.")
    cross_validate(dexp, time_points, ass, lambdas, pH, temp, seq, res1, resn)