"""
Isotopic Envelope Analysis Tool (ExPfact Module)
================================================

Modes:
- Prediction mode (`p`): Generate theoretical isotopic envelopes from protection factors
- Comparison mode (`c`): Apply back-exchange correction and compare experimental data

Usage:
    Mode `p`:
        python isotopic_envelope.py --mode p --ass <assignments> --seq <sequence_file>
               --T_label <temp> --pH_label <pH> --pfact <lnP_file> --times <times_file>
               --pep <peptide_index> --z <charge_state> --prefix <output_prefix>

    Mode `c`:
        python isotopic_envelope.py --mode c --ass <assignments> --seq <sequence_file>
               --T_label <label_temp> --pH_label <label_pH> --T_quench <quench_temp>
               --pH_quench <quench_pH> --pfact <lnP_file> --times <times_file>
               --pep <peptide_index> --z <charge_state> --prefix <exp_prefix>

Author: E. Paci group (GPL-2.0), Upgraded by Mahmoud Shalash
"""

import os
import argparse
import pandas as pd

from logger import log
from read import read_time_points
from isenv_functions import (
    predict_isotopic_envelope,
    generate_back_exchange_time_points,
    sticks_from_experimental_envelope,
    compare_predictions,
)


def corrected_isotopic_envelope_prediction(
    ass_file: str, seq_file: str, T_label: float, pH_label: float,
    T_quench: float, pH_quench: float, lnP_file: str, times_file: str,
    exp_env_pre: str, pep: int, z: int
) -> None:
    """
    Apply forward and back exchange correction to compare with experimental data.
    Saves corrected isotopic envelopes and R² comparisons.
    """
    predict_isotopic_envelope(
        ass_file, seq_file, T_label, pH_label, lnP_file, times_file,
        pep, z, exchange='f', out_file=f"{exp_env_pre}/{pep}"
    )

    generate_back_exchange_time_points()
    times = read_time_points(times_file)

    for exp_idx in range(1, len(times)):
        outfile = f"{exp_env_pre}/{pep}.{exp_idx}.corr"
        pi0file = f"{exp_env_pre}/{pep}.{exp_idx}.isot"

        predict_isotopic_envelope(
            ass_file, seq_file, T_quench, pH_quench, lnP_file,
            "back.times", pep, z, exchange='b', out_file=outfile,
            pi0_file=pi0file
        )

        corr_env = pd.read_csv(f"{outfile}.1.isot", delim_whitespace=True, skiprows=1, header=None)
        exp_env = pd.read_csv(f"{exp_env_pre}.{exp_idx}.txt", delim_whitespace=True, header=None)

        mass, fr = sticks_from_experimental_envelope(exp_env, corr_env, z)

        result_df = compare_predictions(fr, prefix=f"{outfile}*")
        result_df.to_csv(f"{exp_env_pre}/tau.{pep}.{exp_idx}.res", sep='\t', index=False)


def main() -> None:
    """
    Entry point for isotopic_envelope.py.
    Handles both prediction and comparison modes.
    """
    parser = argparse.ArgumentParser(description="Predict or compare HDX isotopic envelopes.")
    parser.add_argument("--mode", choices=["p", "c"], required=True, help="Mode: prediction (p) or comparison (c)")
    parser.add_argument("--ass", required=True, help="Assignment file")
    parser.add_argument("--seq", required=True, help="Protein sequence file")
    parser.add_argument("--T_label", required=True, type=float, help="Labeling temperature")
    parser.add_argument("--pH_label", required=True, type=float, help="Labeling pH")
    parser.add_argument("--pfact", required=True, help="lnP protection factor file")
    parser.add_argument("--times", required=True, help="Exposure times file")
    parser.add_argument("--pep", required=True, type=int, help="Peptide index")
    parser.add_argument("--z", required=True, type=int, help="Peptide charge state")
    parser.add_argument("--prefix", required=True, help="Output prefix or experimental envelope folder")

    # Optional (only required for comparison mode)
    parser.add_argument("--T_quench", type=float, help="Quenching temperature (comparison mode only)")
    parser.add_argument("--pH_quench", type=float, help="Quenching pH (comparison mode only)")

    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.prefix, exist_ok=True)

    if args.mode == 'p':
        log.info("Running isotopic_envelope.py in prediction mode")
        predict_isotopic_envelope(
            args.ass, args.seq, args.T_label, args.pH_label,
            args.pfact, args.times, args.pep, args.z,
            exchange='f', out_file=f"{args.prefix}/{args.pep}"
        )

    elif args.mode == 'c':
        log.info("Running isotopic_envelope.py in comparison mode")
        if args.T_quench is None or args.pH_quench is None:
            raise ValueError("Comparison mode requires --T_quench and --pH_quench.")
        corrected_isotopic_envelope_prediction(
            args.ass, args.seq, args.T_label, args.pH_label,
            args.T_quench, args.pH_quench, args.pfact, args.times,
            args.prefix, args.pep, args.z
        )


if __name__ == "__main__":
    main()
