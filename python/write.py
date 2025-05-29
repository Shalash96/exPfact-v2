"""
Write Utility Functions for ExPfact Outputs
-------------------------------------------

Handles output of protection factors, predicted deuterium uptake,
cost differences, and pooled replicates for HDX-MS data analysis.

Copyright (C) 2019-2020 Emanuele Paci, Simon P. Skinner, Michele Stofella
Upgraded by Mahmoud Shalash
Licensed under GPL-2.0
"""

import numpy as np
from typing import List


def write_pfact(params: np.ndarray, fout_name: str) -> None:
    """
    Write protection factor array to a .pfact file.

    Parameters
    ----------
    params : np.ndarray
        Array of ln(P) values.
    fout_name : str
        Output filename prefix (no extension).
    """
    with open(fout_name + '.pfact', 'w') as f:
        for i, value in enumerate(params):
            f.write(f"{i + 1} {value}\n")


def write_dpred(
    output_file: str,
    dpred: np.ndarray,
    times: List[float],
    eps: float = 0.0,
    suffix: str = ".Dpred"
) -> None:
    """
    Write predicted deuterium uptake to file.

    Parameters
    ----------
    output_file : str
        Output filename prefix.
    dpred : np.ndarray
        Predicted uptake matrix (peptides × timepoints).
    times : List[float]
        List of timepoints.
    eps : float, optional
        Gaussian noise std for simulation. Default is 0.
    suffix : str, optional
        File suffix. Default is '.Dpred'.
    """
    output_array = np.insert(dpred, 0, times, axis=0)
    if eps > 0:
        noise = np.random.normal(scale=eps, size=output_array[1:].shape)
        output_array[1:] += noise
        np.clip(output_array[1:], 0, 1, out=output_array[1:])  # Ensure values ∈ [0, 1]
    np.savetxt(output_file + suffix, output_array.T, fmt='%.7g')


def write_diff(outfile: str, dpred: np.ndarray, dexp: np.ndarray) -> None:
    """
    Write per-peptide RMS difference between prediction and experiment.

    Parameters
    ----------
    outfile : str
        Output file prefix.
    dpred : np.ndarray
        Predicted uptake array.
    dexp : np.ndarray
        Experimental uptake array.
    """
    with open(outfile + '.diff', 'w') as f:
        for i, (pred, exp) in enumerate(zip(dpred, dexp)):
            cost = np.mean((pred - exp) ** 2)
            f.write(f"{i + 1} {cost:.8e}\n")


def write_combined_replicates(files: List[str], out: str) -> None:
    """
    Combine multiple .Dpred replicate files by averaging.

    Parameters
    ----------
    files : List[str]
        List of .Dpred filenames.
    out : str
        Output filename prefix.
    """
    dpred_arrays = [np.loadtxt(f) for f in files]
    comb = np.mean(dpred_arrays, axis=0)

    # Calculate pooled standard deviation across replicates
    all_weights = np.array([arr[:, 1:] for arr in dpred_arrays])  # shape: (replicates, peptides, timepoints)
    pooled_std = np.sqrt(np.mean(np.var(all_weights, axis=0)))
    print(f"Pooled std: {pooled_std:.5f}")

    np.savetxt(out + '.Dpred', comb, fmt='%.7g')
