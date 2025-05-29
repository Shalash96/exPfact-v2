"""
Isotopic envelope utilities for HDX-MS data analysis.

Provides functions to calculate deuterium uptake, isotopic envelope predictions,
and evaluate back-exchange kinetics from protection factor and kint values.

Original Authors: E. Paci group (GPL-2.0)
Updated by Mahmoud Shalash
"""

import pandas as pd
import numpy as np
import glob
import re
from sklearn.metrics import r2_score
from itertools import combinations

from read import read_assignments, read_seq, read_pfact, read_time_points
from kint import calculate_kint_for_sequence
from kback import calculate_kback_for_sequence
from Hisotope import fully_protonated_envelope


def tryint(s):
    """Try to convert a string to integer, fallback to original string."""
    try:
        return int(s)
    except Exception:
        return s


def alphanum_key(s):
    """Split string into a list of strings and integers for natural sorting."""
    return [tryint(c) for c in re.split(r'([0-9]+)', s)]


def natural_sort(lst):
    """Sort list of filenames or strings in human-readable order."""
    lst.sort(key=alphanum_key)


def deuterium_uptake(t, kint, lnp):
    """
    Calculate total fractional deuterium uptake at time t.
    """
    p = np.exp(lnp)
    total = 0.0
    count = 0
    for i in range(len(kint)):
        if kint[i] >= 0:
            count += 1
            total += np.exp(-kint[i] / p[i] * t * 60)
    return (count - total) / count if count else 0.0


def single_deuterium_uptake(t, kint, lnp):
    """Calculate uptake for a single residue."""
    p = np.exp(lnp)
    return 1 - np.exp(-kint / p * t * 60) if kint > 0 else 0


def single_deuterium_downtake(t, kback, lnp):
    """Calculate back-exchange for a single residue."""
    p = np.exp(lnp)
    return 1 - np.exp(-kback / p * t * 60) if kback > 0 else 0


def set_a_group(a, k):
    """Return all k-combinations of elements in list a."""
    return list(combinations(a, k))


def probability_uptake(k, t, kint, lnp):
    """Probability that exactly k residues exchanged at time t (forward)."""
    a = set_a_group(range(len(kint)), k)
    probability = 0.0
    for element in a:
        d1 = [single_deuterium_uptake(t, kint[j], lnp[j]) for j in element]
        d2 = [1 - single_deuterium_uptake(t, kint[j], lnp[j]) for j in range(len(kint)) if j not in element]
        probability += np.prod(d1) * np.prod(d2)
    return probability


def probability_downtake(k, t, kint, lnp):
    """Probability that exactly k residues exchanged at time t (back-exchange)."""
    a = set_a_group(range(len(kint)), k)
    probability = 0.0
    for element in a:
        d1 = [single_deuterium_downtake(t, kint[j], lnp[j]) for j in element]
        d2 = [1 - single_deuterium_downtake(t, kint[j], lnp[j]) for j in range(len(kint)) if j not in element]
        probability += np.prod(d1) * np.prod(d2)
    return probability


def isotopic_envelope(t, kint, lnp, exchange):
    """Returns envelope probability distribution for either forward or back exchange."""
    envelope = []
    for k in range(2 * len(kint)):
        if exchange == 'f':
            envelope.append(probability_uptake(k, t, kint, lnp))
        elif exchange == 'b':
            envelope.append(probability_downtake(k, t, kint, lnp))
    return envelope


def centered_isotopic_envelope(t, kint, lnp, fr0):
    """Convolves fr0 with forward isotopic exchange at time t."""
    fr = isotopic_envelope(t, kint, lnp, exchange='f')
    f = np.zeros(len(fr0))
    for i in range(len(f)):
        for j in range(i + 1):
            if i - j < len(fr0) and j < len(fr):
                f[i] += fr0[i - j] * fr[j]
    return f


def back_centered_isotopic_envelope(t, kint, lnp, fr0):
    """Convolves fr0 with back-exchange envelope at time t."""
    fr = isotopic_envelope(t, kint, lnp, exchange='b')
    f = np.zeros(len(fr0))
    for i in range(len(f)):
        for j in range(len(f) - i):
            if i + j < len(fr0) and j < len(fr):
                f[i] += fr0[i + j] * fr[j]
    return [val / sum(f) * 100 for val in f] if sum(f) != 0 else f


def predict_isotopic_envelope(ass_file, seq_file, temperature, pH,
                              lnp_file, times_file, pep, charge_state,
                              exchange, out_file, pi0_file=''):
    """
    Predict and write isotopic envelopes at all timepoints.
    Handles both forward and back-exchange.
    """
    seq = read_seq(seq_file)
    times = read_time_points(times_file)
    ass = read_assignments(ass_file)
    start_res, end_res = ass[int(pep) - 1][1:3]

    if exchange == 'f':
        kint, _ = calculate_kint_for_sequence(1, len(seq), seq, float(temperature), float(pH))
    else:
        kint, _ = calculate_kback_for_sequence(1, len(seq), seq, float(temperature), float(pH))
    kint = kint[start_res:end_res]
    lnP = read_pfact(lnp_file)[start_res:end_res]

    if exchange == 'f':
        pi0 = fully_protonated_envelope(seq[start_res:end_res + 1], charge_state)
        mass = list(pi0.keys())
        fr0 = list(pi0.values())
        while len(mass) <= 2 * len(kint):
            mass.append((mass[-1] + 1.00627 * charge_state) / charge_state)
            fr0.append(0.0)
    else:
        pi0 = pd.read_csv(pi0_file, skiprows=1, header=None, delim_whitespace=True)
        mass = list(pi0[1])
        u_fr0 = list(pi0[2])
        fr0 = centered_isotopic_envelope(0, kint, lnP, u_fr0)

    for i, t in enumerate(times):
        if exchange == 'f':
            f1 = centered_isotopic_envelope(t, kint, lnP, fr0)
        else:
            f1 = back_centered_isotopic_envelope(t, kint, lnP, fr0)

        f1 = [val / sum(f1) * 100 for val in f1]
        with open(f"{out_file}.{i}.isot", 'w') as f:
            f.write(f"# {seq[start_res:end_res]}\n")
            for j, (mz, intensity) in enumerate(zip(mass, f1)):
                rel = intensity / max(f1) * 100 if max(f1) != 0 else 0
                f.write(f"{j}\t{mz:.5f}\t{intensity:.2f}\t{rel:.2f}\n")


def generate_back_exchange_time_points(start=-5, end=5, num=500):
    """
    Generate logarithmically spaced time points for back exchange simulation.
    """
    times = np.logspace(start, end, num)
    with open("back.times", 'w') as f:
        for t in times:
            f.write(f"{t:.15f}\n")


def sticks_from_experimental_envelope(exp_env, corr_env, z):
    """
    Project experimental data into stick format by aligning with theoretical mass axis.
    """
    mass = list(corr_env[1])
    fr = np.zeros(len(mass))
    for mz, intensity in zip(exp_env[0], exp_env[1]):
        for j, ref_mz in enumerate(mass):
            if ref_mz - 0.5 < mz * z < ref_mz + 0.5:
                fr[j] += intensity
    return mass, (fr / sum(fr)) * 100 if sum(fr) > 0 else fr


def compare_predictions(fr, prefix):
    """
    Compare predictions from .isot files against a reference envelope.

    Returns
    -------
    pd.DataFrame
        DataFrame with filenames, timepoints, and R² scores.
    """
    times = read_time_points("back.times")
    files = glob.glob(prefix)
    natural_sort(files)

    scores = []
    for f in files:
        back_env = pd.read_csv(f, header=None, sep='\t', skiprows=1)
        scores.append(r2_score(fr, back_env[2]))

    return pd.DataFrame({"file": files, "time": times, "r2_score": scores})
