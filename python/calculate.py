"""
HDX-MS P-factor Optimization Utilities
======================================

Provides tools to compute root mean square error, evaluate cost function,
run random search and perform optimization of protection factors (lnP).

Copyright (C) 2019-2020 Emanuele Paci, Simon P. Skinner, Michele Stofella
Licensed under GPL-2.0
"""

import numpy as np
from scipy.optimize import minimize, OptimizeResult
from typing import Optional, List, Dict, Tuple
from calc_dpred import calculate_dpred


def calculate_rms(
    dpred: np.ndarray,
    dexp: np.ndarray,
    nj: int,
    weights: Optional[np.ndarray] = None
) -> float:
    """
    Compute the normalized MSE between predicted and experimental D-uptake.

    Note: the name of the function is misleading; it computes the normalized MSE, not the RMSD.

    Parameters
    ----------
    dpred : np.ndarray
        Predicted uptake values.
    dexp : np.ndarray
        Experimental uptake values.
    nj : int
        Number of peptides (normalization factor).
    weights : Optional[np.ndarray]
        Optional weighting vector per data point.

    Returns
    -------
    float
        Normalized root mean square deviation.
    """
    diff_sq = (dpred - dexp) ** 2
    if weights is not None:
        return np.sum(weights * diff_sq) / nj
    return np.sum(diff_sq) / nj


def harmonic_score(params: List[float], k: float) -> float:
    """
    Harmonic regularization term to smooth neighboring protection factors.

    Parameters
    ----------
    params : list of float
        Protection factor estimates.
    k : float
        Harmonic regularization constant.

    Returns
    -------
    float
        Harmonic penalty score.
    """
    return sum(
        k * (params[i - 1] - 2 * params[i] + params[i + 1]) ** 2
        for i in range(1, len(params) - 1)
        if params[i - 1] >= 0 and params[i] >= 0 and params[i + 1] >= 0
    )


def cost_function(params: List[float], *args) -> float:
    """
    Cost function combining RMSD and harmonic penalty.

    Parameters
    ----------
    params : list of float
        Current estimates of lnP.
    args : tuple
        Arguments required to compute cost:
        (dexp, time_points, assignments, k, kint, weights)

    Returns
    -------
    float
        Total cost score.
    """
    dexp, time_points, assignments, k, kint, weights = args
    dpred = calculate_dpred(np.array(params), time_points, kint, assignments)
    rms = calculate_rms(dpred, dexp, len(assignments), weights)
    return rms + harmonic_score(params, k)


def do_random_search(
    kint: np.ndarray,
    search_steps: int,
    pfactor_filter: set,
    dexp: np.ndarray,
    time_points: np.ndarray,
    assignments: np.ndarray,
    harmonic_term: float,
    prolines: set,
    weights: Optional[np.ndarray],
    seed: Optional[int] = None
) -> Dict[float, List[float]]:
    """
    Perform a brute-force random search for initial protection factors.

    Parameters
    ----------
    kint : np.ndarray
        Intrinsic exchange rates.
    search_steps : int
        Number of random trials to perform.
    pfactor_filter : set
        Set of residue indices to fit.
    dexp : np.ndarray
        Experimental deuterium uptake.
    time_points : np.ndarray
        HDX labeling timepoints.
    assignments : np.ndarray
        Peptide assignment matrix.
    harmonic_term : float
        Harmonic penalty constant.
    prolines : set
        Residue indices of prolines (non-exchanging).
    weights : Optional[np.ndarray]
        Optional data point weights.
    seed : Optional[int]
        Random seed for reproducibility.

    Returns
    -------
    Dict[float, List[float]]
        Dictionary mapping cost scores to pfactor arrays.
    """
    if seed is None and search_steps == 1:
        np.random.seed(42)

    results = {}
    for _ in range(search_steps):
        init = [
            np.random.uniform(0.01, 30.0) if i != 0 and (i + 1) not in prolines and (i + 1) in pfactor_filter else -1
            for i in range(len(kint))
        ]
        score = cost_function(init, dexp, time_points, assignments, harmonic_term, kint, weights)
        results[score] = init

    return results


def fit_pfact(
    init_array: List[float],
    dexp: np.ndarray,
    time_points: np.ndarray,
    assignments: np.ndarray,
    harmonic_term: float,
    kint: np.ndarray,
    bounds: List[Tuple[float, float]],
    tol: float,
    weights: Optional[np.ndarray]
) -> OptimizeResult:
    """
    Perform optimization using L-BFGS-B to fit protection factors.

    Parameters
    ----------
    init_array : List[float]
        Initial guess of lnP values.
    dexp : np.ndarray
        Experimental uptake.
    time_points : np.ndarray
        Labeling timepoints.
    assignments : np.ndarray
        Peptide-to-residue assignment.
    harmonic_term : float
        Smoothing regularization strength.
    kint : np.ndarray
        Intrinsic exchange rates.
    bounds : List[Tuple[float, float]]
        Min-max bounds for each parameter.
    tol : float
        Tolerance for convergence.
    weights : Optional[np.ndarray]
        Optional per-point weights.

    Returns
    -------
    scipy.optimize.OptimizeResult
        Optimization result containing final lnP values.
    """
    result = minimize(
        cost_function,
        init_array,
        args=(dexp, time_points, assignments, harmonic_term, kint, weights),
        method="L-BFGS-B",
        bounds=bounds,
        tol=tol,
        options={
            "disp": False,
            "maxfun": 1_000_000_000,
            "maxiter": 1_000_000_000,
        }
    )
    return result


def predict_dexp(
    pfact: List[float],
    time_points: np.ndarray,
    kint: np.ndarray,
    assignments: np.ndarray
) -> np.ndarray:
    """
    Generate predicted uptake (Dexp) from lnP values.

    Parameters
    ----------
    pfact : List[float]
        Protection factors (lnP).
    time_points : np.ndarray
        HDX time points.
    kint : np.ndarray
        Intrinsic exchange rates.
    assignments : np.ndarray
        Peptide-to-residue map.

    Returns
    -------
    np.ndarray
        Predicted uptake values.
    """
    return calculate_dpred(np.array(pfact), time_points, kint, assignments)

