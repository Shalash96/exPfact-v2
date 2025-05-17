# cython: boundscheck=False, wraparound=False, cdivision=True

from math import exp
import numpy as np
cimport numpy as cnp

# Ensure numpy memory views are declared properly
ctypedef cnp.double_t DTYPE_t
ctypedef cnp.int64_t ITYPE_t


def get_residue_rates(cnp.ndarray[DTYPE_t, ndim=1] kint,
                      cnp.ndarray[DTYPE_t, ndim=1] lnP) -> cnp.ndarray:
    """
    Calculate exchange rate per residue as kint / exp(P)
    - If kint[i] == -1, that residue is skipped (e.g., Prolines)

    Parameters
    ----------
    kint : array of intrinsic exchange rates
    lnP  : array of log protection factors

    Returns
    -------
    rates : array of actual exchange rates per residue
    """
    cdef int n_residues = kint.shape[0]
    cdef cnp.ndarray[DTYPE_t, ndim=1] rates = np.zeros(n_residues)

    cdef int i
    for i in range(n_residues):
        if kint[i] == -1:
            rates[i] = -1.0
        else:
            rates[i] = kint[i] / exp(lnP[i])

    return rates


def compute_peptide_uptake(cnp.ndarray[DTYPE_t, ndim=1] residue_rates,
                           cnp.ndarray[ITYPE_t, ndim=2] assignments,
                           cnp.ndarray[DTYPE_t, ndim=1] time_points) -> cnp.ndarray:
    """
    Calculate deuterium uptake for each peptide at each time point.

    Parameters
    ----------
    residue_rates : array of per-residue exchange rates
    assignments    : [n_fragments x 3] array with [index, start, end] per peptide
    time_points    : exposure times in hours (or seconds)

    Returns
    -------
    dpred : 2D array of predicted uptake [n_peptides x n_timepoints]
    """
    cdef int n_peptides = assignments.shape[0]
    cdef int n_times = time_points.shape[0]

    cdef cnp.ndarray[DTYPE_t, ndim=2] dpred = np.zeros((n_peptides, n_times))
    cdef cnp.ndarray[DTYPE_t, ndim=1] n_amides = np.zeros(n_peptides)

    cdef int i, j, k, start, end

    # Count exchangeable amides for each peptide
    for i in range(n_peptides):
        start = assignments[i, 1]
        end = assignments[i, 2]
        for j in range(start, end):
            if residue_rates[j] >= 0:
                n_amides[i] += 1

    # Compute uptake for each peptide/time
    for i in range(n_peptides):
        start = assignments[i, 1]
        end = assignments[i, 2]
        for k in range(n_times):
            for j in range(start, end):
                if residue_rates[j] >= 0:
                    dpred[i, k] += exp(-residue_rates[j] * time_points[k])
            dpred[i, k] = (n_amides[i] - dpred[i, k]) / n_amides[i]

    return dpred


def calculate_dpred(cnp.ndarray[DTYPE_t, ndim=1] lnP,
                    cnp.ndarray[DTYPE_t, ndim=1] time_points,
                    cnp.ndarray[DTYPE_t, ndim=1] kint,
                    cnp.ndarray[ITYPE_t, ndim=2] assignments) -> cnp.ndarray:
    """
    Main function to compute predicted uptake matrix.

    Parameters
    ----------
    lnP         : natural log protection factors
    time_points : array of exposure times
    kint        : intrinsic exchange rates
    assignments : peptide fragment assignment matrix

    Returns
    -------
    dpred : predicted D-uptake per peptide per time point
    """
    cdef cnp.ndarray[DTYPE_t, ndim=1] rates = get_residue_rates(kint, lnP)
    cdef cnp.ndarray[DTYPE_t, ndim=2] dpred = compute_peptide_uptake(rates, assignments, time_points)
    return dpred
