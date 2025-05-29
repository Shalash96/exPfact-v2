"""
Back-exchange rate (kback) calculator for deuterated proteins in H2O.

Based on the work of Bai et al. (1993), Connelly et al. (1993), Nguyen et al. (2018)

Copyright (C) 2019-2020 Emanuele Paci, Simon P. Skinner, Michele Stofella
GPL-2.0 License
"""

from math import log10
import numpy as np
from typing import List, Tuple
import constants_DH as cst


def calculate_kback_for_sequence(
    first_residue: int,
    last_residue: int,
    seq: str,
    temperature: float,
    pH: float
) -> Tuple[np.ndarray, List[int]]:
    """
    Calculate intrinsic back-exchange rates (kback) for a sequence in H2O.

    Parameters
    ----------
    first_residue : int
        Index of the first residue in the sequence.
    last_residue : int
        Index of the last residue in the sequence.
    seq : str
        Amino acid sequence (one-letter code).
    temperature : float
        Temperature in Kelvin.
    pH : float
        Corrected pD (pD = pD_read + 0.4)

    Returns
    -------
    Tuple[np.ndarray, List[int]]
        kback : array of exchange rates (hr⁻¹), -1 for prolines.
        prolines : list of residue indices corresponding to prolines.
    """
    kback = np.full(last_residue, -1.0)
    prolines: List[int] = []

    res1 = ""
    for idx, res in enumerate(seq, start=1):
        if res1:
            if idx - first_residue == 0:
                kback[idx - 1] = -1
            elif res in {"P", "B"}:
                kback[idx - 1] = -1
                prolines.append(first_residue + idx - 1)
            else:
                kback[idx - 1] = calculate_kback_per_residue(
                    res1, res, idx, len(seq), temperature, pH
                )
        res1 = res

    return kback, prolines


def calculate_kback_per_residue(
    residue1: str,
    residue2: str,
    num: int,
    length: int,
    temperature: float,
    pH: float
) -> float:
    """
    Compute kback for a single residue pair.

    Parameters
    ----------
    residue1 : str
        Residue at position i-1.
    residue2 : str
        Residue at position i.
    num : int
        Position of the current residue.
    length : int
        Length of the full sequence.
    temperature : float
        Temperature in Kelvin.
    pH : float
        Corrected pD value.

    Returns
    -------
    float
        Back-exchange rate for the residue (hr⁻¹).
    """
    lamb1 = acid(residue2, temperature, pH, "lamb")
    rho1 = acid(residue1, temperature, pH, "rho")

    if num == 2:
        rho1 += cst.rho_Nterm_acid
    elif num == length:
        lamb1 += cst.lamb_Cterm_acid

    Fa = 10**(lamb1 + rho1)

    lamb2 = base(residue2, temperature, pH, "lamb")
    rho2 = base(residue1, temperature, pH, "rho")

    if num == 2:
        rho2 += cst.rho_Nterm_base
    elif num == length:
        lamb2 += cst.lamb_Cterm_base

    Fb = 10**(lamb2 + rho2)

    # Convert from s⁻¹ to hr⁻¹ by multiplying with 3600
    kback = (
        Fa * cst.ka * cst.get_D(pH) * cst.get_Fta(temperature) * 3600 +
        Fb * cst.kb * cst.get_OD(pH) * cst.get_Ftb(temperature) * 3600 +
        Fb * cst.kw * cst.get_Ftw(temperature) * 3600
    )
    return kback


def acid(residue: str, temperature: float, pH: float, value: str) -> float:
    """
    Retrieve acid-catalyzed kinetic parameter for a residue.

    Parameters
    ----------
    residue : str
    temperature : float
    pH : float
    value : str
        Either "lamb" or "rho"

    Returns
    -------
    float
        Corresponding kinetic modifier
    """
    if residue == "H":
        pK = cst.get_pK_his(temperature)
        lamb = log10((10**(-0.80 - pH) + 10**(0.00 - pK)) / (10**(-pK) + 10**(-pH)))
        rho = log10((10**(-0.51 - pH) + 10**(0.00 - pK)) / (10**(-pK) + 10**(-pH)))
    elif residue == "D":
        pK = cst.get_pK_asp(temperature)
        lamb = log10((10**(-0.90 - pH) + 10**(0.90 - pK)) / (10**(-pK) + 10**(-pH)))
        rho = log10((10**(-0.12 - pH) + 10**(0.58 - pK)) / (10**(-pK) + 10**(-pH)))
    elif residue == "E":
        pK = cst.get_pK_glu(temperature)
        lamb = log10((10**(-0.60 - pH) + 10**(-0.90 - pK)) / (10**(-pK) + 10**(-pH)))
        rho = log10((10**(-0.27 - pH) + 10**(0.31 - pK)) / (10**(-pK) + 10**(-pH)))
    else:
        lamb, rho = cst.para[residue][:2]

    return lamb if value == "lamb" else rho


def base(residue: str, temperature: float, pH: float, value: str) -> float:
    """
    Retrieve base-catalyzed kinetic parameter for a residue.

    Parameters
    ----------
    residue : str
    temperature : float
    pH : float
    value : str
        Either "lamb" or "rho"

    Returns
    -------
    float
        Corresponding kinetic modifier
    """
    if residue == "H":
        pK = cst.get_pK_his(temperature)
        lamb = log10((10**(0.80 - pH) + 10**(-0.10 - pK)) / (10**(-pK) + 10**(-pH)))
        rho = log10((10**(0.83 - pH) + 10**(0.14 - pK)) / (10**(-pK) + 10**(-pH)))
    elif residue == "D":
        pK = cst.get_pK_asp(temperature)
        lamb = log10((10**(0.69 - pH) + 10**(0.10 - pK)) / (10**(-pK) + 10**(-pH)))
        rho = log10((10**(0.60 - pH) + 10**(-0.18 - pK)) / (10**(-pK) + 10**(-pH)))
    elif residue == "E":
        pK = cst.get_pK_glu(temperature)
        lamb = log10((10**(0.24 - pH) + 10**(-0.11 - pK)) / (10**(-pK) + 10**(-pH)))
        rho = log10((10**(0.39 - pH) + 10**(-0.15 - pK)) / (10**(-pK) + 10**(-pH)))
    else:
        lamb, rho = cst.para[residue][2:]

    return lamb if value == "lamb" else rho
