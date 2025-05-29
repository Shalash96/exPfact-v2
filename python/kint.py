"""
Intrinsic Exchange Rate Calculation for HDX-MS
----------------------------------------------

This module provides functions to compute forward intrinsic exchange rates (kint)
for peptides or proteins under given temperature and pH conditions.

Copyright (C) 2019-2020 Emanuele Paci group, code upgraded by Mahmoud Shalash
Licensed under GPL-2.0
"""

from math import log10
import numpy as np
import constants_HD as cst
from typing import List, Tuple


def calculate_kint_for_sequence(
    first_residue: int,
    last_residue: int,
    seq: str,
    temperature: float,
    pH: float,
    reference: str = "3Ala"
) -> Tuple[np.ndarray, List[int]]:
    """
    Compute the forward intrinsic exchange rate (kint) for each residue.

    Parameters
    ----------
    first_residue : int
        Index of the first residue in the peptide (1-based).
    last_residue : int
        Index of the last residue in the peptide (1-based).
    seq : str
        Amino acid sequence.
    temperature : float
        Temperature in Kelvin.
    pH : float
        pH of the solution.
    reference : str
        Exchange model ("3Ala" or "PDLA"). Default is "3Ala".

    Returns
    -------
    Tuple[np.ndarray, List[int]]
        kint : array of intrinsic exchange rates (hr^-1), -1 for prolines
        prolines : list of residue indices that are proline
    """
    prolines = []
    kint = np.full((last_residue,), -1.0)
    res1 = ""
    for i, res2 in enumerate(seq):
        res_index = first_residue + i
        if not res1:
            res1 = res2
            continue
        if res2 in {"P", "B"} or res_index == first_residue:
            kint[res_index - 1] = -1.0
            prolines.append(res_index)
        else:
            kint[res_index - 1] = calculate_kint_per_residue(
                res1, res2, res_index, len(seq), temperature, pH, reference
            )
        res1 = res2

    return kint, prolines


def calculate_kint_per_residue(
    residue1: str,
    residue2: str,
    num: int,
    length: int,
    temperature: float,
    pH: float,
    reference: str = "3Ala"
) -> float:
    """
    Compute intrinsic exchange rate for one residue.

    Parameters
    ----------
    residue1 : str
        Residue at position i-1.
    residue2 : str
        Residue at position i.
    num : int
        Position in sequence (1-based).
    length : int
        Total sequence length.
    temperature : float
        Temperature in Kelvin.
    pH : float
        pH value.
    reference : str
        Reference model: "3Ala" or "PDLA".

    Returns
    -------
    float
        kint value in hr^-1.
    """
    lamb1 = acid(residue2, temperature, pH, "lamb")
    rho1 = acid(residue1, temperature, pH, "rho")
    if num == 2:
        rho1 += cst.rho_Nterm_acid
    elif num == length:
        lamb1 += cst.lamb_Cterm_acid
    Fa = 10 ** (lamb1 + rho1)

    lamb2 = base(residue2, temperature, pH, "lamb")
    rho2 = base(residue1, temperature, pH, "rho")
    if num == 2:
        rho2 += cst.rho_Nterm_base
    elif num == length:
        lamb2 += cst.lamb_Cterm_base
    Fb = 10 ** (lamb2 + rho2)

    if reference == "PDLA":
        return (
            Fa * cst.ka_pdla * cst.get_D(pH) * cst.get_Fta(temperature) * 3600
            + Fb * cst.kb_pdla * cst.get_OD(pH) * cst.get_Ftb(temperature) * 3600
            + Fb * cst.kw_pdla * cst.get_Ftw(temperature) * 3600
        )
    else:
        return (
            Fa * cst.ka * cst.get_D(pH) * cst.get_Fta(temperature) * 3600
            + Fb * cst.kb * cst.get_OD(pH) * cst.get_Ftb(temperature) * 3600
            + Fb * cst.kw * cst.get_Ftw(temperature) * 3600
        )


def acid(residue: str, temperature: float, pH: float, value: str) -> float:
    """
    Get acid-catalyzed parameter (lambda or rho) for a given residue.
    """
    if residue == "H":
        pK = cst.get_pK_his(temperature)
        den = 10**(-pK) + 10**(-pH)
        lamb = log10(10**(-0.80 - pH) / den + 10**(0.00 - pK) / den)
        rho = log10(10**(-0.51 - pH) / den + 10**(0.00 - pK) / den)
    elif residue == "D":
        pK = cst.get_pK_asp(temperature)
        den = 10**(-pK) + 10**(-pH)
        lamb = log10(10**(-0.90 - pH) / den + 10**(0.90 - pK) / den)
        rho = log10(10**(-0.12 - pH) / den + 10**(0.58 - pK) / den)
    elif residue == "E":
        pK = cst.get_pK_glu(temperature)
        den = 10**(-pK) + 10**(-pH)
        lamb = log10(10**(-0.60 - pH) / den + 10**(-0.90 - pK) / den)
        rho = log10(10**(-0.27 - pH) / den + 10**(0.31 - pK) / den)
    else:
        lamb = cst.para[residue][0]
        rho = cst.para[residue][1]

    return lamb if value == "lamb" else rho



def base(residue: str, temperature: float, pH: float, value: str) -> float:
    """
    Get base-catalyzed parameter (lambda or rho) for a given residue.
    """
    if residue == "H":
        pK = cst.get_pK_his(temperature)
        den = 10**(-pK) + 10**(-pH)
        lamb = log10(10**(0.80 - pH) / den + 10**(-0.10 - pK) / den)
        rho = log10(10**(0.83 - pH) / den + 10**(0.14 - pK) / den)
    elif residue == "D":
        pK = cst.get_pK_asp(temperature)
        den = 10**(-pK) + 10**(-pH)
        lamb = log10(10**(0.69 - pH) / den + 10**(0.10 - pK) / den)
        rho = log10(10**(0.60 - pH) / den + 10**(-0.18 - pK) / den)
    elif residue == "E":
        pK = cst.get_pK_glu(temperature)
        den = 10**(-pK) + 10**(-pH)
        lamb = log10(10**(0.24 - pH) / den + 10**(-0.11 - pK) / den)
        rho = log10(10**(0.39 - pH) / den + 10**(-0.15 - pK) / den)
    else:
        lamb = cst.para[residue][2]
        rho = cst.para[residue][3]

    return lamb if value == "lamb" else rho
