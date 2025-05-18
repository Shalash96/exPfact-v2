"""
HDX kinetic constants and parameters for hydrogen-deuterium exchange
in a protonated protein placed in D₂O solution.

Constants are primarily based on:
- Bai et al., 1993 (for backbone exchange)
- Mori et al., 1997 (for D and E side chains)

This module provides temperature and pH-dependent functions for:
- Deuteron/hydroxide concentration
- Temperature activation scaling (Arrhenius)
- pKa shifts of histidine, aspartic acid, and glutamic acid

Refactored from original code by E. Paci group (GPL-2.0).
"""

from math import exp, log10
from typing import Dict, List

# === Amino acid-specific protection factor coefficients ===
# Format: [Acid N-term, Acid Side Chain, Base N-term, Base Side Chain]
# Source: Bai et al., Mori et al.
para: Dict[str, List[float]] = {
    "A": [0.00,   0.00,  0.00,  0.00],
    "C": [-0.54, -0.46,  0.62,  0.55],
    "F": [-0.52, -0.43, -0.24,  0.06],
    "G": [-0.22,  0.22, -0.03,  0.17],
    "I": [-0.91, -0.59, -0.73, -0.23],
    "K": [-0.56, -0.29, -0.04,  0.12],
    "L": [-0.57, -0.13, -0.58, -0.21],
    "M": [-0.64, -0.28, -0.01,  0.11],
    "N": [-0.58, -0.13,  0.49,  0.32],
    "P": [99999, -0.19, 99999, -0.24],  # Trans proline (no exchange)
    "B": [99999, -0.85, 99999,  0.60],  # Cis proline (rare)
    "Q": [-0.47, -0.27,  0.06,  0.20],
    "R": [-0.59, -0.32,  0.08,  0.22],
    "S": [-0.44, -0.39,  0.37,  0.30],
    "T": [-0.79, -0.47, -0.07,  0.20],
    "V": [-0.74, -0.30, -0.70, -0.14],
    "W": [-0.40, -0.44, -0.41, -0.11],
    "Y": [-0.41, -0.37, -0.27,  0.05],
}

# === Terminal corrections ===
rho_Nterm_acid = -1.32       # Acid-catalyzed N-terminus correction
rho_Nterm_base = 1.62        # Base-catalyzed N-terminus correction
lamb_Cterm_acid = 0.96       # Acid-catalyzed C-terminus correction
lamb_Cterm_base = -1.80      # Base-catalyzed C-terminus correction

# === Physical constants ===
pKD = 15.05                  # Autoionization constant of D₂O
R = 1.987                    # Gas constant [cal/(mol·K)]

# === Exchange rate constants in D₂O (s⁻¹) ===
ka_pdla = 10**1.62 / 60      # Acid-catalyzed rate constant (PDLA method)
kb_pdla = 10**10.18 / 60     # Base-catalyzed rate constant (PDLA method)
kw_pdla = 10**-1.5 / 60      # Water-catalyzed rate constant (PDLA method)

ka = 10**2.04 / 60           # Acid-catalyzed rate constant
kb = 10**10.36 / 60          # Base-catalyzed rate constant
kw = 10**-1.5 / 60           # Water-catalyzed rate constant

# === Activation energies (Arrhenius scaling) in cal/mol ===
Ea = 14000  # Acid-catalyzed
Eb = 17000  # Base-catalyzed
Ew = 19000  # Water-catalyzed

# === Reference temperatures in Kelvin ===
T_REF = 293.15              # Default reference temperature
T_PKA_REF = 278.15          # Reference temperature for pKa measurements


# === pH-Dependent Concentrations ===

def get_D(pH: float) -> float:
    """Returns deuteron concentration [D⁺] from pH."""
    return 10 ** (-pH)


def get_OD(pH: float) -> float:
    """Returns hydroxide ion concentration [OD⁻] from pH and pKD."""
    return 10 ** (pH - pKD)


# === Temperature Activation Scaling (Arrhenius) ===

def get_temperature_normalization(temperature: float) -> float:
    """
    Compute normalized inverse temperature difference for Arrhenius factor.
    
    Parameters
    ----------
    temperature : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Normalized temperature: (1/T - 1/Tref) / R
    """
    return (1 / temperature - 1 / T_REF) / R


def get_Fta(temperature: float) -> float:
    """Returns Arrhenius scaling factor for acid-catalyzed exchange at given temperature."""
    return exp(-Ea * get_temperature_normalization(temperature))


def get_Ftb(temperature: float) -> float:
    """Returns Arrhenius scaling factor for base-catalyzed exchange at given temperature."""
    return exp(-Eb * get_temperature_normalization(temperature))


def get_Ftw(temperature: float) -> float:
    """Returns Arrhenius scaling factor for water-catalyzed exchange at given temperature."""
    return exp(-Ew * get_temperature_normalization(temperature))


# === Temperature-Dependent pKa Values ===

def get_pK_his(temperature: float) -> float:
    """
    Compute temperature-dependent pKa for histidine.

    Reference: pKa = 7.42 at 278.15 K, Ea = 7500 cal/mol
    """
    Ea_his = 7500
    return -log10(10**-7.42 * exp(-Ea_his * (1 / temperature - 1 / T_PKA_REF) / R))


def get_pK_asp(temperature: float) -> float:
    """
    Compute temperature-dependent pKa for aspartic acid.

    Reference: pKa = 4.48 at 278.15 K, Ea = 1000 cal/mol
    """
    Ea_asp = 1000
    return -log10(10**-4.48 * exp(-Ea_asp * (1 / temperature - 1 / T_PKA_REF) / R))


def get_pK_glu(temperature: float) -> float:
    """
    Compute temperature-dependent pKa for glutamic acid.

    Reference: pKa = 4.93 at 278.15 K, Ea = 1083 cal/mol
    """
    Ea_glu = 1083
    return -log10(10**-4.93 * exp(-Ea_glu * (1 / temperature - 1 / T_PKA_REF) / R))
