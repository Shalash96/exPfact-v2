"""
Physical and empirical constants for modeling HDX intrinsic rates
in protonated conditions (H₂O), based on experimental work from Bai et al. (1993).

This module defines:
- Amino acid-dependent intrinsic rate coefficients (`para`)
- Terminal modifications (N-term and C-term)
- Exchange rate constants (ka, kb, kw)
- Activation energies for temperature scaling
- pKa shifts for side chains (His, Asp, Glu)
- Functions to compute D⁺/OD⁻ concentration and temperature correction factors

Refactored from original code by E. Paci group (GPL-2.0).
"""

from math import exp, log10

# === Intrinsic rate modifiers per residue (unitless scaling coefficients)
# Format: [N_acid, S_acid, N_base, S_base]
# Reference: Bai et al., 1993
para = {
    "A": [0.00,   0.00,  0.00,  0.00],
    "C": [-0.54, -0.46,  0.62,  0.55],
    "F": [-0.52, -0.43, -0.24,  0.06],
    "G": [-0.22,  0.22, -0.03,  0.17],
    "I": [-0.91, -0.59, -0.73, -0.23],
    "K": [-0.56, -0.29, -0.04,  0.12],
    "L": [-0.57, -0.13, -0.58, -0.21],
    "M": [-0.64, -0.28, -0.01,  0.11],
    "N": [-0.58, -0.13,  0.49,  0.32],
    "P": [99999, -0.19, 99999, -0.24],  # Prolines excluded from exchange
    "Q": [-0.47, -0.27,  0.06,  0.20],
    "R": [-0.59, -0.32,  0.08,  0.22],
    "S": [-0.44, -0.39,  0.37,  0.30],
    "T": [-0.79, -0.47, -0.07,  0.20],
    "V": [-0.74, -0.30, -0.70, -0.14],
    "W": [-0.40, -0.44, -0.41, -0.11],
    "Y": [-0.41, -0.37, -0.27,  0.05],
}

# === Terminal correction factors
rho_Nterm_acid = -1.32      # N-terminal acid term
rho_Nterm_base = 1.62       # N-terminal base term

lamb_Cterm_acid = 0.08      # C-terminal acid term
lamb_Cterm_base = -1.80     # C-terminal base term

# === Physical constants
pKD = 14.17                 # Water autoionization pKa in H₂O
R = 1.987                   # Gas constant [cal/(mol·K)]

# === Rate constants [s⁻¹] — scaled from literature values
ka = 10**(1.4) / 60         # Acid-catalyzed exchange rate
kb = 10**(10.0) / 60        # Base-catalyzed exchange rate
kw = 10**(-1.6) / 60        # Water-catalyzed exchange rate

# === Activation energies [cal/mol]
Ea = 14000                 # Acid-catalyzed
Eb = 17000                 # Base-catalyzed
Ew = 19000                 # Water-catalyzed

# === Functional forms ===

def get_D(pH: float) -> float:
    """Return [D⁺] concentration (mol/L) given pH."""
    return 10 ** (-pH)


def get_OD(pH: float) -> float:
    """Return [OD⁻] concentration (mol/L) from pH and pK_D."""
    return 10 ** (pH - pKD)


def get_temperature_normalization(temperature: float) -> float:
    """
    Return Arrhenius normalization factor: (1/T - 1/Tref) / R.
    
    Parameters
    ----------
    temperature : float
        Temperature in Kelvin. Reference T = 293 K.
    """
    return (1 / temperature - 1 / 293) / R


def get_pK_his(temperature: float) -> float:
    """
    Temperature-dependent pKa for histidine side chain.
    Reference pKa = 7.00 at 278 K, Ea = 7500 cal/mol.
    """
    Ea_his = 7500
    return -log10(10**(-7.00) * exp(-Ea_his * (1 / temperature - 1 / 278) / R))


def get_pK_asp(temperature: float) -> float:
    """
    Temperature-dependent pKa for aspartic acid side chain.
    Reference pKa = 3.87 at 278 K, Ea = 960 cal/mol.
    """
    Ea_asp = 960
    return -log10(10**(-3.87) * exp(-Ea_asp * (1 / temperature - 1 / 278) / R))


def get_pK_glu(temperature: float) -> float:
    """
    Temperature-dependent pKa for glutamic acid side chain.
    Reference pKa = 4.33 at 278 K, Ea = 1083 cal/mol.
    """
    Ea_glu = 1083
    return -log10(10**(-4.33) * exp(-Ea_glu * (1 / temperature - 1 / 278) / R))


def get_Fta(temperature: float) -> float:
    """
    Return acid-catalyzed Arrhenius temperature scaling factor.
    Based on Ea = 14 kcal/mol.
    """
    return exp(-Ea * get_temperature_normalization(temperature))


def get_Ftb(temperature: float) -> float:
    """
    Return base-catalyzed Arrhenius temperature scaling factor.
    Based on Eb = 17 kcal/mol.
    """
    return exp(-Eb * get_temperature_normalization(temperature))


def get_Ftw(temperature: float) -> float:
    """
    Return water-catalyzed Arrhenius temperature scaling factor.
    Based on Ew = 19 kcal/mol.
    """
    return exp(-Ew * get_temperature_normalization(temperature))
