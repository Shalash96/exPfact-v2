"""
Compute coarse isotopic envelope of fully protonated peptide using pyOpenMS.

Author: E. Paci group (GPL-2.0), upgraded by Mahmoud Shalash
"""

import os
import sys
from typing import Dict, Optional
from pyopenms import AASequence, CoarseIsotopePatternGenerator
from logger import log
import argparse


def fully_protonated_envelope(sequence: str, z: int, write: bool = True) -> Optional[Dict[float, float]]:
    """
    Compute the coarse isotopic envelope of a fully protonated peptide.

    Parameters
    ----------
    sequence : str
        Amino acid sequence of the peptide (one-letter code).
    z : int
        Charge state of the peptide.
    write : bool, optional
        Whether to write the isotopic envelope to a file named '<sequence>.txt'. Default is True.

    Returns
    -------
    Optional[Dict[float, float]]
        Dictionary mapping m/z values to relative intensities (sum-normalized to 100).
        Returns None if an error occurs.
    """
    if not sequence:
        log.error("Input sequence is empty.")
        return None
    if z <= 0:
        log.error(f"Charge state must be positive. Got: {z}")
        return None

    try:
        seq = AASequence.fromString(sequence)
        seq_formula = seq.getFormula()
    except Exception as e:
        log.error(f"Invalid peptide sequence '{sequence}': {e}")
        return None

    try:
        isotope_generator = CoarseIsotopePatternGenerator(2 * len(sequence))
        isotopes = seq_formula.getIsotopeDistribution(isotope_generator)
    except Exception as e:
        log.error(f"Error generating isotopic distribution: {e}")
        return None

    try:
        isotopic_envelope = {
            (iso.getMZ() + z) / z: iso.getIntensity() * 100
            for iso in isotopes.getContainer()
        }
        if not isotopic_envelope:
            log.warning("No isotopic peaks were generated.")
            return {}
    except Exception as e:
        log.error(f"Error constructing isotopic envelope: {e}")
        return None

    if write:
        file_name = f"{sequence}.txt"
        try:
            max_intensity = max(isotopic_envelope.values()) or 1.0
            with open(file_name, 'w') as f:
                f.write(f"# {sequence}\n")
                for i, (mz, intensity_sum_norm) in enumerate(sorted(isotopic_envelope.items())):
                    intensity_max_norm = intensity_sum_norm / max_intensity * 100
                    f.write(f"{i} {mz:.5f} {intensity_sum_norm:.2f} {intensity_max_norm:.2f}\n")
            log.info(f"Fully protonated envelope saved in file: {file_name}")
        except IOError as e:
            log.error(f"Failed to write output file '{file_name}': {e}")
        except Exception as e:
            log.error(f"Unexpected error during file write: {e}")

    return isotopic_envelope


def main() -> None:
    """
    Main entry point when running this script directly.
    Parses command-line arguments for peptide sequence and charge.
    """
    parser = argparse.ArgumentParser(description="Compute isotopic envelope for a protonated peptide.")
    parser.add_argument("--seq", required=True, help="Peptide sequence")
    parser.add_argument("--z", required=True, type=int, help="Charge state (positive integer)")
    args = parser.parse_args()

    log.info("Running Hisotope.py")
    result = fully_protonated_envelope(sequence=args.seq, z=args.z)
    if result is None:
        log.error("Isotopic envelope calculation failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
