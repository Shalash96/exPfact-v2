"""
MD2Pfact

Calculate protection factors from a Molecular Dynamics trajectory using the
Best-Vendruscolo model:

    < ln(P) > = beta_c * < N_c > + beta_h * < N_h >

N_c: heavy atom contacts
N_h: hydrogen bonds

Usage:
    python MD2Pfact.py --pdb input.pdb --dcd input.dcd --out output_folder
                       [--bc beta_c --bh beta_h --step step --no-plot]
"""

import os
import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis
from MDAnalysis.analysis.distances import distance_array
from Bio.SeqUtils import IUPACData
from typing import Tuple, List


def contacts_smooth_function(dist: np.ndarray) -> np.ndarray:
    """Smooth sigmoid function for heavy atom contacts"""
    return 1 / (1 + np.exp(5 * (dist - 6.5)))


def hydrogen_bonds_smooth_function(dist: np.ndarray) -> np.ndarray:
    """Smooth sigmoid function for hydrogen bond distances"""
    return 1 / (1 + np.exp(10 * (dist - 2.4)))


def evaluate_pfact(Nc: float, Nh: float, A: float = 0.35, B: float = 2.00) -> float:
    return A * Nc + B * Nh


def calculate_pfact_from_trajectory(
    PDB: str,
    DCD: str,
    OUT: str,
    step: int = 100,
    bh: float = 2.00,
    bc: float = 0.35
) -> Tuple[List[int], np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate protection factors from a trajectory using the Best-Vendruscolo model.
    """
    u = MDAnalysis.Universe(PDB, DCD)
    chains = list(set(u.select_atoms("protein").segids))
    monomer = len(chains) == 1
    residues = sorted(set(u.select_atoms("protein").residues.resnums))
    sequence = u.select_atoms("protein").residues.resnames
    seq = "".join(IUPACData.protein_letters_3to1.get(res if res != "HSD" else "HIS", "X") for res in sequence)

    with open(PDB.replace(".pdb", ".seq"), 'w') as f:
        f.write(seq)

    pfact = []

    for ts in u.trajectory[::step]:
        pfact_t = []
        for chain in chains:
            for r in residues:
                try:
                    amide_nitrogen = u.select_atoms(f"protein and segid {chain} and name N and resid {r}")
                    exclude = u.select_atoms(
                        f"not resid {r} and not resid {r - 1} and not resid {r + 1} and "
                        f"not resid {r - 2} and not resid {r + 2}"
                    )
                    O_list = exclude.select_atoms("type O")
                    H_list = exclude.select_atoms("not type H")

                    dist_O = distance_array(amide_nitrogen, O_list)[0]
                    dist_H = distance_array(amide_nitrogen, H_list)[0]

                    Nc = np.sum(contacts_smooth_function(dist_H))
                    Nh = np.sum(hydrogen_bonds_smooth_function(dist_O))

                    pfact_t.append(evaluate_pfact(Nc, Nh, bc, bh))
                except Exception:
                    pfact_t.append(0.0)  # Fallback if atom missing
        pfact.append(pfact_t)

    p_all = np.array(pfact)
    p_avg = np.mean(p_all, axis=0)
    p_std = np.std(p_all, axis=0)

    np.savetxt(f"{OUT}_all.txt", p_all, fmt="%.5e")

    if not monomer:
        residues = list(range(len(residues) * len(chains)))
    else:
        residues = [i for i in range(len(p_avg))]

    with open(f"{OUT}_avg.txt", 'w') as f:
        for i in range(len(residues)):
            f.write(f"{residues[i] + 1} {p_avg[i]:.5e} {p_std[i]:.5e}\n")

    return residues, p_avg, p_std, p_all


def main():
    parser = argparse.ArgumentParser(description="Calculate ln(P) from MD trajectory.")
    parser.add_argument("--pdb", required=True, help="Input PDB file.")
    parser.add_argument("--dcd", required=True, help="Input DCD trajectory.")
    parser.add_argument("--out", required=True, help="Output filename prefix.")
    parser.add_argument("--step", type=int, default=100, help="Step size in trajectory frames.")
    parser.add_argument("--bc", type=float, default=0.35, help="Coefficient for heavy contacts.")
    parser.add_argument("--bh", type=float, default=2.00, help="Coefficient for hydrogen bonds.")
    parser.add_argument("--no-plot", action="store_true", help="Suppress final plot.")

    args = parser.parse_args()

    start = time.time()
    residues, p_avg, p_std, _ = calculate_pfact_from_trajectory(
        PDB=args.pdb, DCD=args.dcd, OUT=args.out,
        step=args.step, bh=args.bh, bc=args.bc
    )
    end = time.time()

    print(f"Elapsed time: {end - start:.2f} s")

    if not args.no_plot:
        plt.figure()
        plt.title(f"Protection Factors from {args.dcd}", fontsize=15)
        plt.errorbar(residues, p_avg, p_std, fmt='o-', color='black', capsize=2)
        plt.xlabel("Residue Index", fontsize=15)
        plt.ylabel("<ln(P)>", fontsize=15)
        plt.ylim(-1, 21)
        plt.tight_layout()
        plt.savefig(f"{args.out}_avg.png", dpi=300)
        plt.show()


if __name__ == '__main__':
    main()
