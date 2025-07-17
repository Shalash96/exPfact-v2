# python/generate_dataset.py

import argparse
import os
import random
import numpy as np
import pandas as pd
from datetime import datetime

def generate_random_protein_sequence(length: int, seed: int = None) -> str:
    amino_acids = 'ACDEFGHIKLMNQRSPTVWY'
    if seed is not None:
        random.seed(seed)
    return ''.join(random.choice(amino_acids) for _ in range(length))

def assign_random_ln_p(sequence: str, seed: int = None) -> pd.DataFrame:
    """Assigns random ln(p) values, with NaN for Prolines."""
    if seed is not None:
        np.random.seed(seed)
    
    ln_p_values = np.full(len(sequence), np.nan)
    
    non_proline_mask = np.array([res != 'P' for res in sequence])
    num_non_prolines = np.sum(non_proline_mask)
    random_values = np.random.uniform(0, 20, num_non_prolines)
    
    ln_p_values[non_proline_mask] = random_values
    
    df = pd.DataFrame({'Residue': list(sequence), 'ln(p)': ln_p_values})
    df.index += 1 # 1-based indexing for residues
    return df
    
def generate_random_peptides(sequence: str, num_peptides: int, min_length: int, max_length: int, seed: int = None) -> pd.DataFrame:
    """Generates a set of unique random peptides from the sequence."""
    if seed is not None:
        random.seed(seed)
    peptides = set()
    max_attempts = num_peptides * 100 
    attempts = 0
    while len(peptides) < num_peptides and attempts < max_attempts:
        start_idx = random.randint(1, len(sequence) - min_length + 1)
        max_possible_length = min(max_length, len(sequence) - start_idx + 1)
        if max_possible_length < min_length:
            attempts += 1
            continue
        length = random.randint(min_length, max_possible_length)
        end_idx = start_idx + length - 1
        peptide = sequence[start_idx - 1:end_idx]
        peptides.add((start_idx, end_idx, peptide))
        attempts += 1

    if len(peptides) < num_peptides:
        print(f"Warning: Only generated {len(peptides)} unique peptides out of the {num_peptides} requested.")

    df = pd.DataFrame(list(peptides), columns=['Start', 'End', 'Sequence'])
    df = df.sort_values(by=['Start', 'End']).reset_index(drop=True)
    df.index += 1 # 1-based indexing for peptides
    return df

def generate_overlapping_di_peptides(sequence: str) -> pd.DataFrame:
    """Generates all possible overlapping dipeptides from the sequence."""
    data = [(i + 1, i + 2, sequence[i:i+2]) for i in range(len(sequence) - 1)]
    df = pd.DataFrame(data, columns=['Start', 'End', 'Sequence'])
    df.index += 1 # 1-based indexing
    return df

def main():
    parser = argparse.ArgumentParser(description="Generate a random dataset for ExPfact testing.")
    parser.add_argument("--length", type=int, default=100, help="Length of the protein sequence.")
    parser.add_argument("--seed", type=int, help="Random seed for reproducibility. (Optional)")
    parser.add_argument("--mode", choices=['random', 'dipeptide', 'both'], default='random', help="Peptide generation mode.")
    parser.add_argument("--min_len", type=int, default=5, help="Min peptide length for random mode.")
    parser.add_argument("--max_len", type=int, default=15, help="Max peptide length for random mode.")
    parser.add_argument("--num_pep", type=int, default=100, help="Number of peptides for random mode.")
    
    args = parser.parse_args()

    # --- Collect all configuration settings ---
    generation_time = datetime.now()
    config = {
        "Generation Time": generation_time.strftime('%Y-%m-%d %H:%M:%S'),
        "Sequence Length": args.length,
        "Random Seed": args.seed if args.seed is not None else "None (random)",
        "Peptide Generation Mode": args.mode.capitalize()
    }
    
    # Add random peptide settings if applicable
    if args.mode in ["random", "both"]:
        config["Min Peptide Length"] = args.min_len
        config["Max Peptide Length"] = args.max_len
        config["Number of Peptides"] = args.num_pep

    # --- Generate Data ---
    sequence = generate_random_protein_sequence(config["Sequence Length"], args.seed)
    # The script is run from the 'results' dir, so 'out_dir' will be created inside it.
    out_dir = f"random_dataset_{generation_time.strftime('%Y-%m-%d_%H-%M-%S')}"
    os.makedirs(out_dir, exist_ok=True)
    print(f"Generating dataset in directory: {os.path.abspath(out_dir)}")

    # --- Write Configuration File ---
    config_path = os.path.join(out_dir, 'config.txt')
    with open(config_path, 'w') as f:
        f.write("--- Data Generation Configuration ---\n")
        for key, value in config.items():
            f.write(f"{key}: {value}\n")
    
    # --- Write Data Files ---
    with open(os.path.join(out_dir, 'sequence.seq'), 'w') as f:
        f.write(sequence)

    ln_p_df = assign_random_ln_p(sequence, args.seed)
    ln_p_df.to_csv(
        os.path.join(out_dir, 'ln_p_values.pfact'), 
        index=True, 
        header=False, 
        columns=['ln(p)'], 
        sep=' ',
        na_rep='nan' # Use 'nan' for prolines
    )

    if args.mode in ["random", "both"]:
        pep_df = generate_random_peptides(
            sequence, 
            config["Number of Peptides"], 
            config["Min Peptide Length"], 
            config["Max Peptide Length"], 
            args.seed
        )
        pep_df.to_csv(os.path.join(out_dir, 'random_peptides.ass'), index=True, header=False, sep=' ')

    if args.mode in ["dipeptide", "both"]:
        dipep_df = generate_overlapping_di_peptides(sequence)
        dipep_df.to_csv(os.path.join(out_dir, 'di_peptide.ass'), index=True, header=False, sep=' ')
        
    print("Dataset generation complete.")

if __name__ == '__main__':
    main()
