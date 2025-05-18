"""
Clean .pfact files by removing lines with IDs listed in 'none.dat'.

For each .pfact file in the current directory:
- Remove any lines where the first column matches an entry in 'none.dat'
- Save the cleaned file with a .cpfact extension

Author: 
    - Original author:
        - Emanuele Paci, Simon P. Skinner, Michele Stofella
    - Updated by:
        Mahmoud Shalash
"""

import os
from pathlib import Path


def load_exclusion_list(filepath: str = "none.dat") -> set:
    """
    Load excluded residue IDs from a file.

    Parameters
    ----------
    filepath : str
        Path to the file containing residue IDs to exclude.

    Returns
    -------
    Set[str]
        A set of string identifiers to be excluded.
    """
    try:
        with open(filepath, "r") as file:
            return {line.strip() for line in file if line.strip()}
    except FileNotFoundError:
        print(f"Error: '{filepath}' not found.")
        return set()


def clean_pfact_files(exclude_ids: set, directory: str = ".") -> None:
    """
    Clean .pfact files in the given directory by removing excluded IDs.

    Parameters
    ----------
    exclude_ids : set
        Set of string residue IDs to exclude.
    directory : str
        Directory containing .pfact files.
    """
    for pfact_path in Path(directory).glob("*.pfact"):
        output_path = pfact_path.with_suffix(".cpfact")

        with pfact_path.open("r") as fin, output_path.open("w") as fout:
            for line in fin:
                fields = line.strip().split()
                if fields and fields[0] not in exclude_ids:
                    fout.write(line)

        print(f"Cleaned: {pfact_path.name} → {output_path.name}")


if __name__ == "__main__":
    exclusion_set = load_exclusion_list("none.dat")
    clean_pfact_files(exclusion_set)
