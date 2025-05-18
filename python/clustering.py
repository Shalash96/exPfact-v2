"""
Clustering of HDX-MS contiguous peptide regions using Mclust (via Rscript).

- Reads peptide assignments from a file
- Identifies contiguous regions with no peptide coverage
- Calls `multi.r` R script to cluster each region using Mclust

Author: Code updated by Mahmoud Shalash based on original code by E. Paci group
"""

import argparse
import subprocess
from pathlib import Path
import numpy as np
from read import read_assignments
from typing import List

import numpy as np
from typing import List

def find_covered_regions_between_gaps(assignments: np.ndarray) -> List[str]:
    """
    Identifies contiguous covered regions based on gaps between uncovered residues.
    Optimized for computational efficiency.

    Parameters
    ----------
    assignments : np.ndarray
        Peptide assignment array with shape (N, 3), where columns are:
        [peptide_index, start_residue (1-based), end_residue (1-based)]

    Returns
    -------
    List[str]
        List of contiguous covered regions (as 'start-end' strings), inferred between uncovered blocks.
    """
    # Use more efficient vectorized operation to find max value
    max_residue = int(np.max(assignments[:, 2]))
    
    # Use boolean array instead of integers for coverage - more memory efficient
    coverage = np.zeros(max_residue, dtype=bool)
    
    # Vectorized approach to fill coverage
    # Extract start and end positions
    starts = assignments[:, 1].astype(int) - 1  # Convert to 0-based
    ends = assignments[:, 2].astype(int)        # No subtraction needed for end
    
    # Fill coverage using a single loop instead of nested loops
    for start, end in zip(starts, ends):
        coverage[start:end] = True
    
    # Find uncovered residue indices
    uncovered_indices = np.where(~coverage)[0]  # Using ~ is faster than == False
    
    # Quick return for full coverage
    if len(uncovered_indices) == 0:
        return [f"1-{max_residue}"]
    
    # Convert to 1-based residue numbers
    uncovered_residues = uncovered_indices + 1
    
    # Pre-allocate results list with estimated size to avoid resizing
    # Worst case: alternating covered/uncovered regions
    regions = []
    
    # Handle first region if exists
    if uncovered_residues[0] > 1:
        regions.append(f"1-{uncovered_residues[0] - 1}")
    
    # Find gaps between uncovered residues
    # Use numpy operations to identify consecutive differences > 1
    if len(uncovered_residues) > 1:
        # Calculate gaps between consecutive uncovered residues
        gaps = np.diff(uncovered_residues)
        # Find where gaps are larger than 1
        gap_indices = np.where(gaps > 1)[0]
        
        # Process each gap to create region strings
        for i in gap_indices:
            start = uncovered_residues[i] + 1
            end = uncovered_residues[i + 1] - 1
            regions.append(f"{start}-{end}")
    
    # Handle last region if exists
    if uncovered_residues[-1] < max_residue:
        regions.append(f"{uncovered_residues[-1] + 1}-{max_residue}")
    
    return regions




def run_mclust_on_regions(regions: list[str], r_script: str = "../R/multi.r") -> None:
    """
    Run Mclust clustering via Rscript for each uncovered region.

    Parameters
    ----------
    regions : List[str]
        List of residue ranges in format 'start-end'.
    r_script : str
        Path to the R clustering script.
    """
    for region in regions:
        start, end = region.split('-')
        command = ["Rscript", r_script, start, end]
        try:
            subprocess.run(command, check=True)
            print(f"Clustered region: {start}-{end}")
        except subprocess.CalledProcessError:
            print(f"Failed clustering for region: {start}-{end}")


def main():
    parser = argparse.ArgumentParser(description="Run Mclust clustering on uncovered HDX-MS regions.")
    parser.add_argument("--ass", required=True, help="Path to peptide assignment file (.list)")
    args = parser.parse_args()

    assignments = read_assignments(args.ass)
    regions = find_covered_regions_between_gaps(assignments)
    run_mclust_on_regions(regions)


if __name__ == "__main__":
    main()
