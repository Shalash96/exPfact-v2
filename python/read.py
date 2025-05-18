"""
File readers for HDX-MS input files.

Includes utilities to read assignment files, sequences, protection factors,
intrinsic exchange rates, time points, and configuration files.

Copyright (C) 2019–2020 Emanuele Paci, Simon P. Skinner, Michele Stofella
Upgraded by Mahmoud Shalash
Licensed under GPL-2.0
"""

import json
import numpy as np
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union


def read_assignments(assignment_file: Union[str, Path]) -> np.ndarray:
    """
    Reads peptide assignment data.

    Parameters
    ----------
    assignment_file : str or Path
        Path to the assignment file.

    Returns
    -------
    np.ndarray
        2D array of assignments [index, start, end].
    """
    assignments = []
    with open(assignment_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                assignments.append([int(parts[0]), int(parts[1]), int(parts[2])])
    return np.array(assignments, dtype=int)


def read_seq(seq_file: Union[str, Path]) -> str:
    """
    Reads a protein sequence from a file.

    Parameters
    ----------
    seq_file : str or Path
        Path to the sequence file.

    Returns
    -------
    str
        Amino acid sequence as a string.
    """
    with open(seq_file, 'r') as f:
        return f.read().strip()


def read_pfact(pfact_file: Union[str, Path]) -> np.ndarray:
    """
    Reads protection factors.

    Parameters
    ----------
    pfact_file : str or Path
        Path to the protection factor file.

    Returns
    -------
    np.ndarray
        1D array of ln(P) values.
    """
    values = []
    with open(pfact_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                values.append(float(parts[1]))
    return np.array(values, dtype=float)


def read_kint(kint_file: Union[str, Path]) -> np.ndarray:
    """
    Reads intrinsic exchange rates.

    Parameters
    ----------
    kint_file : str or Path
        Path to the kint file.

    Returns
    -------
    np.ndarray
        1D array of intrinsic rates.
    """
    values = []
    with open(kint_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                values.append(float(parts[1]))
    return np.array(values, dtype=float)


def read_time_points(time_points_file: Union[str, Path]) -> np.ndarray:
    """
    Reads time points from a text file.

    Parameters
    ----------
    time_points_file : str or Path

    Returns
    -------
    np.ndarray
        1D array of time points.
    """
    return np.loadtxt(time_points_file, dtype=float)


def read_dexp(dexp_file: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reads D-uptake values and associated time points from a .dexp file.

    Parameters
    ----------
    dexp_file : str or Path

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Transposed uptake matrix [peptides x timepoints],
        and corresponding 1D array of time points.
    """
    uptake_data = []
    time_data = []

    with open(dexp_file, 'r') as f:
        for line in f:
            tokens = line.strip().split()
            if len(tokens) >= 2:
                time_data.append(float(tokens[0]))
                uptake_values = [float(val) for val in tokens[1:]]
                uptake_data.append(uptake_values)

    return np.array(uptake_data).T, np.array(time_data, dtype=float)


def read_configuration(config_file: Union[str, Path]) -> Dict[str, Any]:
    """
    Reads a JSON config file and loads associated data files.

    Parameters
    ----------
    config_file : str or Path
        Path to a configuration JSON file.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing loaded configuration and data arrays.
    """
    with open(config_file, 'r') as f:
        config_json = json.load(f)

    base_path = Path(config_json['base'])
    get_path = lambda key: base_path / config_json[key]

    config: Dict[str, Any] = {
        'base': str(base_path),
        'assignments': read_assignments(get_path('assignments_file')),
        'dexp': None,
        'times': None,
        'harmonic_factor': config_json['harmonic_term'],
        'kint': read_kint(get_path('kint_file')),
        'output': config_json['output_file'],
        'pfact': '',
        'predict': config_json['predict'],
        'random_search': config_json['do_random_search'],
        'random_search_steps': config_json['random_search_steps'],
        'time_points': config_json['time_points_file'],
        'tolerance': config_json['tolerance']
    }

    config['dexp'], config['times'] = read_dexp(get_path('dexp_file'))

    if config_json.get('pfact_file'):
        config['pfact'] = read_pfact(get_path('pfact_file'))

    return config
