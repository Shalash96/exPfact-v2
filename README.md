# ExPfact Suite: Advanced HDX-MS Protection Factor Analysis

**Current Version: 2.0 (Python 3.13 Upgrade)**

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
<!-- Optional: Add build status, DOI, etc. badges here if applicable -->
<!-- e.g., [![DOI](https://zenodo.org/badge/DOI/your-zenodo-doi.svg)](https://doi.org/your-zenodo-doi) -->

**Note:** This repository contains a significantly modernized version of the original [*ExPfact*](https://github.com/pacilab/exPfact) software suite for Hydrogen-Deuterium Exchange Mass Spectrometry (HDX-MS) data analysis. It has been fully ported to Python 3.13 and extensively refactored for improved clarity, maintainability, and usability. While the core scientific algorithms aim to replicate the original's functionality, this version is under active development. Rigorous testing and validation against established experimental benchmarks and the original ExPfact results are ongoing.

---

## Overview

**ExPfact** is a computational framework designed to estimate protection factors (PFs) at single-residue resolution from Hydrogen-Deuterium Exchange Mass Spectrometry (HDX-MS) data. A key strength of ExPfact is its ability to handle sparse and underdetermined datasets, which are common in HDX-MS experiments.

The suite applies statistical modeling and optimization techniques to infer one or multiple plausible P-factor profiles (typically as ln(P)) consistent with experimental deuterium uptake kinetics measured at the peptide level. This allows for detailed insights into protein structure, dynamics, and interactions.

This upgraded version includes a PyQt6-based Graphical User Interface (GUI) to streamline the workflow across its various analytical modules.

---

## Citation

If you use the ExPfact suite or its underlying methodologies in your research, please cite the relevant original publications:

1.  Skinner, S. P., Radou, G., Tuma, R., Houwing-Duistermaat, J. J., & Paci, E. (2019). Estimating Constraints for Protection Factors from HDX-MS Data. *Biophysical Journal*, 116(7), 1194–1203.
    [DOI: 10.1016/j.bpj.2019.02.024](https://doi.org/10.1016/j.bpj.2019.02.024)

2.  Stofella, M., Skinner, S. P., Sobott, F., Houwing-Duistermaat, J., & Paci, E. (2022). High-Resolution Hydrogen–Deuterium Protection Factors from Sparse Mass Spectrometry Data Validated by Nuclear Magnetic Resonance Measurements. *Journal of the American Society for Mass Spectrometry*, 33(5), 813–822.
    [DOI: 10.1021/jasms.2c00005](https://pubs.acs.org/doi/full/10.1021/jasms.2c00005)

---

## Key Features & Updates in This Version

This modernized ExPfact suite offers:

*   **Python 3.13 Compatibility:** Fully migrated from earlier Python versions.
*   **Comprehensive Code Refactoring:**
    *   **Type Hinting:** Enhanced code clarity and maintainability with Python type hints.
    *   **Docstrings:** NumPy/reStructuredText style documentation for all modules and functions.
    *   **Improved Modularity:** Better separation of concerns across different Python scripts.
    *   **Enhanced Error Handling:** More robust error checking and informative logging.
    *   **Modern Libraries:** Compatibility with current versions of core scientific Python libraries (NumPy, SciPy, Pandas, MDAnalysis, pyOpenMS, BioPython).
*   **PyQt6 Graphical User Interface (GUI):**
    *   A user-friendly interface to access and run all modules of the ExPfact suite.
    *   Streamlined workflow for typical HDX-MS analysis tasks.
*   **Simplified Installation:** Primarily through `conda` using a provided `environment.yml` file.
*   **Core Functionality Preserved:** Includes modules for:
    *   P-factor fitting (`exPfact.py`)
    *   Calculating protection factors from MD simulations (`MD2Pfact.py`)
    *   Predicting deuterium uptake from P-factors (`pfact2dpred.py`)
    *   Calculating theoretical isotopic envelopes (`Hisotope.py`, `isotopic_envelope.py`)
    *   Data processing and utility scripts (e.g., `process_DnXcluster.py`, `clustering.py`, `descriptive_statistics.py`, `cross_validation.py`)

### Planned Enhancements / Future Work

*   **Advanced Visualization:** Integration of tools for visualizing deuterium uptake mapped onto protein structures (PDB files).
*   **Custom Plotting:** More flexible and customizable plotting options within the GUI and scripts.
*   **Reproducible Deployment:** Options for Docker-based deployment to ensure consistent environments.
*   **Tutorial Datasets:** Expanded set of example and tutorial datasets to facilitate learning and testing.
*   **Expanded Data Import:** Support for a wider range of common HDX-MS data formats and other relevant scientific file types.

---

## Installation

Follow these steps to set up the ExPfact suite:

### Prerequisites

*   [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda) installed.
*   A C/C++ compiler compatible with your Python version (for compiling Cython extensions).
    *   **Linux:** Typically `gcc` (install via your package manager, e.g., `sudo apt-get install build-essential`).
    *   **macOS:** Xcode Command Line Tools (`xcode-select --install`).
    *   **Windows:** Microsoft C++ Build Tools (available with Visual Studio Installer).

### Steps

1.  **Clone the Repository (if applicable):**
    ```bash
    git clone https://github.com/Shalash96/exPfact-v2.git
    cd exPfact-v2 # Or your repository name
    ```

2.  **Create and Activate the Conda Environment:**
    We provide an `env.yml` file to create an environment with all necessary dependencies.
    ```bash
    conda env create -f env.yml
    conda activate exPfact-v2 # The environment name is defined in env.yml
    ```
    *Alternatively, to create the environment manually (ensure all dependencies from `env.yml` are included):*
    ```bash
    # conda create -n exPfact-v2 python=3.13 numpy scipy pandas matplotlib mdanalysis pyopenms biopython pyqt6 cython
    # conda activate exPfact-v2
    # conda install -c conda-forge r-base # If not already in env.yml
    ```

3.  **Compile Cython Modules:**
    The core `calc_dpred` module is written in Cython for performance. Navigate to the directory containing the Cython source files which is `python/` and compile:
    ```bash
    cd python/
    python setup_calc_dpred.py build_ext --inplace
    ```
    *(Ensure `setup_calc_dpred.py` correctly points to `calc_dpred.pyx`)*

4.  **Install R and `mclust` Package (for Clustering):**
    The clustering functionality relies on R and the `mclust` package.
    *   If `r-base` was not included in your `env.yml`, install it via conda:
        ```bash
        conda install -c conda-forge r-base
        ```
    *   Then, start an R session from your activated conda environment and install `mclust`:
        ```R
        # Inside an R console
        install.packages("mclust")
        ```
        (You might be prompted to choose a CRAN mirror.)

---

## Getting Started

1.  **Activate the Conda Environment:**
    ```bash
    conda activate exPfact-v2
    ```

2.  **Launch the Graphical User Interface (GUI):**
    Navigate to the directory containing `expfact_gui.py` which is `guiApp/` (or adjust path) and run:
    ```bash
    python expfact_gui.py
    ```
    The GUI provides access to all modules and their parameters. It's recommended to start here.


3.  **Tutorial / Test Data:**
    Explore the `testing/` directory for sample datasets and walkthroughs to familiarize yourself with the workflow and expected file formats.

---

## Current Status & Contribution

This modernized version of ExPfact is currently in an **active development / beta stage**. While core functionalities are implemented, comprehensive validation against a wide range of experimental datasets and comparison with results from the original ExPfact version are ongoing.


**Key areas for future development and validation include:**
*   Rigorous benchmarking against diverse HDX-MS datasets.
*   Comparison with outputs from the original ExPfact Fortran/Python2 versions.
*   Implementation of the "Features in Development" listed above.

---
