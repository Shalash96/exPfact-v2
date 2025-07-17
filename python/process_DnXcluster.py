"""
DynamX Cluster File to exPfact Format Converter
-----------------------------------------------

This script processes a DynamX cluster CSV file and generates:
  - An assignment file (.list) representing peptide coverage
  - A normalized uptake file (.dexp) compatible with exPfact

Normalization options:
  1. Theoretical maximum uptake (uses MaxUptake column × D2O %)
  2. Fully deuterated sample as a reference (last timepoint per peptide)

Author: Original by E. Paci group, upgraded by Mahmoud Shalash
"""

import argparse
import pandas as pd
import numpy as np
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="Convert DynamX cluster file to exPfact format.")
    parser.add_argument("infile", help="Path to input CSV file (DynamX cluster export)")
    parser.add_argument("--norm_method", required=True, choices=["1", "2"],
                        help="Normalization method: 1 = theoretical max uptake; 2 = fully deuterated sample")
    parser.add_argument("--d2o_perc", type=float, default=None,
                        help="D2O percentage (only required for norm_method=1)")
    args = parser.parse_args()

    infile = args.infile
    norm_method = args.norm_method
    d2o_perc = args.d2o_perc

    if not os.path.isfile(infile):
        print(f"Error: file '{infile}' not found.")
        sys.exit(1)

    alldata = pd.read_csv(infile)
    print(f"Converting DynamX cluster '{infile}' into exPfact input files.\n")

    # Mass correction
    alldata['Dt'] = alldata['Center'] * alldata['z'] - alldata['z']
    alldata['Time'] = round(alldata['Exposure'], 2)
    alldata['CorrDt'] = 0

    # Extract clean state labels and times
    states = sorted(set(str(s) for s in alldata['State'] if pd.notna(s)))
    times = alldata['Time'].dropna().unique()
    times = np.sort(times)

    print("Experimental conditions (States):")
    for i, s in enumerate(states, 1):
        print(f"  {i}) {s}")

    print("\nLabeling times (minutes):")
    for i, t in enumerate(times, 1):
        print(f"  {i}) {t:5.1f}")

    if norm_method == "1":
        if d2o_perc is None:
            print("Error: --d2o_perc must be specified when using norm_method 1.")
            sys.exit(1)
    elif norm_method == "2":
        back_exs = []
        peptides = set(alldata['Sequence'])
        final_time = times[-1]
        for pep in peptides:
            sub = alldata[(alldata["Sequence"] == pep) & (alldata["Exposure"] == final_time)]
            if not sub.empty and sub["MHP"].notna().all() and sub["MaxUptake"].notna().all():
                fd = sub["Dt"].dropna().mean()
                mhp = sub["MHP"].dropna().mean()
                maxu = sub["MaxUptake"].dropna().mean()
                if maxu != 0:
                    back_ex = (fd - mhp) / maxu
                    back_exs.append(back_ex)
        if back_exs:
            back_ex_perc = np.mean(back_exs) * 100
            print(f"\nAverage back-exchange plateau: {back_ex_perc:.2f}%")
        else:
            print("\nWarning: No valid peptides found for estimating back-exchange.")
            back_ex_perc = 0

    for state in states:
        data = alldata[alldata["State"] == state].reset_index()

        peptides = []
        sequences = []
        for i in range(len(data)):
            pair = (data['Start'][i], data['End'][i])
            if pair not in peptides:
                peptides.append(pair)
                sequences.append(data['Sequence'][i])

        print(f"\nState: {state} — {len(peptides)} peptides")

        # Assignment file
        basename = os.path.splitext(os.path.basename(infile))[0]
        ass_path = f"{basename}_{state}.list"        
        with open(ass_path, 'w') as f:
            for i, (start, end) in enumerate(peptides):
                f.write(f"{i+1} {start} {end} {sequences[i]}\n")

        dexp = np.zeros((len(peptides)+1, len(times)))
        dexp[0] = times / 60  # Convert minutes to hours

        for i, (start, end) in enumerate(peptides, 1):
            sub = data[(data['Start'] == start) & (data['End'] == end)]
            control = sub[sub['Time'] == times[0]]['Dt'].mean()

            if norm_method == "1":
                maxUptake = d2o_perc / 100 * sub['MaxUptake'].iloc[0]
            else:
                fd = sub[sub['Time'] == times[-1]]['Dt'].max(skipna=True)
                if not np.isnan(fd):
                    maxUptake = fd
                else:
                    maxUptake = back_ex_perc / 100 * sub['MaxUptake'].iloc[0]

            for j, t in enumerate(times):
                data_t = sub[sub['Time'] == t]
                if not data_t.empty:
                    if norm_method == "1":
                        dt = (data_t['Dt'].mean() - control) / maxUptake
                    else:
                        dt = (data_t['Dt'] - control) / (maxUptake - control)
                        dt = dt.mean()
                else:
                    print(f"  Missing time {t} min for peptide {sequences[i-1]} in {state}")
                    dt = np.nan

                dt = min(max(dt, 0), 1) if not np.isnan(dt) else 0
                dexp[i][j] = dt

        dexp_out = f"{basename}_{state}.dexp"
        np.savetxt(dexp_out, dexp.T[:-1] if norm_method == "2" else dexp.T, fmt="%5.5f")

        print(f"  ➤ Saved: {ass_path}, {dexp_out}")

if __name__ == '__main__':
    main()
