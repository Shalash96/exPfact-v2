"""
DynamX Cluster File to exPfact Format Converter
-----------------------------------------------

This script processes a DynamX cluster CSV file and generates:
  - An assignment file (.list) representing peptide coverage
  - A normalized uptake file (.dexp) compatible with exPfact

Normalization options:
  1. Theoretical maximum uptake (uses MaxUptake column x D2O %)
  2. Fully deuterated sample as a reference (last timepoint per peptide)

Author: Original by E. Paci group, upgraded by Mahmoud Shalash
"""

import pandas as pd
import numpy as np
import sys
import os

def main():
    if len(sys.argv) != 2:
        print("Usage: process_DnXcluster.py <filename.csv>")
        sys.exit(1)

    infile = sys.argv[1]
    if not os.path.isfile(infile):
        print(f"Error: file {infile} not found.")
        sys.exit(1)

    alldata = pd.read_csv(infile)
    print(f"Converting DynamX cluster '{infile}' into exPfact input files.\n")

    # Mass correction
    alldata['Dt'] = alldata['Center'] * alldata['z'] - alldata['z']
    alldata['Time'] = round(alldata['Exposure'], 2)
    alldata['CorrDt'] = 0

    # Identify states and timepoints
    states = sorted(set(alldata['State']))
    times = np.sort(alldata['Exposure'].unique())

    print("Experimental conditions (States):")
    for i, s in enumerate(states, 1):
        print(f"  {i}) {s}")

    print("\nLabeling times (minutes):")
    for i, t in enumerate(times, 1):
        print(f"  {i}) {t:5.1f}")

    print("\nNormalization Method:")
    print("  1) Theoretical max uptake (uses MaxUptake × D2O%)")
    print("  2) Fully deuterated sample (uses final exposure time as max)")
    norm_method = input("Select [1/2]: ").strip()

    if norm_method == "1":
        d2o_perc = float(input("Enter D2O percentage [0-100]: ").strip())
    elif norm_method == "2":
        back_exs = []
        peptides = set(alldata['Sequence'])
        final_time = times[-1]
        for pep in peptides:
            sub = alldata[(alldata["Sequence"] == pep) & (alldata["Exposure"] == final_time)]
            if not sub.empty:
                fd = np.average(sub["Dt"])
                mhp = np.average(sub["MHP"])
                maxu = np.average(sub["MaxUptake"])
                back_ex = (fd - mhp) / maxu
                back_exs.append(back_ex)
        back_ex_perc = np.mean(back_exs) * 100
        print(f"\nAverage back-exchange plateau: {back_ex_perc:.2f}%")
    else:
        print("Invalid normalization method selected.")
        sys.exit(1)

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
        ass_path = infile.replace(".csv", f"_{state}.list")
        with open(ass_path, 'w') as f:
            for i, (start, end) in enumerate(peptides):
                f.write(f"{i+1} {start} {end} {sequences[i]}\n")

        times = np.sort(data['Time'].unique())
        dexp = np.zeros((len(peptides)+1, len(times)))
        dexp[0] = times / 60  # Convert minutes to hours

        for i, (start, end) in enumerate(peptides, 1):
            sub = data[(data['Start'] == start) & (data['End'] == end)]
            control = sub[sub['Time'] == times[0]]['Dt'].mean()

            if norm_method == "1":
                maxUptake = d2o_perc / 100 * sub['MaxUptake'].iloc[0]
            else:
                fd = sub[sub['Time'] == times[-1]]['Dt']
                if not fd.empty:
                    maxUptake = fd.mean()
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

        dexp_out = infile.replace(".csv", f"_{state}.dexp")
        np.savetxt(dexp_out, dexp.T[:-1] if norm_method == "2" else dexp.T, fmt="%5.5f")

        print(f"  ➤ Saved: {ass_path}, {dexp_out}")


if __name__ == '__main__':
    main()
