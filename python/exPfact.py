"""
Main ExPfact script: fits protection factors to HDX-MS experimental data.

This script performs optimization of protection factors using optional
harmonic restraints and random initialization. Supports multiprocessing
for repeated minimizations.

Required:
    --temp : Temperature in Kelvin
    --pH   : pH of the system
    --dexp : Experimental deuterium uptake file
    --ass  : Assignment file

Optional:
    --base, --out, --pfact, --rand, --tol, --harm, --rep, --ncores
"""

import os
import sys
import time
import argparse
from multiprocessing import Pool

from calc_dpred import calculate_dpred
from calculate import cost_function, do_random_search, fit_pfact
from kint import calculate_kint_for_sequence
from read import (
    read_assignments,
    read_configuration,
    read_dexp,
    read_pfact,
    read_seq,
)
from write import write_diff, write_dpred, write_pfact
from logger import log


def run(base_dir, dexp, assignments, pfact, random_steps, time_points,
        harmonic_term, output_file, tolerance, weights, pH, temperature,
        seq, res1, resn):
    """
    Runs a single protection factor optimization process.
    """
    log.info(f"Running ExPfact, output name: {output_file}")

    assignment_set = {x for ass in assignments for x in range(int(ass[1]), int(ass[2]) + 1)}
    pfactor_filter = {x for ass in assignments for x in range(int(ass[1] + 1), int(ass[2]) + 1)}
    for ass in assignments:
        if ass[1] < min(pfactor_filter):
            pfactor_filter.add(ass[1])

    kint, prolines = calculate_kint_for_sequence(res1, resn, seq, temperature, pH)

    if not pfact:
        if random_steps:
            rand_output = do_random_search(
                kint, random_steps, pfactor_filter, dexp, time_points,
                assignments, harmonic_term, prolines, weights, seed=None
            )
            init_array = rand_output[min(rand_output.keys())]
        else:
            init_array = [
                1 if ii not in prolines or ii == 0 or ii + 1 in pfactor_filter else -1
                for ii in range(len(seq))
            ]
    else:
        init_array = read_pfact(pfact)

    bounds = [
        (0.00001, 30) if x >= 0 else (-1, -1) if x == -1 else (0, 0)
        for x in init_array
    ]

    pfit = fit_pfact(
        init_array, dexp, time_points, assignments,
        harmonic_term, kint, bounds, tolerance, weights
    )

    write_pfact(pfit.x, output_file)
    dpred = calculate_dpred(pfit.x, time_points, kint, assignments)
    write_dpred(output_file, dpred, time_points)
    write_diff(output_file, dpred, dexp)

    score_harm = cost_function(pfit.x, dexp, time_points, assignments, harmonic_term, kint, weights)
    score_no_harm = cost_function(pfit.x, dexp, time_points, assignments, 0.0, kint, weights)

    print(f"Final cost function (with harm): {score_harm}")
    print(f"Final cost function (no harm): {score_no_harm}")


def main(argv):
    log.info("Running exPfact.py")

    parser = argparse.ArgumentParser()
    parser.add_argument("--base")
    parser.add_argument("--dexp")
    parser.add_argument("--ass")
    parser.add_argument("--pfact")
    parser.add_argument("--weights")
    parser.add_argument("--out")
    parser.add_argument("--predict")
    parser.add_argument("--times")
    parser.add_argument("--tol")
    parser.add_argument("--harm")
    parser.add_argument("--rand")
    parser.add_argument("--temp")
    parser.add_argument("--pH")
    parser.add_argument("--seq")
    parser.add_argument("--rep")
    parser.add_argument("--ncores")

    if len(sys.argv) > 1 and sys.argv[1].endswith(".json"):
        config = read_configuration(sys.argv[1])
    else:
        opts = parser.parse_args()
        config = {}

        config["base"] = opts.base if opts.base else os.getcwd()
        if opts.dexp:
            config["dexp"], config["times"] = read_dexp(opts.dexp)
        if opts.ass:
            config["assignments"] = opts.ass
        if opts.temp:
            config["temperature"] = float(opts.temp)
        if opts.pH:
            config["pH"] = float(opts.pH)
        if opts.seq:
            config["sequence"] = read_seq(opts.seq)
            config["res1"] = 1
            config["resn"] = len(config["sequence"])
        config["pfact"] = opts.pfact if opts.pfact else None
        config["output"] = opts.out if opts.out else None
        config["do_random_search"] = bool(opts.rand)
        config["random_search_steps"] = int(opts.rand) if opts.rand else None
        config["tolerance"] = float(opts.tol) if opts.tol else None
        config["harmonic_factor"] = float(opts.harm) if opts.harm else 0
        config["weights"] = read_dexp(opts.weights)[0] if opts.weights else None
        n_rep = int(opts.rep) if opts.rep else 1
        n_cores = int(opts.ncores) if opts.ncores else 1

    assignments = read_assignments(config["assignments"])

    tic = time.time()
    log.info("ExPfact starts")

    args = []
    for i in range(1, n_rep + 1):
        out_name = f"{config['output']}{i}" if n_rep > 1 else config["output"]
        args.append((
            config["base"],
            config["dexp"],
            assignments,
            config["pfact"],
            config["random_search_steps"],
            config["times"],
            config["harmonic_factor"],
            out_name,
            config["tolerance"],
            config["weights"],
            config["pH"],
            config["temperature"],
            config["sequence"],
            config["res1"],
            config["resn"]
        ))

    with Pool(n_cores) as pool:
        pool.starmap(run, args)

    log.info(f"ExPfact ends, total time: {time.time() - tic:.5f} seconds")


if __name__ == "__main__":
    try:
        sys.argv[1]
    except IndexError:
        print(__doc__)
        log.error("Missing arguments for ExPfact.py")
        sys.exit(1)
    main(sys.argv[1:])
