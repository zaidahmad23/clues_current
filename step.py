import numpy as np
import sys
from matplotlib import pyplot as plt
from pathlib import Path
import argparse

# Simulation Pipeline
def simulate_selected_backwards(p0, s, N):
    delta = 1/(4*N)
    traj = [p0]
    a = s*N*2
    while traj[-1] != 1 and traj[-1] != 0:
        curr = traj[-1]
        nextFreq = np.random.normal(-a*curr*(1-curr)/np.tanh(a*curr)*delta + curr, np.sqrt(delta) * np.sqrt(curr*(1-curr)) )
        if nextFreq > 1:
            nextFreq = 1
        if nextFreq < 0:
            nextFreq = 0
        traj.append(nextFreq)
    return traj[1:]


def simulate_selected_forwards(p0, s, tOn, tOff, N, eps=0.001):
    delta = 1/(4*N)
    traj = [p0]

    for t in range(tOn,0,-1):
        if t > tOff:
            a = s*N*2
        else:
            a = 0.0

        curr = traj[-1]
        nextFreq = np.random.normal(a*curr*(1-curr)*delta + curr, np.sqrt(delta) * np.sqrt(curr*(1-curr)) )
        if nextFreq > 1:
            nextFreq = 1
        if nextFreq < 0:
            nextFreq = 0
        traj.append(nextFreq)

    if np.min([traj[-1],1-traj[-1]]) < eps:
        print('Warning: MAF less than %.4f'%(eps))

    return traj



def simulate_traj(p0, s, tOn, tOff, N):
    bwd = simulate_selected_backwards(p0, 1e-8, N)
    fwd = simulate_selected_forwards(p0, s, tOn, tOff, N)
    traj = fwd[::-1] + bwd
    print(traj[50:100:10])
    return traj


def save_mssel_input(N, output_file_path, traj):
    delta = float(1 / (4 * N))
    with open(output_file_path, "w") as f:
        f.writelines(["ntraj: 1\n", "npop: 1\n", f"n: {len(traj)}\n"])
        for (i, x) in enumerate(traj):
            f.write(f"{i * delta:.9f} {traj[i]:.6f}\n")


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "-p0",
        "--initial-allele-freq",
        type=float,
        required=True,
        help="Present-day allele frequency.",
    )
    argparser.add_argument(
        "-s",
        "--selection-coefficient",
        type=float,
        required=True,
        help="Selection coefficient.",
    )
    argparser.add_argument(
        "-n",
        "--effective-population-size",
        type=float,
        required=True,
        help="Effective population size.",
    )
    argparser.add_argument(
        "--ton",
        type=int,
        required=True,
        help="Time before present that selection starts.",
    )
    argparser.add_argument(
        "--toff",
        type=int,
        required=True,
        help="Time before present that selection ends.",
    )
    argparser.add_argument("--output-file-path", type=str, default="mssel.traj")
    args = argparser.parse_args()

    p0 = args.initial_allele_freq
    s = args.selection_coefficient
    N = args.effective_population_size
    tOn = args.ton
    tOff = args.toff

    output_file_path = args.output_file_path

    save_mssel_input(
        N, traj=simulate_traj(p0, s, tOn, tOff, N), output_file_path=output_file_path
    )


if __name__ == "__main__":
    main()
