import argparse
import numpy as np
import sys
from matplotlib import pyplot as plt
from pathlib import Path

#Generate traj
def simulate_selected_forwards(p0,s,tOn,tOff,N=10000,eps=0.001):

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

def simulate_traj(p0,s,tOn,tOff,N):
    bwd = simulate_selected_backwards(p0,1e-8,N)
    fwd = simulate_selected_forwards(p0,s,tOn,tOff,N)
    traj = fwd[::-1]+bwd
    return traj

#Simulating aDNA samples from trajectory
def simulate_gls_from_traj(gens, nsamp, traj):
    '''
    Simulate genotype likelihoods of ancient samples

    Input:
        gens: number of generations back the ancient samples go
            (assuming they are uniformly distributed between [0,gens])
        nsamp: number of ancient samples

    Output:
        ancGLs: matrix of ancient genotype likelihoods, *in log space*;
                1st column is sampling time (in [0,gens]), columns 2-4 are
                genotype likelihoods of the genotypes.
        genos: vector of actual genotypes; each element is in {0,1,2}
    '''
    epochs = np.linspace(0,gens,gens)
    times = np.random.uniform(0,gens,size=nsamp)
    ancGLs = np.zeros((len(times),4))
    ancGLs[:,1:] = -np.inf
    genos = []
    for i,t in enumerate(np.sort(times)):
        if t >= gens - 1:
            ancGLs[i, 1]=0.0
            ancGLs[i, 0]=t
        else:
            geno = np.random.binomial(2,traj[np.digitize(t,epochs)+1])
            genos.append(geno)
            ancGLs[i,geno+1] = 0.0
            ancGLs[i,0] = t

    return ancGLs, genos


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
    argparser.add_argument("--output-file-path", type=str, default="ancientSamples.txt")
    argparser.add_argument("--ancient-samples-generation-gap", type=int, help="Number of generations back the ancient samples go", required=True)
    argparser.add_argument("--number-of-ancient-samples", type=int, help="Number of ancient samples", required=True)
    args = argparser.parse_args()

    p0 = args.initial_allele_freq
    s = args.selection_coefficient
    N = args.effective_population_size
    tOn = args.ton
    tOff = args.toff
    gens = args.ancient_samples_generation_gap
    nsamp = args.number_of_ancient_samples
    output_file_path = args.output_file_path

    ancGLs, genos= simulate_gls_from_traj(gens,nsamp,traj=simulate_traj(p0,s,tOn,tOff,N))
    np.savetxt(args.output_file_path, ancGLs)

if __name__ == '__main__':
    main()
