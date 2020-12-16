from matplotlib import pyplot as plt
import numpy as np
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('inputPrefix',type=str)
parser.add_argument('figurePrefix',type=str)
parser.add_argument('msselTrajOutput', type=str)
parser.add_argument('effectivePopulationSize', type=int)
parser.add_argument('--ext',type=str,default='pdf')
args = parser.parse_args()

with open(args.msselTrajOutput) as f:
    lines = f.readlines()

### n = int(re.search("\d+", lines[2]).group(0))
data = np.genfromtxt(args.msselTrajOutput, delimiter=" ", skip_header=3)
### data = np.array(lines[-1].strip().split(" ")).astype(np.float).reshape((n, 2)
data[:,0] *= 4 * args.effectivePopulationSize

fig, ax = plt.subplots(1, 1, figsize=(20, 10))
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

ax.plot(data[:,0], data[:,1], 'r')

epochs = np.load(f'{args.inputPrefix}.epochs.npy')
freqs = np.load(f'{args.inputPrefix}.freqs.npy')
logpost = np.load(f'{args.inputPrefix}.post.npy')

im = ax.pcolormesh(epochs[:-1],freqs,np.exp(logpost)[:,:], shading="auto", vmax=np.max(np.exp(logpost)))
ax.axis((0,len(epochs[:-1]),0,1.0))
ax.set_ylabel('Allele frequency',fontsize=20)
ax.set_xlabel('Generations before present',fontsize=20)

cbar = fig.colorbar(im, ax=ax)
cbar.ax.set_ylabel('Posterior prob.\n\n',rotation=270, fontsize=20, labelpad=40)
cbar.ax.tick_params(labelsize=18)

fig.savefig('%s.%s'%(args.figurePrefix,args.ext),format=args.ext)
