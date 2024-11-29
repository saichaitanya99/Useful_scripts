import numpy as np
import matplotlib.pyplot as plt
import libstempo as LT
import libstempo.plot as LP
import argparse
import os
from matplotlib import gridspec

parser = argparse.ArgumentParser(description='Residuals plot')

parser.add_argument('-parfile', type=str, required=True,help='par file')
parser.add_argument('-timfile', type=str, required=True,help='tim file')
parser.add_argument('-recfile', type=str, required=False,help='File with Recovered ToAs with columns MJD, TOA, TOA error. The columns should be separated by space and the comments should start with #')

args = parser.parse_args()

if not args.parfile or not args.timfile or not args.recfile:
    print("Required files not provided")
    exit()

parfile = args.parfile
timfile = args.timfile
recfile = args.recfile

if not os.path.isfile(parfile) or not os.path.isfile(timfile) or not os.path.isfile(recfile):
    print("Files not found")
    exit()  

# Load the pulsar and TOAs
    
psr = LT.tempopulsar(parfile,timfile, maxobs=100000)
#psr.fit(iters=2)

# Load the recovered TOAs
rec_toas = np.loadtxt(recfile, comments="#", usecols=(0,1,2), unpack=True)

# Plot the residuals
plt.figure(figsize=(10,6),)
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.1)

mask=np.argsort(psr.toas())
ax1 = plt.subplot(gs[0])
ax1.errorbar(psr.toas()[mask], psr.residuals()[mask], psr.toaerrs[mask]*1e-6, color = 'b', fmt='o', alpha=0.5, label='Original TOAs',markersize=2)
ax1.errorbar(rec_toas[0][mask], rec_toas[1][mask], rec_toas[2][mask], color = 'r', fmt='o', label='Recovered TOAs',markersize=2)
ax1.axhline(0, color='k', linestyle='--')
ax1.legend()
ax1.set_ylabel('Residuals (s)')

ax2 = plt.subplot(gs[1], sharex=ax1)
ax2.plot(psr.toas()[mask], rec_toas[1][mask]-psr.residuals()[mask], color = 'k', label='Recovered-Original',markersize=2,alpha=0.5)
ax2.axhline(0, color='k', linestyle='--')
ax2.set_xlabel('MJD')
ax2.set_ylabel('Recovered-Original (s)')
ax2.legend()
plt.show()







