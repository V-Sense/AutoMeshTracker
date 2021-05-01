import matplotlib.pyplot as plt
import numpy as np
import sys
import os

fname = sys.argv[1]
e = np.loadtxt(fname)

plt.plot(e)
plt.ylabel('Feasibility Score')
plt.xlabel('Frame Number')
outname = fname[0:len(fname)-4] + '.png'
plt.savefig(outname)