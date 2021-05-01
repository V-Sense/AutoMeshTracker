import matplotlib.pyplot as plt
import numpy as np
import sys

fname = sys.argv[1]
e = np.loadtxt(fname)

plt.plot(e)
bottom, top = plt.ylim()
plt.ylim(0,top)
plt.ylabel('val')
plt.xlabel('neighbor')
outname = fname[0:len(fname)-4] + '.png'
plt.savefig(outname)