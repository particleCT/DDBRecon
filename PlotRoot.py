import numpy as np
from root_numpy import hist2array
from ROOT import *
import sys
import matplotlib.pyplot as plt
f=TFile(sys.argv[1])
hist = f.Get("ddb_out")
hist = hist2array(hist)
hist = hist*np.pi/(2*360)
plt.imshow(hist[:,:,190])
plt.show()
plt.plot(hist[:,150,190])
plt.show()
print(f.ls())
