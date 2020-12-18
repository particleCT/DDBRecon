import numpy as np
import sys
import matplotlib.pyplot as plt

A= np.load(sys.argv[1])['arr_0']
print(A.shape)
plt.imshow(A[:,:,100])
plt.show()

