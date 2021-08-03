import glob
from ROOT import *
import numpy as np
from scipy.fftpack import fft, ifft, fftshift
import matplotlib.pyplot as plt
from root_numpy import hist2array
from scipy import ndimage
import sys
import matplotlib.pyplot as plt

#Not used in current reconstruction but maintained for reference
def get_fourier_filter(size, filter_name):
    n = np.concatenate((np.arange(1, size / 2 + 1, 2, dtype=int),
                        np.arange(size / 2 - 1, 0, -2, dtype=int)))
   
    f = np.zeros(size)
    f[0] = 0.25
    f[1::2] = -1 / (np.pi * n) ** 2
   
    # Computing the ramp filter from the fourier transform of its frequency domain representation lessens artifacts and removes a small bias
    fourier_filter = 2 * np.real(fft(f))         # ramp filter
    if filter_name == "ramp":
        pass
    elif filter_name == "shepp-logan":
        # Start from first element to avoid divide by zero
        omega = np.pi * fftmodule.fftfreq(size)[1:]
        fourier_filter[1:] *= np.sin(omega) / omega
    elif filter_name == "cosine":
        freq = np.linspace(0, np.pi, size, endpoint=False)
        cosine_filter = fftshift(np.sin(freq))
        fourier_filter *= cosine_filter
    elif filter_name == "hamming":
        fourier_filter *= fftmodule.fftshift(np.hamming(size))
    elif filter_name == "hann":
        fourier_filter *= fftmodule.fftshift(np.hanning(size))
    elif filter_name is None:
        fourier_filter[:] = 1
    return fourier_filter


i=0
FileList = sorted(glob.glob(sys.argv[1]+"*.root"))
NProj=len(FileList)
print(NProj)
repeat_list = []
for filename in FileList:


##Load the data
    f = TFile(filename)
    print(filename)
    proj      = f.Get("ddb_filtered")
    array     = hist2array(proj)
    NX,NY,NZ  = array.shape
    img_filtered   =  array
    if i ==0:
        img_final = np.zeros((NX,NY,NZ))

    ang = float(filename.split('_')[-1][0:-5])
    print(ang, img_filtered.max(), img_filtered.min())

    if img_filtered.max() < 10000 or img_filtered.min()>-10000:
        img_final = img_final + ndimage.rotate(img_filtered,ang,reshape = False,order=1,axes=(0,1))
    else:
        print("bad projection")
        repeat_list.append(ang)
    i+=1
FOV = 200.  
nbins = 320.
img_final = img_final*(nbins/FOV)*(np.pi/(2*NProj))
np.savez("recon_out.npz",img_final)

# This section results from a bug somewhere that means the ddb matrix is filled with nan as part of the fft.
# A second run through of the problem projections usually resolves this.
if len(repeat_list) !=0:
    print("please add the below list to teh secondrun.py file and run this script with the same parameters as before")
    print(repeat_list)
    
else:
    print("Complete")
