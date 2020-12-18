from ROOT import *
import numpy as np
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt
from root_numpy import hist2array
from scipy import ndimage
import sys

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
        cosine_filter = fftmodule.fftshift(np.sin(freq))
        fourier_filter *= cosine_filter
    elif filter_name == "hamming":
        fourier_filter *= fftmodule.fftshift(np.hamming(size))
    elif filter_name == "hann":
        fourier_filter *= fftmodule.fftshift(np.hanning(size))
    elif filter_name is None:
        fourier_filter[:] = 1
    return fourier_filter

##Load the data
f = TFile(sys.argv[1])
proj      = f.Get("ddb_0")
proj      = proj.ProjectionXYZ("","")
array     = hist2array(proj)
array     = ndimage.gaussian_filter(array,sigma=(3,3,3))
NX,NY,NZ  = array.shape
theta     = 360
img_shape = array.shape[1]
projection_size_padded = max(64, int(2 ** np.ceil(np.log2(2 * img_shape)))) ## Pad in the Y direction, this is where we cumulate the Radon DDB. Particles travel along X, and Rotation is around Z
pad_width = ((0,0),(0, projection_size_padded - img_shape), (0, 0))
img       = np.pad(array, pad_width, mode='constant', constant_values=0)
fourier_filter = get_fourier_filter(projection_size_padded, "ramp")
fourier_filter = fourier_filter[np.newaxis,:, np.newaxis] ## Adjust the axis to multiply corectly
projection     = fft(img,axis=1)*fourier_filter ## Fourier transform along the Y axis
img_filtered   = np.real(ifft(projection,axis=1)[:,:img_shape, :]) #Inverse fourier also alogn the Y axis
img_final      = np.zeros((NX,NY,NZ))

Arc   = 360
NProj = 24
Step  = int(Arc/NProj)
for ang in range(0,360,Step):
    print(ang)
    img_final = img_final + ndimage.rotate(img_filtered, ang,reshape=False,order=1, axes=(0,1))
img_final = img_final * np.pi/(2*NProj)
np.savez("recon.npz",img_final)
