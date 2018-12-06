#!/usr/bin/env python

import sys
EXPT_NAME = sys.argv[1]

if EXPT_NAME == 'physionet':
    from physionet_experiment_setup import *
elif EXPT_NAME == 'sim':
    from sim_experiment_setup import *

import numpy as np
import matplotlib.pyplot as plt

# from sklearn import linear_model
from lasso import cd_lasso
from scipy.signal import lombscargle

def half_flip(phi):
    M, N = phi.shape
    phi_flipped = np.zeros((M,N))
    phi_flipped[:M//2,:N//2] = phi[M//2:,N//2:]
    phi_flipped[M//2:,N//2:] = phi[:M//2,:N//2]
    return phi_flipped

def zcr(X, t):
    X_centred = X - np.mean(X)
    N, = X.shape
    Xl = X_centred[:N-1]
    Xr = X_centred[1:]
    Xprod = Xl * Xr
    Xsum = np.abs(Xl) + np.abs(Xr)
    nzc = sum(Xprod <= 0)
    n_bothzero = sum(Xsum == 0)
    nzc -= n_bothzero
    return nzc / (2 * (t[-1] - t[0]))

def damp(x, p):
    y = np.copy(x)
    for i in range(1, y.size):
        y[i] = p * y[i] + (1 - p) * y[i-1]
    return y

# Simulation to make sure we can run the compressive sampling algorithm
# on some data at least.
original_data = np.loadtxt(original_samples_filename, delimiter=',')

# Take only as many samples as set in experiment_setup
original_data = original_data[:N0,:]

subsampled_data = np.loadtxt(subsamples_filename, delimiter=',')
phi_flags = np.loadtxt(phi_flags_filename, delimiter=',')
phi = np.zeros((M, N_window))
for i, v in enumerate(np.nonzero(phi_flags)[0]):
    phi[i, v] = 1
phi_flipped = half_flip(phi)

t0 = original_data[:,0]
X0 = original_data[:,1]

ts = subsampled_data[:,0]
ys = subsampled_data[:,1]

# Estimate period directly.
freqs = np.arange(20, 300, 3)   # in bpm
freqs_Hz = freqs / 60.0         # in Hz
pgram_direct = lombscargle(ts, ys, freqs_Hz, precenter=True)
freq_max_direct = freqs[np.argmax(pgram_direct)]
print('f_max (direct): {}'.format(freq_max_direct))

psi = np.empty((N_window, N_window))
n, k = np.meshgrid(np.arange(N_window), np.arange(N_window))
psi = np.cos((np.pi / N_window) * (n + 0.5) * k)
psi[0,:] *= (1.0 / np.sqrt(2))
psi *= np.sqrt(2.0 / N_window)

A = np.dot(phi, psi.T)
A_flipped = np.dot(phi_flipped, psi.T)
Xr = np.zeros(X0.shape)

# Filling in the 1st quarter-window
Y = ys[:M]
s = cd_lasso(Y, A)
xr = np.dot(psi.T, s)
Xr[:N_window//4] = xr[:N_window//4]

for i in range(0, ys.size - M + 1, M//2):
    Y = ys[i:i+M]
    if i % M == 0:
        s = cd_lasso(Y, A)
    else:
        s = cd_lasso(Y, A_flipped)
    xr = np.dot(psi.T, s)

    k = i * (N_window // M)
    Xr[k+N_window//4:k+3*N_window//4] = xr[N_window//4:3*N_window//4]

    # Estimate period from reconstruction.
    freq_reconstruct = zcr(np.diff(Xr[k+N_window//4:k+3*N_window//4]),
                           t0[k+N_window//4+1:k+3*N_window//4]) * 60
    # freq_reconstruct = zcr(Xr[k+N_window//4:k+3*N_window//4],
    #                        t0[k+N_window//4:k+3*N_window//4])
    # pgram_window = lombscargle(t0[k+N_window//4:k+3*N_window//4],
    #                            Xr[k+N_window//4:k+3*N_window//4],
    #                            freqs_Hz,
    #                            precenter=True)
    # freq_reconstruct = freqs[np.argmax(pgram_window)]
    print('Frequency: {} for [{}, {}]'.format(freq_reconstruct,
                                              t0[k+N_window//4],
                                              t0[k+3*N_window//4]))

# Filling in the last quarter-window
Y = ys[-M:]
if ys.size % M == 0:
    s = cd_lasso(Y, A)
else:
    s = cd_lasso(Y, A_flipped)
xr = np.dot(psi.T, s)
Xr[-N_window//4:] = xr[-N_window//4:]

plt.plot(t0, X0, label='original')
plt.plot(t0, Xr, label='reconstructed')
plt.legend()
plt.show()

print('Corrcoef: ', np.corrcoef(X0, Xr))
print('Saving to {}...'.format(reconstruction_filename))
out_data = np.hstack((t0.reshape(-1, 1), Xr.reshape(-1, 1)))
np.savetxt(reconstruction_filename, out_data, fmt='%f', delimiter=',')
