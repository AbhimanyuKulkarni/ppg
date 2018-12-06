#!/usr/bin/env python

from experiment_setup import *
import numpy as np

def half_flip(phi):
    M, N = phi.shape
    phi_flipped = np.zeros((M,N))
    phi_flipped[:M//2,:N//2] = phi[M//2:,N//2:]
    phi_flipped[M//2:,N//2:] = phi[:M//2,:N//2]
    return phi_flipped

phi_flags = np.zeros((N_window,), dtype=np.int)
chosen = np.random.choice(N_window, size=M, replace=False)
phi_flags[chosen] = 1
np.savetxt(phi_flags_filename, phi_flags, fmt='%d', delimiter='\n')

# Dump psi and A too
phi_mat = np.zeros((M, N_window))
for i, v in enumerate(np.nonzero(phi_flags)[0]):
    phi_mat[i, v] = 1
phi_flipped = half_flip(phi_mat)
psi = np.empty((N_window, N_window))
n, k = np.meshgrid(np.arange(N_window), np.arange(N_window))
psi = np.cos((np.pi / N_window) * (n + 0.5) * k)
psi[0,:] *= (1.0 / np.sqrt(2))
psi *= np.sqrt(2.0 / N_window)

A = np.dot(phi_mat, psi.T)
A_flipped = np.dot(phi_flipped, psi.T)

np.savetxt(psiT_mat_filename, psi.T, delimiter=',')
np.savetxt(A_mat_filename, A, delimiter=',')
np.savetxt(A_flip_mat_filename, A_flipped, delimiter=',')
