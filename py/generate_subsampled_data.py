#!/usr/bin/env python

""" generate_subsampled_data.py
    ---------------------------

"""

import numpy as np
from experiment_setup import *

original_data = np.loadtxt(original_samples_filename, delimiter=',')

# Take only as many samples as set in experiment_setup
original_data = original_data[:N0,:]
t0, x0 = original_data[:,0], original_data[:,1]

phi_flags = np.loadtxt(phi_flags_filename, delimiter=',')
phi_mat = np.zeros((M, N_window))
for i, v in enumerate(np.nonzero(phi_flags)[0]):
    phi_mat[i, v] = 1

# Generate subsampled windows.
t0_subsamples, x0_subsamples = np.empty((0, 1)), np.empty((0, 1))
for t in range(0, N0, N_window):
    t_window = t0[t:t+N_window]
    x_window = x0[t:t+N_window]
    subsamples_t = np.dot(phi_mat, t_window)
    subsamples = np.dot(phi_mat, x_window)

    t0_subsamples = np.vstack((t0_subsamples,
                               subsamples_t.reshape(-1, 1)))
    x0_subsamples = np.vstack((x0_subsamples,
                               subsamples.reshape(-1, 1)))
out_data_subsampled = np.hstack((t0_subsamples, x0_subsamples))
np.savetxt(subsamples_filename, out_data_subsampled, delimiter=',')
