#!/usr/bin/env python

from experiment_setup import *
import numpy as np

phi_flags = np.zeros((N_window,), dtype=np.int)
chosen = np.random.choice(N_window, size=M, replace=False)
phi_flags[chosen] = 1
np.savetxt(phi_flags_filename, phi_flags, fmt='%d', delimiter='\n')
