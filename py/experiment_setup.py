""" experiment_setup.py
    
    Setup parameters in one place instead of having to edit them
    tediously in each file.
"""

EXPT_NAME = 'sim'

if EXPT_NAME == 'physionet':
    T = 60                      # s, total length of signal
    f0 = 256                    # Hz, original sampling rate
    T_window = 60               # s
    CF = 32                     # compression ratio
elif EXPT_NAME == 'sim':
    T = 600
    f0 = 4
    T_window = 60
    CF = 12
    
T0 = 1.0 / f0               # s
# samples per window, making sure it's a multiple of 4
N_window = (T_window * f0) - ((T_window * f0) % 4)
M = N_window // CF          # how many random samples to take per
                            # window
N0 = (T * f0) - ((T * f0) % (N_window // 2))
                            # total number of samples, throwing away
                            # incomplete windows


# Paths relative to project root, which is where we should run scripts
# from.
# TODO: do this more nicely.
original_samples_filename = 'data/{}/samples.csv'.format(EXPT_NAME)
phi_flags_filename \
        = 'data/{}/phi_flags_{}_{}.csv'.format(EXPT_NAME,
                                                N_window, M)
subsamples_filename \
        = 'data/{}/subsamples_{}_{}.csv'.format(EXPT_NAME,
                                                N_window, M)
reconstruction_filename \
        = 'data/{}/reconstructed_py.csv'.format(EXPT_NAME)
