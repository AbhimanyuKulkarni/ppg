""" physionet_experiment_setup.py
    
    Setup parameters only for PhysioNet. Import this if you want those.
"""

EXPT_NAME = 'physionet'

T = 60                  # in s, total length of signal
f0 = 256                # in Hz, original sampling rate
T_window = 60           # in s
CF = 32                 # compression factor

T0 = 1.0 / f0           # s

# samples per window, making sure it's a multiple of 4
N_window = (T_window * f0) - ((T_window * f0) % 4)
M = N_window // CF      # how many random samples to take per window

N0 = (T * f0) - ((T * f0) % (N_window // 2))
                        # total number of samples, throwing away
                        # incomplete windows

# Paths relative to project root. So always run main.py from there only.
# TODO: do this more nicely.
original_samples_filename = 'data/{}/samples.csv'.format(EXPT_NAME)
phi_flags_filename \
        = 'data/{}/phi_flags_{}_{}.csv'.format(EXPT_NAME,
                                                N_window, M)
psiT_mat_filename \
        = 'data/{}/psiT_{}_{}.csv'.format(EXPT_NAME, N_window, M)
A_mat_filename \
        = 'data/{}/A_{}_{}.csv'.format(EXPT_NAME, N_window, M)
A_flip_mat_filename \
        = 'data/{}/A_flip_{}_{}.csv'.format(EXPT_NAME, N_window, M)
subsamples_filename \
        = 'data/{}/subsamples_{}_{}.csv'.format(EXPT_NAME,
                                                N_window, M)
reconstruction_filename \
        = 'data/{}/reconstructed_py.csv'.format(EXPT_NAME)
