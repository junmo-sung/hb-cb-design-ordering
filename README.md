# Hybrid Beamformer Codebook Design and Ordering for Compressive mmWave Channel Estimation

The Matlab codebase in this repo can be used to reproduce the simulation results provided in the paper. The file 'main.m' is the top-level Matlab script that calls subfunctions provided in the 'functions' folder.

The pre-run data is also provided in the 'data' folder so that readers need not run the entire simulations to reproduce the results. In order to plot Figures 3, 4 and 5 in the paper, please run 'load_and_plot.m' in the folder.

To plot Fig. 2, run 'random\_permutation\_mut\_distribution.m' which performs Algorithm 1, stores the results, and plot the distribution. 