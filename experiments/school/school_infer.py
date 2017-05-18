#!/usr/bin/env python3
import numpy as np
import pystan

import sys
import math
import random
import pickle
from scipy.io import loadmat

def load_school_data(interval_length):
    # List of timestep sizes
    ALLDELTA = [interval_length]
    # Fontsize parameter for plotting
    FS = 28

    # Dataset
    NETDIR = 'primaryschool.csv'
    # ID file
    IDDIR = 'id_class_gender'

    # Import data
    ndat = open(NETDIR, 'r').readlines()
    ndat = np.asarray([np.asarray([elem for elem in line.strip().split('\t')[:-2]]+line.strip().split('\t')[-2:]) for line in ndat])

    # Get all time points
    all_t = np.asarray([int(elem) for elem in ndat[:,0]])
    # Zero time points
    all_t -= min(all_t)
    min_t = min(all_t)
    max_t = max(all_t)
    assert(min_t==0)

    # Get all node ID's
    all_n = list(set([elem for elem in ndat[:,1]]+[elem for elem in ndat[:,2]]))
    # Get all class ID's
    all_d = list(set([elem for elem in ndat[:,3]]+[elem for elem in ndat[:,4]]))
    nnodes = len(all_n)
    ndept = len(all_d)

    all_nodes = []
    for delta in ALLDELTA:
        
        curr_ndat = np.copy(ndat)
        curr_max_t = max_t
        curr_all_t = np.copy(all_t)
        
        # Coarse-grain time resolution
        # Units of time steps are seconds
        # Time steps are given in 20 second chunks
        curr_all_t //= 60*(delta+1)
        curr_max_t //= 60*(delta+1)
        
        max_min_diff = []

        all_A_d = []
        for ts in range(curr_max_t):

            curr_filt_ts = curr_all_t==ts

            ccurr_ndat = curr_ndat[curr_filt_ts]
            if len(ccurr_ndat.shape)<2:
                ccurr_ndat = np.asarray([ccurr_ndat])
            assert(len(ccurr_ndat.shape)==2)

            # Compute adjacency matrix of interactions that take place during time ts
            A_d = np.zeros((ndept+1,ndept+1), dtype=float)
            
            for row in ccurr_ndat:
                
                nrind1, nrind2 = all_n.index(row[1]), all_n.index(row[2])
                drind1, drind2 = all_d.index(row[3]), all_d.index(row[4])
                
                A_d[drind1][drind2] += 1.
                if drind1 != drind2:
                    A_d[drind2][drind1] = A_d[drind1][drind2]
                
            all_A_d.append(A_d.flatten())
        all_A_d = np.asarray(all_A_d).reshape((len(all_A_d),12,12)).astype(int)
        all_nodes.append(all_A_d)
    return all_nodes[0]

interval_length = 10
data_matrices = load_school_data(interval_length)

sm = pystan.StanModel(file = 'school_model.stan')

# Hyperparameters for rate matrix decomposition
n = 12
p = 2
d = int(n * p - p * (p+1)/2)

# Hyperparameters for HMM
n_states = 3
beta = [1/n_states for _ in range(n_states)]

# Parameters for Stan Inference
iter = 10
data = {'n' : n, 'p' : p, 'd' : d, 'K' : n_states, 'N' : data_matrices.shape[0], 'interactions' : data_matrices, 'beta' : beta}
chains = 1

fit = sm.sampling(data = data, iter = iter, chains = chains)

with open('fit.pkl', 'wb') as f:
    pickle.dump(fit.extract(), f)
