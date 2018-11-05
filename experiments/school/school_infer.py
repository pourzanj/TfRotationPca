#!/usr/bin/env python3
import numpy as np
# import pystan
# import pickle
from school_utils import load_school_data

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

with open('school_fit.pkl', 'wb') as f:
    pickle.dump(fit.extract(), f)
