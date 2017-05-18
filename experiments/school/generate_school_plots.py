#!/usr/bin/env python3
import pickle
import numpy as np

from school_utils import load_school_data
import matplotlib.pyplot as plt

interval_length = 10
data_matrices = load_school_data(interval_length)

# Read the fit
with open('school_fit.pkl', 'rb') as f:
    fit = pickle.load(f)

angles = fit['theta']
transprobs = fit['phi']
y_star = fit['y_star']
Ws = fit['Ws']
Lambda = fit['Lambda']

# Plot the heatmap matrices
font = {'size'   : 30}
plt.rc('font', **font)

colorscheme = 'afmhot'
width = 20
height = 10

fig, axes = plt.subplots(nrows = 2, ncols = 3)
fig.set_size_inches(width, height)

im00 = axes[0,0].imshow(np.mean(np.exp(Ws[:,0,:,:]), axis = 0), vmin = 0, vmax = 200, cmap=plt.get_cmap(colorscheme))
axes[0,0].axis('off')
axes[0,0].grid(color='b', linestyle='-', linewidth=2)

im10 = axes[1,0].imshow(data_matrices[2], vmin = 0, vmax = 140, cmap=plt.get_cmap(colorscheme))
axes[1,0].axis('off')

im01 = axes[0,1].imshow(np.mean(np.exp(Ws[:,1,:,:]), axis = 0), vmin = 0, vmax = 200, cmap=plt.get_cmap(colorscheme))
axes[0,1].axis('off')

im11 = axes[1,1].imshow(data_matrices[23], vmin = 0, vmax = 140, cmap=plt.get_cmap(colorscheme))
axes[1,1].axis('off')

im02 = axes[0,2].imshow(np.mean(np.exp(Ws[:,2,:,:]), axis = 0), vmin = 0, vmax = 200, cmap=plt.get_cmap(colorscheme))
axes[0,2].axis('off')

im12 = axes[1,2].imshow(data_matrices[109], vmin = 0, vmax = 140, cmap=plt.get_cmap(colorscheme))
axes[1,2].axis('off')

fig.colorbar(im02, ax = axes[0,2])
fig.colorbar(im12, ax = axes[1,2])
plt.savefig('heatmap.svg')
