#!/usr/bin/env python3

import pickle
import numpy as np
import plotly.offline as pyo
import plotly.graph_objs as go

with open('fit.pkl', 'rb') as f:
    fit = pickle.load(f)

angles = fit['theta']
transprobs = fit['phi']
y_star = fit['y_star']
Ws = fit['Ws']
Lambda = fit['Lambda']
