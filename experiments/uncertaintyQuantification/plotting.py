import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import plotly.offline as pyo
import plotly.graph_objs as go

data = pd.read_csv('x.csv', index_col = 0)
data = data.values

pca = PCA(n_components = 3)
pca.fit(data)

pca_score = pca.explained_variance_ratio_
V = pca.components_
x_pca_axis, y_pca_axis, z_pca_axis = V.T * pca_score / pca_score.min()

x_pca_axis, y_pca_axis, z_pca_axis = 3 * V.T
x_pca_plane = np.r_[x_pca_axis[:2], -x_pca_axis[1::-1]]
y_pca_plane = np.r_[y_pca_axis[:2], -y_pca_axis[1::-1]]
z_pca_plane = np.r_[z_pca_axis[:2], -z_pca_axis[1::-1]]

x_pca_plane.shape = (2,2)
y_pca_plane.shape = (2,2)
z_pca_plane.shape = (2,2)

xx, yy = np.meshgrid(np.linspace(-2.25, 2.25, 10), np.linspace(-2.25, 2.25, 10))
zz = np.zeros(xx.shape)


mpl = False
if mpl:
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.scatter(xs = data[:,0], ys = data[:, 1], zs = data[:, 2])
    ax.plot_surface(X = xx, Y = yy, Z = zz, color = 'red')
    ax.plot_surface(X = x_pca_plane, Y = y_pca_plane, Z = z_pca_plane, color = 'blue')

    plt.show()
else:

    camera = dict(
            up = dict(x = 0, y = 0, z = 1),
            center = dict(x = 0, y = 0, z = 0),
            eye = dict(x = -0.5, y = 2, z = 0.5)
    )


    points = go.Scatter3d(x = data[:,0], y = data[:, 1], z = data[:, 2], mode = 'markers')
    true_surface = go.Surface(x = xx, y = yy, z = zz, showscale = False, opacity = 0.95)
    pca_surface = go.Surface(x = x_pca_plane, y = y_pca_plane, z = z_pca_plane, showscale = False, opacity = 0.95)

    layout = go.Layout(scene=dict(camera = camera))

    fig = go.Figure(data = [points, true_surface, pca_surface], layout = layout)
    pyo.plot(fig, image = 'svg')
