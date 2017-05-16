import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mayavi.mlab import *

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

N = 6
x_lin = np.linspace(-2.25, 2.25, N)
y_lin = np.linspace(-2.25, 2.25, N)
xx, yy = np.meshgrid(x_lin, y_lin)
zz = np.zeros(xx.shape)
xx2 = np.zeros(xx.shape)
yy2 = np.zeros(xx.shape)
zz2 = np.zeros(xx.shape)
for i in range(N):
    for j in range(N):
        res = pca.transform(np.array([x_lin[i],y_lin[j],0]).reshape(1,-1))
        xx2[i,j] = res[0,0]
        yy2[i,j] = res[0,1]
        zz2[i,j] = res[0,2]

mpl = False
mayavi = False
if mayavi:
    s = surf(xx, yy, zz)
    show()
elif mpl:

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
        
    p1 = data[0,:]
    p2 = data[1,:]
    p3 = data[2,:]
    
    v1 = p3 - p1
    v2 = p2 - p1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = np.dot(cp, p3)

    zzp = (d - a * xx - b * yy) / c

    ax.scatter(xs = data[:,0], ys = data[:,1], zs = data[:,2], color = 'orange')
    ax.plot_surface(X = xx, Y = yy, Z = zzp, color = 'red')
    ax.plot_surface(X = xx, Y = yy, Z = zz, color = 'blue')
   
    plt.savefig('uncertainty.svg')
    plt.show()
else:

    camera = dict(
            up = dict(x = 0, y = 0, z = 1),
            center = dict(x = 0, y = 0, z = 0),
            eye = dict(x = -0.5, y = 2, z = 0.9)
    )


    points = go.Scatter3d(x = data[:,0], y = data[:, 1], z = data[:, 2], mode = 'markers', marker=dict(color = 'orange', size = 30))
    true_surface = go.Surface(x = xx, y = yy, z = zz, showscale = False, opacity = 0.75, colorscale = 'Greys')
    pca_surface = go.Surface(x = x_pca_plane, y = y_pca_plane, z = z_pca_plane, showscale = False, opacity = 0.75, colorscale = 'Greys')
    
    xaxis = go.XAxis(title = '', showaxeslabels = False, showticklabels = False)
    yaxis = go.YAxis(title = '', showaxeslabels = False, showticklabels = False)
    zaxis = go.ZAxis(title = '', showaxeslabels = False, showticklabels = False)
    scene = go.Scene(camera = camera, xaxis = xaxis, yaxis = yaxis, zaxis = zaxis)
    layout = go.Layout(scene=scene, width = 2000, height = 2000)

    fig = go.Figure(data = [points, true_surface, pca_surface], layout = layout)
    pyo.plot(fig)
