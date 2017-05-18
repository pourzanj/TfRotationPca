import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

data = pd.read_csv('x.csv', index_col = 0)
data = data.values

# Perform PCA on the data
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

mpl = False 
if mpl:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # Create the ground truth plane with z = 0
    # using n N x N grid
    N = 6
    xxt, yyt = np.meshgrid(np.linspace(-2.25, 0, N), np.linspace(-2.25, 2.25, N))
    xxb, yyb = np.meshgrid(np.linspace(0, 2.25, N), np.linspace(-2.25, 2.25, N))
 
    # If we use mpl, we have issues with intersecting planes
    # due to the way that matplotlib draws in 3d 

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
        
    p1 = np.array([x_pca_plane[0,0], y_pca_plane[0,0], z_pca_plane[0,0]])
    p2 = np.array([x_pca_plane[0,1], y_pca_plane[0,1], z_pca_plane[0,1]])
    p3 = np.array([x_pca_plane[1,0], y_pca_plane[1,0], z_pca_plane[1,0]])
    
    v1 = p3 - p1
    v2 = p2 - p1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = np.dot(cp, p3)
    
    zzt = (d - a * xxt - b * yyt) / c
    zzb = (d - a * xxb - b * yyb) / c

    ax.scatter(xs = data[:,0], ys = data[:,1], zs = data[:,2], color = 'orange')
    ax.plot_surface(X = xxt, Y = yyt, Z = zzt, color = 'red')
    ax.plot_surface(X = xxt, Y = yyt, Z = np.zeros(xxt.shape), color = 'blue')
    ax.plot_surface(X = xxb, Y = yyb, Z = np.zeros(xxb.shape), color = 'blue')
    ax.plot_surface(X = xxb, Y = yyb, Z = zzb, color = 'red')
   
    #plt.savefig('uncertainty.svg')
    plt.show()
else:

    # Plotly works but saving as an svg is useless
    # so we resort to saving as a png
    import plotly.offline as pyo
    import plotly.graph_objs as go

    N = 6
    x_lin = np.linspace(-2.25, 2.25, N)
    y_lin = np.linspace(-2.25, 2.25, N)
    xx, yy = np.meshgrid(x_lin, y_lin)
    zz = np.zeros(xx.shape)
 
    camera = dict(
            up = dict(x = 0, y = 0, z = 1),
            center = dict(x = 0, y = 0, z = 0),
            eye = dict(x = -0.5, y = 2, z = 0.9)
    )

    points = go.Scatter3d(x = data[:,0], y = data[:, 1], z = data[:, 2], mode = 'markers', marker=dict(color = 'orange', size = 10))
    true_surface = go.Surface(x = xx, y = yy, z = zz, showscale = False, opacity = 0.75, colorscale = 'Greys', showlegend = True)
    pca_surface = go.Surface(x = x_pca_plane, y = y_pca_plane, z = z_pca_plane, showscale = False, opacity = 0.75, surfacecolor = 255 * np.ones(x_pca_plane.shape), colorscale = 'orange')
    
    xaxis = go.XAxis(title = '', showaxeslabels = False, showticklabels = False)
    yaxis = go.YAxis(title = '', showaxeslabels = False, showticklabels = False)
    zaxis = go.ZAxis(title = '', showaxeslabels = False, showticklabels = False)
    scene = go.Scene(camera = camera, xaxis = xaxis, yaxis = yaxis, zaxis = zaxis)
    layout = go.Layout(scene=scene, width = 2000, height = 2000)

    fig = go.Figure(data = [points, true_surface, pca_surface], layout = layout)
    pyo.plot(fig)
