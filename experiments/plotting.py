import plotly.offline as pyo
import plotly.graph_objs as go
import numpy as np
from sklearn.decomposition import PCA

def surface(x, y, z):
    return go.Surface(x = x, y = y, z = z, opacity = 0.5, colorscale = "Viridis")

def scatter(points):
    return go.Scatter3d(x = points[:,0], y = points[:,1], z = points[:,2], mode = 'markers')

pca = np.loadtxt('pca/pca.txt')
model = PCA(n_components = 3)
model.fit(pca)
c = model.components_

x_pca_axis, y_pca_axis, z_pca_axis = 10 * c.T

x_pca_plane = np.r_[x_pca_axis[:2], - x_pca_axis[1::-1]]
y_pca_plane = np.r_[y_pca_axis[:2], - y_pca_axis[1::-1]]
z_pca_plane = np.r_[z_pca_axis[:2], - z_pca_axis[1::-1]]
x_pca_plane.shape = (2,2)
y_pca_plane.shape = (2,2)
z_pca_plane.shape = (2,2)

z2_pca_plane = z_pca_plane + 10 

data = [scatter(pca), surface(x_pca_plane, y_pca_plane, z_pca_plane), surface(x_pca_plane, y_pca_plane, z2_pca_plane)]
#layout = go.Layout()
fig = go.Figure(data = data)
pyo.plot(fig, filename = 'scatter_surface_plot.html')
