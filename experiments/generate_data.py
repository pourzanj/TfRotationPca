import numpy as np
import plotly as py
import plotly.graph_objs as go

def generate_ppca(W, L, n_points, noise = 1):

    n, p = W.shape

    z = np.random.multivariate_normal(np.zeros(p), np.eye(p), size = n_points)
    x =  np.dot(np.dot(W, L), z.T).T + np.random.normal(0, noise, size = (n_points, n))

    return x

def pca():
    W = np.array([[0.5, np.sqrt(2)/2],[-0.5, np.sqrt(2)/2],[np.sqrt(2)/2, 0]])
    L = np.diag([5,3])
    x = generate_ppca(W, L, 100)

    return x

def cross(numA, numB):
    WA = np.array([[0.5, np.sqrt(2)/2],[-0.5, np.sqrt(2)/2],[np.sqrt(2)/2, 0]])
    WB = np.array([[0.5, -np.sqrt(2)/2],[0.5, np.sqrt(2)/2],[np.sqrt(2)/2, 0]])
    L = np.diag([5, 3])

    xA = generate_ppca(WA, L, numA, noise = 0.5)
    xB = generate_ppca(WB, L, numB, noise = 0.5)

    return np.vstack([xA, xB])


pca_data = pca()
cross_data = cross(75, 50)

np.savetxt('pca.txt', pca_data)
np.savetxt('cross.txt', cross_data)
