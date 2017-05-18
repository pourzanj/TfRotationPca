import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set random seed for reproducibility
np.random.seed(10)

# Sample uniformly on plane for theta, phi
def sample_angle_points(n, r):
    theta = np.random.uniform(0, 2*np.pi, size = n)
    phi = np.random.uniform(0, np.pi, size = n)
    
    return theta, phi
    
# Translate plane to spherical coordinates
def angles_to_spherical_coordinates(theta, phi, r):
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    
    return x, y, z

def generate_plot_plane_sphere(n, r = 1, width = 8, height = 6,
                               point_color='black', grid_color='orange'):
    theta, phi = sample_angle_points(n = n, r = r)
    x, y, z = angles_to_spherical_coordinates(theta, phi, r = r)
    
    # Plot sphere with points
    fig = plt.figure(figsize = (width, height))
    ax = fig.gca(projection = '3d')
    u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
    xs, ys, zs = angles_to_spherical_coordinates(u, v, r = r)

    ax.scatter(x, y, z, color = point_color, s = 3)
    ax.plot_wireframe(xs, ys, zs, color = grid_color)
    plt.axis('off')
    plt.savefig('sphere.svg')
    plt.close()
    
    # Plot plane
    fig = plt.figure(figsize = (width, height))
    ax = fig.gca()

    plt.grid(True)
    plt.rc('grid', linestyle='-', color= grid_color)
    plt.scatter(theta, phi, color = marker_color)
    plt.savefig('plane.svg', transparent = True)

generate_plot_plane_sphere(1000)
