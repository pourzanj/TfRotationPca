import numpy as np
import autograd as ag
import autograd.numpy as agnp
import pandas as pd
import pickle

from scipy.linalg import expm

from multiprocessing import Process, Manager

def Embedded_Manifold_HMC_MT(initial_params, n_iters, epsilon, L, U, grad_U, flow, projection, proc_num, return_dict, print_iter = 10000):
    np.random.seed()
    chain = [initial_params]
    accepts = 0
    rejects = 0

    for it in range(1, n_iters):
        
        q = chain[it - 1].copy()
        p = np.random.normal(0, 1, size = q.shape)
        p = projection(p, q)

        # Compute current potential and kinetic energy
        current_U  = U(q)
        current_K = np.sum(p**2) / 2
       
        for i in range(L):
            p = p - epsilon * grad_U(q) / 2
            p = projection(p, q)
            q, p = flow(q, p, epsilon)
            p = p - epsilon * grad_U(q) / 2
            p = projection(p,q)

        proposed_U = U(q)
        proposed_K = np.sum(p**2) / 2

        if np.random.uniform(0,1) < np.exp(current_U - proposed_U + current_K - proposed_K):
            chain.append(q)
            accepts += 1
        else:
            chain.append(chain[it-1])
            rejects += 1

        if it % print_iter == 0:
            print(it)

    return_dict[proc_num] = chain

# The projection onto the stiefel manifold
# Input: 
#   V: Velocity
#   X: Position
def projection_Stiefel(V,X):
    inner = np.dot(X.T,V) + np.dot(V.T, X)
    return V - (1/2) * np.dot(X,inner)

# The geodesic flow of the particle on the Stiefel Manifold
# Input:
#   X: Position
#   V: Velocity
#   epsilon: The time-step.
def flow_Stiefel(X,V,epsilon):
    A = np.dot(X.T,V)
    S = np.dot(V.T,V)
    mexp = expm(-epsilon * A)
    t1 = expm(epsilon * np.bmat([[A, -S], [np.eye(A.shape[0]), A]]))
    t2 = np.bmat([[mexp, np.zeros(mexp.shape)],[np.zeros(mexp.shape), mexp]])
    F0 = np.bmat([X, V])
    R = np.dot(F0, np.dot(t1,t2))
    X_up, V_up = np.hsplit(R, 2)
    return np.array(X_up), np.array(V_up)

"""
Helper functions to transform angles to the Given's representation
and vice versa
"""
def left_rotate_counter_clockwise(A, angle, i, j):
    AR = A.copy()
    AR[i,:] = np.cos(angle)*A[i,:] - np.sin(angle)*A[j,:]
    AR[j,:] = np.sin(angle)*A[i,:] + np.cos(angle)*A[j,:]
    
    return AR

def right_rotate_counter_clockwise(A, angle, i, j):
    AR = A.copy()
    AR[:,i] = np.cos(angle)*A[:,i] + np.sin(angle)*A[:,j]
    AR[:,j] = -np.sin(angle)*A[:,i] + np.cos(angle)*A[:,j]
    
    return AR

def inverse_givens_transform(angles, n, p):
    G = np.eye(n)
    idx = 0
    for i in range(p):
        for j in range(i+1, n):
            G = right_rotate_counter_clockwise(G, angles[idx], i, j)
            idx = idx + 1
    return G[:,0:p]

def givens_transform(W):
    n, p = W.shape
    angles = [0 for _ in range(int(n*p-p*(p+1)/2))]
    idx = 0
    for i in range(p):
        for j in range(i+1, n):
            angle = np.arctan2(W[j,i],W[i,i])
            W = left_rotate_counter_clockwise(W, -angle, i, j)
            angles[idx] = angle
            idx = idx + 1
    
    return [a/np.pi for a in angles]

data = pd.read_csv('x.csv', index_col = 0)
data = data.values
N = data.shape[0]

# The model is y ~ N(Wz + \epsilon)
# This is the PPCA model and the log-likelihood is specified by Tipping and Bishop
# We build the potential energy = -log-likelihood and use autograd to find the derivative
sigma_hat = (1/N) * np.dot(data.T,data)
z = np.diag(np.array([1.,1.]))
def U(Q):
    C = agnp.dot(Q, agnp.dot(z,Q.T)) + np.eye(Q.shape[0])
    C_inv = agnp.linalg.inv(C)
    return N/2 * (agnp.log(agnp.linalg.det(C)) + agnp.trace(agnp.dot(C_inv,sigma_hat)))
grad_U = ag.grad(U)

# Setting up the params for the sampling
n_dim = 3
p_dim = 2
initial_params = np.eye(n_dim)[:,:p_dim]
n_iters = 10
epsilon = 5e-2
L = 50
flow = flow_Stiefel
projection = projection_Stiefel
n_chains = 3

jobs = []
m = Manager()
return_dict = m.dict()
for i in range(n_chains):
    p = Process(target = Embedded_Manifold_HMC_MT, args = (initial_params, n_iters, epsilon, L, U, grad_U, flow, projection, i, return_dict))
    jobs.append(p)
    p.start()
    
for i in range(n_chains):
    p.join()

samples0 = []
for i in range(n_chains):
    samples0 += return_dict[i]

with open('fit.pkl', 'wb') as f:
    pickle.dump(samples0, f)
