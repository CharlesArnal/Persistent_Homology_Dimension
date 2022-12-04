import numpy                 as np
import matplotlib.pyplot     as plt
import pandas                as pd
import gudhi                 as gd
from tqdm                    import tqdm

import gudhi.wasserstein


import matplotlib


np.random.seed(0)

# Samples points on two interlaced 2-spheres in R^5, and computes the total persistance of their alpha filtration w.r.t the Wasserstein distance of a certain order


def compute_total_persistence(N,epsilon,homology_degree,order):
    print(f"\nN = {N}")
    print(f"Sampling points")
    # centre première sphère = [1,0,0,0,0], plan x1,x2,x3
    X_1 =  np.array(np.concatenate((np.random.randn(N,3), np.zeros((N,2))), axis=1))

    X_1 = (X_1.T/np.linalg.norm(X_1,axis=1)).T
    centers_1 = np.concatenate((np.ones((N,1)),np.zeros((N,4))), axis = 1)
    X_1 = X_1 + centers_1

    # Se restreindre à des poins proches du centre de la seconde sphère
    X_1_refined = []
    for x in X_1:
        if np.linalg.norm(x)<epsilon:
            X_1_refined.append(x)
    X_1_refined = np.array(X_1_refined)


    # centre seconde sphère = [0,0,0,0,0], plan x1, x4, x5

    X_2 =  np.array(np.concatenate((np.random.randn(N,1), np.zeros((N,2)),np.random.randn(N,2)), axis=1))

    X_2 = (X_2.T/np.linalg.norm(X_2,axis=1)).T

    # Se restreindre à des poins proches du centre de la seconde première
    X_2_refined = []
    for x in X_2:
        if np.linalg.norm(x-[1,0,0,0,0])<epsilon:
            X_2_refined.append(x)
    X_2_refined = np.array(X_2_refined)


    X_total = np.concatenate((X_1_refined,X_2_refined), axis = 0)


    print(f"Number of sampled points = {np.shape(X_total)[0]}")

    print("Computing persistence diagram")

    st = gd.AlphaComplex(points=X_total).create_simplex_tree()
    dgm = st.persistence()

    dgm_in_correct_degree = []

    for point in dgm :
        if point[0] == homology_degree:
            dgm_in_correct_degree.append(point[1])
    dgm_in_correct_degree = np.array(dgm_in_correct_degree)
    print("Computing Wasserstein distance")
    #return gd.wasserstein.wasserstein_distance(dgm,dgm, order = order, keep_essential_parts=True)
    return gd.wasserstein.wasserstein_distance(dgm_in_correct_degree,np.array([]), order = order, keep_essential_parts=True)




epsilon = 0.1
homology_degree = 1
order = 2



for N in [100, 1000, 10000, 100000, 200000]:
    
    print(f"--- persistance = {compute_total_persistence(N,epsilon, homology_degree, order)}")

