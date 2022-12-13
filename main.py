import numpy                 as np

from points_clouds_generation import points_on_interlaced_2_circles_in_R3, points_on_interlaced_2_spheres_in_R5, points_on_parallel_hyperplanes_in_Rd
from compute_total_persistence import compute_total_persistence_alpha
from plot_crit_points import plot_crit_points_3D_points_cloud

np.random.seed(0)


# Experience 1 : plot the critical points

N = 800
d = 3
epsilon = 0.5
X = points_on_parallel_hyperplanes_in_Rd(N, d)

positive_or_negative = "positive"
dimension = 1

plot_crit_points_3D_points_cloud(X, dimension, positive_or_negative)


# Experience 2 : compute the persistent homology
 

epsilon = 0.1
homology_degree = 2
order = 2
precision = 1e-5

for N in [10, 100, 300, 1000]:
    X = points_on_interlaced_2_circles_in_R3(N, epsilon)    
    print(f"--- persistance = {compute_total_persistence_alpha(X, homology_degree, order, precision, plot = False)}")




