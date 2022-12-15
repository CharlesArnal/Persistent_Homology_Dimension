import numpy                 as np

from points_clouds_generation import points_on_interlaced_2_circles_in_R3, points_on_interlaced_2_spheres_in_R5, points_on_parallel_cubes_in_RD,\
    points_on_regular_grid_on_parallel_cubes_in_RD
from persistence_computations import compute_total_persistence, compute_alpha_complex, dgm_square_to_correct_values
from explore_critical_simplices import plot_crit_points_3D_points_cloud
import matplotlib.pyplot as plt

np.random.seed(0)


# Experiment 1 : plot the critical points

N = 800
d = 2
D = 3
epsilon = 0.5
X = points_on_parallel_cubes_in_RD(N, d, D)
Acomplex = compute_alpha_complex(X)
st = Acomplex.create_simplex_tree()
dgm = st.persistence()
dgm = dgm_square_to_correct_values(dgm)

positive_or_negative = "positive"
dimension = 1

fig = plt.figure( num = f"Critical points for {positive_or_negative} simplices of dimension {dimension}")
ax = plt.axes(projection='3d')
plot_crit_points_3D_points_cloud(st, Acomplex, dimension, ax, positive_or_negative)
ax.scatter3D(X[:,0],X[:,1],X[:,2])
plt.show()

# Experiment 2 : compute the persistent homology
 

epsilon = 0.1
homology_degree = 2
order = 2
precision = 1e-5

for N in [10, 100, 300, 1000]:
    X = points_on_interlaced_2_circles_in_R3(N, epsilon)    
    Acomplex = compute_alpha_complex(X, precision)
    st = Acomplex.create_simplex_tree()
    dgm = st.persistence()
    print(f"--- persistance = {compute_total_persistence(dgm, homology_degree, order)}")




