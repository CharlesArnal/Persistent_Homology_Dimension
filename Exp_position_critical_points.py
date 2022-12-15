import numpy                 as np
import matplotlib.pyplot as plt

from points_clouds_generation import points_on_J_T_parabolas_example, points_on_generalized_critical_configuration, points_on_interlaced_2_spheres_in_R5,\
 points_on_interlaced_2_circles_in_R3, points_on_parallel_cubes_in_RD, points_on_regular_grid_on_parallel_cubes_in_RD, points_on_two_segments_in_R3, \
    points_on_torus_in_R3, points_on_lissajous_curve_in_R2


from explore_critical_simplices import plot_crit_points_3D_points_cloud, plot_crit_circumspheres_3D_points_cloud
from persistence_computations import  plot_diagrams
from basic_functions import compute_alpha_complex, dgm_square_to_correct_values

# Experiments to study the position of critical points (i.e. the centers of the circumspheres of critical simplices)


seed = np.random.randint(0,10000)
print(f"Random seed = {seed}")
np.random.seed(seed)

# Choose which experiment to run
experiment_number = 1


if experiment_number == 1:


    N = 1000
    d = 2
    D = 3
    precision = 1e-15
    length_intervals = 1e-10
    length = 3


    X = points_on_parallel_cubes_in_RD(N,d,D,length)

    positive_or_negative = "positive"
    dimension = 1

    Acomplex = compute_alpha_complex(X, precision)
    st = Acomplex.create_simplex_tree()
    dgm= st.persistence()
    dgm = dgm_square_to_correct_values(dgm)

    fig = plt.figure( num = f"Critical points for {positive_or_negative} simplices of dimension {dimension} on two parallel squares")
    ax = plt.axes(projection='3d')
    plot_crit_points_3D_points_cloud(st, Acomplex, dimension, ax, positive_or_negative)
    ax.scatter3D(X[:,0],X[:,1],X[:,2])
    plt.show()
    plot_diagrams(dgm)
   
if experiment_number == 2:

    N = 1000
    d = 2
    D = 3
    precision = 1e-15
    length_intervals = 1e-10
    length = 3
    epsilon = 0.1


    X = points_on_regular_grid_on_parallel_cubes_in_RD(epsilon,d,D,length = 3)

    positive_or_negative = "positive"
    dimension = 1

    Acomplex = compute_alpha_complex(X, precision)
    st = Acomplex.create_simplex_tree()
    dgm = st.persistence()
    dgm = dgm_square_to_correct_values(dgm)

    fig = plt.figure( num = f"Critical points for {positive_or_negative} simplices of dimension {dimension} on two parallel squares with regular sampling")
    ax = plt.axes(projection='3d')
    plot_crit_points_3D_points_cloud(st, Acomplex, dimension, ax, positive_or_negative)
    ax.scatter3D(X[:,0],X[:,1],X[:,2])
    plt.show()

    plot_diagrams(dgm)

if experiment_number == 3:

    N = 1000
    d = 2
    D = 3
    precision = 1e-15
    length_intervals = 1e-10
    length = 3
    epsilon = 0.1


    X = points_on_regular_grid_on_parallel_cubes_in_RD(epsilon,d,D,length = 3)

    positive_or_negative = "positive"
    dimension = 1

    Acomplex = compute_alpha_complex(X, precision)
    st = Acomplex.create_simplex_tree()
    dgm = st.persistence()
    dgm = dgm_square_to_correct_values(dgm)

    fig = plt.figure( num = f"Critical points for {positive_or_negative} simplices of dimension {dimension} on two parallel squares with regular sampling")
    ax = plt.axes(projection='3d')
    plot_crit_points_3D_points_cloud(st, Acomplex, dimension, ax, positive_or_negative)
    ax.scatter3D(X[:,0],X[:,1],X[:,2])
    plt.show()

    plot_diagrams(dgm)