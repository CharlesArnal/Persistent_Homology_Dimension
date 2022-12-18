import numpy                 as np
import matplotlib.pyplot as plt

from points_clouds_generation import points_on_J_T_parabolas_example, points_on_generalized_critical_configuration, points_on_interlaced_2_spheres_in_R5,\
 points_on_interlaced_2_circles_in_R3, points_on_parallel_cubes_in_RD, points_on_regular_grid_on_parallel_cubes_in_RD, points_on_two_segments_in_R3, \
    points_on_torus_in_R3, points_on_lissajous_curve_in_R2


from explore_critical_simplices import plot_crit_points_3D_points_cloud, plot_crit_circumspheres_3D_points_cloud, h_index_of_points
from persistence_computations import  plot_diagrams
from basic_functions import compute_alpha_complex, dgm_square_to_correct_values

# Experiments to study the position of critical points (i.e. the centers of the circumspheres of critical simplices)


seed = np.random.randint(0,10000)
print(f"Random seed = {seed}")
np.random.seed(seed)

# Choose which experiment to run
experiment_number = 4


if experiment_number == 1:


    N = 2000
    d = 2
    D = 3
    precision = 1e-15
    length_intervals = 1e-10
    length = 2


    X = points_on_parallel_cubes_in_RD(N,d,D,length)

    positive_or_negative = "positive"
    dimension = 2

    Acomplex = compute_alpha_complex(X, precision)
    st = Acomplex.create_simplex_tree()
    dgm= st.persistence()
    dgm = dgm_square_to_correct_values(dgm)

    fig = plt.figure( num = f"Critical points for {positive_or_negative} simplices of dimension {dimension} on two parallel squares")
    ax = plt.axes(projection='3d')
    plot_crit_points_3D_points_cloud(st, Acomplex, dimension, ax, positive_or_negative, min_len_intervals=length_intervals)
    ax.scatter3D(X[:,0],X[:,1],X[:,2])
    plt.show()
    print(f"h-index des points {h_index_of_points(st, Acomplex, dimension, positive_or_negative, min_len_intervals=length_intervals)}")

    plot_diagrams(dgm)
   
if experiment_number == 2:

    N = 2000
    d = 2
    D = 3
    precision = 1e-15
    length_intervals = 1e-10
    length = 2
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
    plot_crit_points_3D_points_cloud(st, Acomplex, dimension, ax, positive_or_negative, min_len_intervals=length_intervals)
    ax.scatter3D(X[:,0],X[:,1],X[:,2])
    plt.show()

    print(f"h-index des points {h_index_of_points(st, Acomplex, dimension, positive_or_negative, min_len_intervals=length_intervals)}")

    plot_diagrams(dgm)

if experiment_number == 3:

    N = 1000
    precision = 1e-15
    length_intervals = 1e-10
    min_birth_value = 0.2

    X = points_on_interlaced_2_circles_in_R3(N, 0.5)

    positive_or_negative = "positive"
    dimension = 1

    Acomplex = compute_alpha_complex(X, precision)
    print("Creating simplex tree")
    st = Acomplex.create_simplex_tree()
    print("Computing diagram")
    dgm = st.persistence()
    dgm = dgm_square_to_correct_values(dgm)

    fig = plt.figure( num = f"Critical points for {positive_or_negative} simplices of dimension {dimension} on the interlaced circles")
    ax = plt.axes(projection='3d')
    plot_crit_points_3D_points_cloud(st, Acomplex, dimension, ax, positive_or_negative, min_len_intervals=length_intervals, min_birth_value=min_birth_value)
    ax.scatter3D(X[:,0],X[:,1],X[:,2])
    plt.show()

    print(f"h-index des points {h_index_of_points(st, Acomplex, dimension, positive_or_negative, min_len_intervals=length_intervals)}")

    #plot_diagrams(dgm)

if experiment_number == 4:

    N = 1000
    # Careful, epsilon is not the same at all as above
    epsilon = 0.001
    precision = 1e-15
    length_intervals = 1e-10
    min_birth_value = 0.2

    X = points_on_J_T_parabolas_example(N, 0.005)

    positive_or_negative = "positive"
    dimension = 1

    Acomplex = compute_alpha_complex(X, precision)
    print("Creating simplex tree")
    st = Acomplex.create_simplex_tree()
    print("Computing diagram")
    dgm = st.persistence()
    dgm = dgm_square_to_correct_values(dgm)

    fig = plt.figure( num = f"Critical points for {positive_or_negative} simplices of dimension {dimension} on the interlaced parabolas")
    ax = plt.axes(projection='3d')
    plot_crit_points_3D_points_cloud(st, Acomplex, dimension, ax, positive_or_negative, min_len_intervals=length_intervals, min_birth_value=min_birth_value)
    ax.scatter3D(X[:,0],X[:,1],X[:,2])
    plt.show()

    print(f"h-index des points {h_index_of_points(st, Acomplex, dimension, positive_or_negative, min_len_intervals=length_intervals)}")

    plot_diagrams(dgm)