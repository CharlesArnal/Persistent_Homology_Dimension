import numpy                 as np

from points_clouds_generation import points_on_J_T_parabolas_example, points_on_generalized_critical_configuration, points_on_interlaced_2_spheres_in_R5,\
 points_on_interlaced_2_circles_in_R3, points_on_parallel_cubes_in_RD, points_on_regular_grid_on_parallel_cubes_in_RD, points_on_two_segments_in_R3, \
    points_on_torus_in_R3, points_on_lissajous_curve_in_R2


from explore_critical_simplices import plot_crit_points_3D_points_cloud, plot_crit_circumspheres_3D_points_cloud
from persistence_computations import  plot_diagrams, get_number_of_intervals, get_mean_and_std_of_intervals_lengths
from basic_functions import compute_alpha_complex, dgm_square_to_correct_values

datasets_functions_with_fixed_dimensions = [ (2,points_on_interlaced_2_spheres_in_R5),(1, points_on_interlaced_2_circles_in_R3), (2,points_on_torus_in_R3)]


np.random.seed(0)
experiment_number = 1

# First batch of experiments : for each dataset, what is the maximum degree in which there is persistent homology ?

if experiment_number == 1:

    precision = 1e-15
    

    positive_or_negative = "positive"


    def get_max_degree(X):
        Acomplex = compute_alpha_complex(X,precision)
        st = Acomplex.create_simplex_tree()
        print("Computing persistence diagram")
        dgm = st.persistence()
        dgm = dgm_square_to_correct_values(dgm)
        max_degree = max([point[0] for point in dgm])
        plot_diagrams(dgm)
        return max_degree

    number_of_points = [10,30,100,300,1000]

    
    for N in number_of_points:
        X = points_on_interlaced_2_circles_in_R3(N, 0.5)
        print(f"For N = {N} and the two circles, the max degree is {get_max_degree(X)}")
        print("\n")

    """
    # Result : 2

    for N in number_of_points[:-1]:
        X = points_on_interlaced_2_spheres_in_R5(N, 0.5)
        print(f"For N = {N} and the two interlaced spheres, the max degree is {get_max_degree(X)}")
        print("\n")

    # Result : 4
    
    for N in number_of_points:
        X = points_on_J_T_parabolas_example(N, 0.05)
        print(f"For N = {N} and the two parabolas, the max degree is {get_max_degree(X)}")
        print("\n")

    # Result : 2
    

    for N in number_of_points:
        X = points_on_torus_in_R3(N, 1, 5)
        print(f"For N = {N} and the torus, the max degree is {get_max_degree(X)}")
        print("\n")
    
    # Result : 2
    
    for N in number_of_points+[2000]:
        X = points_on_parallel_cubes_in_RD(N,2,3,length = 2)
        print(f"For N = {N} and the two 2-squares in R^3, the max degree is {get_max_degree(X)}")
        print("\n")

    # bien vérifier si l'homologie persistente en degré 2 avec exposant 2 reste bornée (elle paraît très grande)

    # Result : 2

    for N in number_of_points+[2000]:
        X = points_on_parallel_cubes_in_RD(N,3,6,length = 2)
        print(f"For N = {N} and the two 3-cubes in R^6, the max degree is {get_max_degree(X)}")
        print("\n")

    # Result : 3
    
    for N in number_of_points+[2000]:
        epsilon = np.sqrt(8.0/float(N))
        X = points_on_regular_grid_on_parallel_cubes_in_RD(epsilon,2,3,length = 2)
        print(f"For N = {N} and the two 2-cubes in R^3 with regular grid, the max degree is {get_max_degree(X)}")
        print("\n")

    # Result : 2 (mais dur à dire avec la multiplicité)
   
   
    for N in [10,30,50,100][2:]:
        X = points_on_generalized_critical_configuration(N,2)
        print(f"For N = {N} and the generalised configuration with a 2 simplex (3 circles), the max degree is {get_max_degree(X)}")
        print("\n")

    # Result : 4 (behaviour more interesting in degree 3 and 4 than 0, 1 and 2)
"""
    for N in [10,30,50]:
        X = points_on_generalized_critical_configuration(N,3)
        print(f"For N = {N} and the generalised configuration with a 3 simplex (4 circles), the max degree is {get_max_degree(X)}")
        print("\n")

    # Result : 6 (behaviour more interesting in degree 5 and 6)
  
"""

points_on_J_T_parabolas_example, points_on_generalized_critical_configuration, points_on_interlaced_2_spheres_in_R5,\
 points_on_interlaced_2_circles_in_R3, points_on_parallel_cubes_in_RD, points_on_regular_grid_on_parallel_cubes_in_RD, points_on_two_segments_in_R3, \
    points_on_torus_in_R3, points_on_lissajous_curve_in_R2

# Second batch of experiments : estimate the number of non-zero intervals and their mean lengths for a few datasets and dimensions

if experiment_number == 2:
        
    epsilon = 0.5
    precision = 1e-15
    

    positive_or_negative = "positive"
    dimension = 1
    number_of_points = [10,30,100,300,1000,2000]

    number_of_intervals = dict()
    mean_lengths_intevals = dict()
    std_lengths_intervals = dict()


    for N in number_of_points:
        X = points_on_interlaced_2_circles_in_R3(N, epsilon)
        Acomplex = compute_alpha_complex(X,precision)
        st = Acomplex.create_simplex_tree()
        print("Computing persistence diagram")
        dgm = st.persistence()
        dgm = dgm_square_to_correct_values(dgm)
        max_degree = max([point[0] for point in dgm])
        for degree in max_degree : 
            number_of_intervals[degree] 







if exper
d = 2
D = 3
epsilon = 0.5
precision = 1e-15
X = points_on_parallel_cubes_in_RD(N, d, D)

positive_or_negative = "positive"
dimension = 1

print()

wasserstein_coeffs = [0.5,1,1.25,1.5,1.75,2]

for N in [10,30,100,300,1000]:

    Acomplex = compute_alpha_complex(X,precision)
    st = Acomplex.create_simplex_tree()
    print("Computing persistence diagram")
    dgm = st.persistence()
    dgm = dgm_square_to_correct_values(dgm)
    max_degree = max([point[0] for point in dgm])
    for degree in max_degree : 




"""