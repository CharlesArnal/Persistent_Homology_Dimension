
import numpy                 as np
import matplotlib.pyplot     as plt
import gudhi                 as gd
import random 

from points_clouds_generation import points_on_parallel_cubes_in_RD, points_on_interlaced_2_circles_in_R3
from basic_functions import test_criticality, circumcenter_of_simplex, compute_alpha_complex, get_critical_pairs, get_lengths_intervals



def h_index_of_points(simplex_tree_with_persistence_already_computed, Acomplex, dimension : int, positive_or_negative = "positive", min_len_intervals = 1e-10, min_birth_value = 0.0):
    """
    Takes as input a simplex_tree st with persistence already computed, and its associated alpha complex Acomplex

    The function returns the greatest integer h such that there are at least h points of the poins cloud that each belong to at least h critical simplices of the given dimension

    If positive_or_negative = "positive", these simplices will correspond to birth point of intervals in the homology of degree dimension

    If positive_or_negative = "positive", they will correspond to death point of intervals in the homology of degree (dimension - 1)

    A simplex is considered critical if it gives birth or death to an interval of length at least len_intervals

    Only simplices born after min_birth_value are considered
    """
    
    st = simplex_tree_with_persistence_already_computed

    if positive_or_negative == "positive":
        homology_degree = dimension
    else:
        homology_degree = dimension - 1

    # persistence_pairs are pairs of negative and negative simplices giving birth and death to an interval
    # simplices are represented as sets of integers (each corresponding to a point in the points cloud)
    critical_pairs_in_correct_degree = get_critical_pairs(st, homology_degree, min_len_intervals, min_birth_value, False)
    critical_simpl_in_correct_degree = []
    if positive_or_negative == "positive":
        critical_simpl_in_correct_degree = [pair[0] for pair in critical_pairs_in_correct_degree]
    else:
        critical_simpl_in_correct_degree = [pair[1] for pair in critical_pairs_in_correct_degree]

    # Go through all the simplices, count the number of appearances of each point 
    # The points are stored (as integers) in a list
    # Their number of appearances as integers in another list (such that the i-th entries of both lists correspond to each other)
    # Not the most efficient way of doing it
    if len(critical_simpl_in_correct_degree) == 0:
        return 0
    points = []
    appearances = []
    for simplex in critical_simpl_in_correct_degree:
        for point in simplex:
            if point in points:
                index = points.index(point)
                appearances[index] += 1
            else:
                points.append(point)
                appearances.append(1)

    appearances.sort(reverse=True)
    h_index = 0
    while(appearances[h_index]>=h_index+1):
        h_index+=1
    return h_index



def plot_crit_points_3D_points_cloud(simplex_tree_with_persistence_already_computed, Acomplex, dimension : int, axes_3D, positive_or_negative = "positive", min_len_intervals = 1e-10, min_birth_value = 0.0):
    """
    Takes as input a simplex_tree st with persistence already computed, and its associated alpha complex Acomplex

    axes_3D is a matplotlib axes obtained with fig = plt.figure(),  axes_3D = plt.axes(projection='3d')

    The function plots on axes_3D the centers of the critical positive (resp, negative) simplices of specified dimension in the persistence of the simplex_tree

    The function does NOT call matplotlib.pyplot.show()

    If positive_or_negative = "positive", these simplices will correspond to birth point of intervals in the homology of degree dimension

    If positive_or_negative = "positive", they will correspond to death point of intervals in the homology of degree (dimension - 1)

    A simplex is considered critical if it gives birth or death to an interval of length at least len_intervals

    Only simplices born after min_birth_value are considered
    
    """
    
    print(f"Plotting the critical points for {positive_or_negative} simplices of dimension {dimension} :")
    st = simplex_tree_with_persistence_already_computed

    if positive_or_negative == "positive":
        homology_degree = dimension
    else:
        homology_degree = dimension - 1

    # persistence_pairs are pairs of negative and negative simplices giving birth and death to an interval
    # simplices are represented as sets of integers (each corresponding to a point in the points cloud)
    critical_pairs_in_correct_degree = get_critical_pairs(st, homology_degree, min_len_intervals, min_birth_value, False)
    critical_simpl_in_correct_degree = []
    if positive_or_negative == "positive":
        critical_simpl_in_correct_degree = [pair[0] for pair in critical_pairs_in_correct_degree]
    else:
        critical_simpl_in_correct_degree = [pair[1] for pair in critical_pairs_in_correct_degree]
    # the criticality of a simplex is simply the length of the interval to which it gives birth (respectively death)
    criticality = get_lengths_intervals(critical_pairs_in_correct_degree, st)

    critical_points = np.array([circumcenter_of_simplex(np.array([Acomplex.get_point(point) for point in simplex])) for simplex in critical_simpl_in_correct_degree])
    criticality = np.array(criticality)
    axes_3D.scatter3D(critical_points[:,0],critical_points[:,1],critical_points[:,2],c = criticality)



def plot_crit_circumspheres_3D_points_cloud(N_spheres : int, simplex_tree_with_persistence_already_computed, Acomplex, dimension : int, axes_3D,\
    positive_or_negative = "positive", min_radius_simplex = 0.8, len_intervals = 1e-10):
    """
    Takes as input a simplex_tree st with persistence already computed, and its associated alpha complex Acomplex

    axes_3D is a matplotlib axes obtained with fig = plt.figure(),  axes_3D = plt.axes(projection='3d')

    The function plots on axes_3D N_spheres circumspheres of the critical positive (resp, negative) simplices of specified dimension in the persistence of the simplex_tree

    The simplices whose circumspheres are plotted are chosen at random; only simplices whose circumsphere has radius greater than min_radius_simplex can be selected

    The function does NOT call matplotlib.pyplot.show()

    If positive_or_negative = "positive", these simplices will correspond to birth point of intervals in the homology of degree dimension

    If positive_or_negative = "positive", they will correspond to death point of intervals in the homology of degree (dimension - 1)

    A simplex is considered critical if it gives birth or death to an interval of length at least len_intervals
    
    """
    print(f"Plotting the circumspheres of some {positive_or_negative} critical simplices of dimension {dimension} :")

    st = simplex_tree_with_persistence_already_computed

    if positive_or_negative == "positive":
        homology_degree = dimension
    else:
        homology_degree = dimension - 1

    # persistence_pairs are pairs of negative and negative simplices giving birth and death to an interval
    # simplices are represented as sets of integers (each corresponding to a point in the points cloud)
    critical_pairs_in_correct_degree = get_critical_pairs(st, homology_degree, len_intervals, min_radius_simplex, False)
    critical_simpl_in_correct_degree = []
    if positive_or_negative == "positive":
        critical_simpl_in_correct_degree = [pair[0] for pair in critical_pairs_in_correct_degree]
    else:
        critical_simpl_in_correct_degree = [pair[1] for pair in critical_pairs_in_correct_degree]
    
   
    # Randomly select (at most) N_spheres simplices whose circumsphere will be plotted :
    chosen_simplices = random.choices(critical_simpl_in_correct_degree, k = min(N_spheres, len(critical_simpl_in_correct_degree)))
    """
    chosen_simplices = []
    random.shuffle(critical_simpl_in_correct_degree)
    i = 0
    while i<len(critical_simpl_in_correct_degree) and len(chosen_simplices) < N_spheres :
        simplex = critical_simpl_in_correct_degree[i]
        points = np.array([Acomplex.get_point(point) for point in simplex])
        center_circumsphere = circumcenter_of_simplex(points)
        radius = np.linalg.norm(center_circumsphere-points[0])
        if radius >= min_radius_simplex:
            chosen_simplices.append(simplex)
        i += 1 
    """

    for simplex in chosen_simplices:
        points = np.array([Acomplex.get_point(point) for point in simplex])
        center_circumsphere = circumcenter_of_simplex(points)
        radius = np.linalg.norm(center_circumsphere-points[0])
        u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
        x = np.cos(u)*np.sin(v)*radius + center_circumsphere[0]
        y = np.sin(u)*np.sin(v)*radius + center_circumsphere[1]
        z = np.cos(v)*radius + center_circumsphere[2]
        # alpha controls opacity
        axes_3D.plot_surface(x, y, z, color="g", alpha=0.3)
        axes_3D.scatter3D(points[:,0],points[:,1],points[:,2], s=50 , alpha=1, color = "red")




if __name__ == "__main__":
    #np.random.seed(0)
    N = 300
    d = 2
    D = 3
    epsilon = 0.5
    precision = 1e-15

    #X = points_on_parallel_cubes_in_RD(N, d, D)
    X = points_on_interlaced_2_circles_in_R3(N, epsilon)


    positive_or_negative = "positive"
    dimension = 1

    Acomplex = compute_alpha_complex(X, precision)
    print("Computing persistence diagram")
    st = Acomplex.create_simplex_tree()
    st.compute_persistence()

    fig = plt.figure( num = f"Critical points for {positive_or_negative} simplices of dimension {dimension}")
    ax = plt.axes(projection='3d')
    plot_crit_points_3D_points_cloud(st, Acomplex, dimension, ax, positive_or_negative)
    N_spheres = 2
    min_radius_simplex = 0.2
    plot_crit_circumspheres_3D_points_cloud(N_spheres, st, Acomplex, dimension, ax, positive_or_negative, min_radius_simplex)
    ax.scatter3D(X[:,0],X[:,1],X[:,2])
    plt.show()


    
