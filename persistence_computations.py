import numpy                 as np
import gudhi                 as gd
import gudhi.wasserstein
import matplotlib.pyplot     as plt

from points_clouds_generation import points_on_interlaced_2_spheres_in_R5
from basic_functions import compute_alpha_complex, dgm_square_to_correct_values
from explore_critical_simplices import get_critical_pairs, get_lengths_intervals






def plot_diagrams(dgm):
    max_degree = max([point[0] for point in dgm])
    #fig, axes = plt.subplots(1, max_dimension)
    figs = []
    for degree in range(0,max_degree+1):
        dgm_in_chosen_degree = []
        for point in dgm:
            if point[0] == degree:
                dgm_in_chosen_degree.append(point[1])
        _, ax = plt.subplots(num = f"Degree {degree}" )
        gd.plot_persistence_diagram(dgm_in_chosen_degree, axes = ax)
    plt.show()


def get_number_of_intervals(simplex_tree_with_computed_persistence, homology_degree : int, min_len_intervals = 1e-10, min_birth_value = 0.3, essential_part_included = False):
    relevant_pairs = get_critical_pairs(simplex_tree_with_computed_persistence, homology_degree, min_len_intervals, min_birth_value, essential_part_included)
    return len(relevant_pairs)


def get_mean_and_std_of_intervals_lengths(simplex_tree_with_computed_persistence, homology_degree : int, min_len_intervals = 1e-10, min_birth_value = 0.3, essential_part_included = False):
    relevant_pairs = get_critical_pairs(simplex_tree_with_computed_persistence, homology_degree, min_len_intervals, min_birth_value, essential_part_included)
    lengths = np.array(get_lengths_intervals(relevant_pairs, simplex_tree_with_computed_persistence))
    return np.mean(lengths), np.std(lengths)



def compute_total_persistence(dgm, homology_degree : int, order : int) -> float:
    """
    dgm must be the product of a gudhi simplex_tree persistence() method

    Returns the total persistence of the homology of degree "homology_degree" of the alpha filtration of X
    w.r.t the Wasserstein distance of order "order"

    If plot = True, the function also plots the persistence diagram of given degree
    """
    dgm_in_correct_degree = []

    for point in dgm :
        if point[0] == homology_degree:
            dgm_in_correct_degree.append(point[1])
    dgm_in_correct_degree = np.array(dgm_in_correct_degree)

    print("Computing Wasserstein distance")
    return gudhi.wasserstein.wasserstein_distance(dgm_in_correct_degree,np.array([]), order = order, keep_essential_parts=True)



if __name__ == "__main__":
    N = 80
    homology_degree = 1
    order = 2
    epsilon = 1
    precision = 1e-15
    X = points_on_interlaced_2_spheres_in_R5(N, epsilon)

    Acomplex = compute_alpha_complex(X,precision)
    st = Acomplex.create_simplex_tree()
    print("Computing persistence diagram")
    dgm = st.persistence()
    dgm = dgm_square_to_correct_values(dgm)
    plot_diagrams(dgm)
    print(f"Wasserstein distance = {compute_total_persistence(dgm, homology_degree, order)}")

