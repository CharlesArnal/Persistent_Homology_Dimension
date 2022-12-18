
import numpy                 as np
import gudhi                 as gd



def dgm_square_to_correct_values(dgm):
    """
    Gudhi returns the squared birth and death values, for reasons unknown
    """
    corrected_dgm = []
    for interval in dgm:
        (degree, (birth,death)) = interval
        corrected_dgm.append((degree,(np.sqrt(birth),np.sqrt(death))))
    return corrected_dgm


def test_criticality(center, radius, points):
    """
    tests if the (circum)sphere of given center and radius contains any of the points
    """
    for point in points:
        if np.linalg.norm(point-center)<radius-0.000000001:
            print("Criticality breached !")

def Cayley_Menger_matrix(simplex):
    """
    simplex is a numpy array of shape (n+1,d) representing an n simplex in Rd

    Returns its Cayley-Menger matrix as a numpy array of shape (n+2,n+2)
    """
    n = np.shape(simplex)[0]-1
    M = np.ones((n+2,n+2), dtype= float)
    M[0,0] = 0
    for i in range(1,n+2):
        for j in range(1,n+2):
            M[i,j] = np.linalg.norm(simplex[i-1]-simplex[j-1])**2
    return M


def circumcenter_of_simplex(simplex):
    """
    simplex is a numpy array of shape (n+1,d) representing an n simplex in Rd

    Returns the center of its (smallest) circumsphere as a numpy array of shape (d,)
    """
    # Not as efficient as it could be (we only use the first column of Q, would be easy to fix)
    M = Cayley_Menger_matrix(simplex)
    Q = -2* np.linalg.inv(M)
    partial_Q = Q[0,1:]/np.sum(Q[0,1:])
    center = np.matmul(simplex.T,partial_Q)
    return center


def compute_alpha_complex(X : np.ndarray, precision = 1e-15):
    """
    X should be a numpy array of shape (N,d)

    Returns the alpha complex of X with precisions "precision"

    """
    Acomplex = gd.AlphaComplex(points=X)
    Acomplex.set_float_relative_precision(precision)
    print(f"Precision = {Acomplex.get_float_relative_precision()}")
    return Acomplex

    

def get_critical_pairs(simplex_tree_with_computed_persistence, homology_degree : int, min_len_intervals = 1e-10, min_birth_value = 0.0, essential_part_included = False):
    """
    simplex_tree_with_computed_persistence is a simplex_tree on which compute_persistence() has already been called

    Returns all the pairs of simplices (as pairs of lists of integers) giving birth and death to intervals of length longer than len_intervals
    in degree homology_degree and born after min_birth_value

    If essential_parts_included is True, infinite intervals are also considered (corresponding to pairs similar to ([1,6,12],[]) )

    """
    st = simplex_tree_with_computed_persistence
    # persistence_pairs are pairs of negative and negative simplices giving birth and death to an interval
    # simplices are represented as sets of integers (each corresponding to a point in the points cloud)
    persistence_pairs = st.persistence_pairs()
    critical_pairs_in_correct_degree = []

    for pair in persistence_pairs:
        # only consider homology in the correct degree
        if len(pair[0]) == homology_degree+1:
            # only consider intervals of length longer than len_intervals
            if abs(np.sqrt(st.filtration(pair[0])) - np.sqrt(st.filtration(pair[1])))>=min_len_intervals:
                # only consider intervals born after min_birth_value
                if np.sqrt(st.filtration(pair[0])) >= min_birth_value:
                    # some intervals are infinite
                    if essential_part_included == True or pair[1] != []:
                        critical_pairs_in_correct_degree.append(pair)

    return critical_pairs_in_correct_degree

def get_lengths_intervals(persistence_pairs, simplex_tree) -> list:
    """
    Takes as input a list of persistence pairs and the simplex_tree whose persistence is being considered

    Returns the length of the interval corresponding to each pair
    """

    criticality = []
    for pair in persistence_pairs:
        criticality.append(np.sqrt(simplex_tree.filtration(pair[1])) - np.sqrt(simplex_tree.filtration(pair[0])))
    return criticality