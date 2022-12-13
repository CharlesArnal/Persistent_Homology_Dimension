
import numpy                 as np
import matplotlib.pyplot     as plt
import gudhi                 as gd

from points_clouds_generation import points_on_parallel_hyperplanes_in_Rd, points_on_interlaced_2_circles_in_R3



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



def plot_crit_points_3D_points_cloud(X : np.ndarray, dimension: int, positive_or_negative = "positive", len_intervals = 1e-10, precision = 1e-15):
    """
    X should be a numpy array of shape (N,3)

    The function plots the centers of the critical positive (resp, negative) simplices of specified dimension for the alpha filtration of X

    If positive_or_negative = "positive", these simplices will correspond to birth point of intervals in the homology of degree dimension

    If positive_or_negative = "positive", they will correspond to death point of intervals in the homology of degree (dimension - 1)

    A simplex is considered critical if it gives birth or death to an interval of length at least len_intervals
    
    """

    print("Plotting the critical points of the alpha-filtration of the points cloud :")

    Acomplex = gd.AlphaComplex(points=X)
    Acomplex.set_float_relative_precision(precision)
    print(f"Precision = {Acomplex.get_float_relative_precision()}")

    print("Computing persistence diagram")
    st = Acomplex.create_simplex_tree()
    st.compute_persistence()
    # persistence_pairs are pairs of negative and negative simplices giving birth and death to an interval
    # simplices are represented as sets of integers (each corresponding to a point in the points cloud)
    persistence_pairs = st.persistence_pairs()
    critical_simpl_in_correct_degree = []
    # the criticality of a simplex is simply the length of the interval to which it gives birth (respectively death)
    criticality = []

    if positive_or_negative == "positive":
        homology_degree = dimension
    else:
        homology_degree = dimension - 1

    for pair in persistence_pairs:
        # only consider homology in the correct degree
        if len(pair[0]) == homology_degree+1:
            # only consider intervals of length longer than len_intervals
            if abs(st.filtration(pair[0]) - st.filtration(pair[1]))>len_intervals:
                if positive_or_negative == "positive":
                    critical_simpl_in_correct_degree.append(pair[0])
                    criticality.append(st.filtration(pair[0]) - st.filtration(pair[1]))
                else:
                    # some intervals are infinite
                    if pair[1] != []:
                        critical_simpl_in_correct_degree.append(pair[1])
                        criticality.append(st.filtration(pair[0]) - st.filtration(pair[1]))


    critical_points = np.array([circumcenter_of_simplex(np.array([Acomplex.get_point(point) for point in simplex])) for simplex in critical_simpl_in_correct_degree])
    criticality = np.array(criticality)


    fig = plt.figure( num = f"Critical points for {positive_or_negative} simplices of dimension {dimension}")
    ax = plt.axes(projection='3d')
    print(np.shape(critical_points))
    ax.scatter3D(critical_points[:,0],critical_points[:,1],critical_points[:,2],c = criticality)
    ax.scatter3D(X[:,0],X[:,1],X[:,2])
    plt.show()


if __name__ == "__main__":
    np.random.seed(0)
    N = 800
    d = 3
    epsilon = 0.5
    X = points_on_parallel_hyperplanes_in_Rd(N, d)
    #X = points_on_interlaced_2_circles_in_R3(N, epsilon)


    positive_or_negative = "positive"
    dimension = 1
    plot_crit_points_3D_points_cloud(X, dimension, positive_or_negative)
    
    
