import numpy                 as np
import gudhi                 as gd
import gudhi.wasserstein
import matplotlib.pyplot     as plt

from points_clouds_generation import points_on_interlaced_2_spheres_in_R5


def compute_total_persistence_alpha(X : np.ndarray, homology_degree : int, order : int, precision = 1e-05, plot = False) -> float:
    """
    X must be an (N,d) numpy array (representing a poins cloud)

    Returns the total persistence of the homology of degree "homology_degree" of the alpha filtration of X
    w.r.t the Wasserstein distance of order "order"

    If plot = True, the function also plots the persistence diagram of given degree
    """

    print("Computing persistence diagram")

    Acomplex = gd.AlphaComplex(points=X)
    Acomplex.set_float_relative_precision(precision)
    print(f"Precision = {Acomplex.get_float_relative_precision()}")
    st = Acomplex.create_simplex_tree()
    dgm = st.persistence()

    dgm_in_correct_degree = []

    for point in dgm :
        if point[0] == homology_degree:
            dgm_in_correct_degree.append(point[1])
    dgm_in_correct_degree = np.array(dgm_in_correct_degree)


    if plot:
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


    print("Computing Wasserstein distance")
    return gudhi.wasserstein.wasserstein_distance(dgm_in_correct_degree,np.array([]), order = order, keep_essential_parts=True)



if __name__ == "__main__":
    N = 80
    homology_degree = 1
    order = 2
    epsilon = 1
    plot = True
    precision = 1e-8
    X = points_on_interlaced_2_spheres_in_R5(N, epsilon)
    print(f"Wasserstein distance = {compute_total_persistence_alpha(X, homology_degree, order, precision, plot)}")

