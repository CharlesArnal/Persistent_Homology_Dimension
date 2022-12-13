import numpy                 as np
from math                    import floor


def points_on_parallel_hyperplanes_in_Rd(N : int, d: int, length = 10.0) -> np.ndarray:
    """
    Samples N points on two parallel cubes of dimension d-1 and side "length" in Rd, facing each other and at distance 1

    The points are sampled uniformly (floor(N/2) on one cube and N-floor(N/2) on the other)

    Returns an (N,d) numpy array
    """
    print(f"Sampling points on two parallel cubes of codimension 1 in R^{d}")

    # Sampling points on the first cube, in the hyperplane {x_d = 0}
    M_1 = floor(N/2)
    X_1 = np.concatenate((np.random.rand(M_1,d-1)*length,np.zeros((M_1,1))), axis = 1)
    
    M_2 = N - M_1
    X_2 = np.concatenate((np.random.rand(M_2,d-1)*length,np.ones((M_2,1))), axis = 1)
    
    X_total = np.concatenate((np.array(X_1),np.array(X_2)), axis = 0)
    print(f"Number of sampled points = {np.shape(X_total)[0]}")
    
    return X_total



def points_on_interlaced_2_circles_in_R3(N : int, epsilon: float) -> np.ndarray:
    """
    Samples N points on two interlaced circles of radius 1 in R^3 (floor(N/2) on one and N-floor(N/2) on the other)

    epsilon is a concentration parameter

    The points are first uniformly sampled on the circles, then the points at distance more than epsilon
    from the center of the circle on which they are not are rejected (points keep being sampled until we have enough of them)

    Returns an (N,3) numpy array
    """
    # Rmk : the sampling method is not very efficient, but the sampling isn't the computational bottleneck of the experiments anyway

    
    print("Sampling points on two interlaced cicles in R3")

    # center first circle = [1,0,0], plan x1,x3
    center_1 = np.array([1,0,0])
    # center second circle = [0,0,0], plan x1,x2
    center_2 = np.zeros(3)

    # Sampling points on the first circle
    M_1 = floor(N/2)
    X_1 = []
    while len(X_1)<M_1:
        # points uniformly sampled on a circle centered in center_1
        candidate = np.concatenate((np.random.randn(1),np.zeros(1),np.random.randn(1)))
        candidate =  candidate/np.linalg.norm(candidate)
        candidate = candidate + center_1
        # only keep those at distance less than epsilon from center_2
        if np.linalg.norm(candidate-center_2)<epsilon:
            X_1.append(candidate)

    
    # Sampling points on the second circle
    M_2 = N - M_1
    X_2 = []
    while len(X_2)<M_2:
        # points uniformly sampled on a sphere centered in center_2
        candidate = np.concatenate((np.random.randn(2),np.zeros(1)))
        candidate =  candidate/np.linalg.norm(candidate)
        candidate = candidate + center_2
        # only keep those at distance less than epsilon from center_1
        if np.linalg.norm(candidate-center_1)<epsilon:
            X_2.append(candidate)


    X_total = np.concatenate((np.array(X_1),np.array(X_2)), axis = 0)
    print(f"Number of sampled points = {np.shape(X_total)[0]}")
    return X_total



def points_on_interlaced_2_spheres_in_R5(N : int, epsilon : float) -> np.ndarray:
    """
    Samples N points on two interlaced 2-spheres of radius 1 in R^5 (floor(N/2) on one and N-floor(N/2) on the other)

    epsilon is a concentration parameter

    The points are first uniformly sampled on the spheres, then the points at distance more than epsilon
    from the center of the sphere on which they are not are rejected (points keep being sampled until we have enough of them)

    Returns an (N,5) numpy array
    """
    # Rmk : the sampling method is not very efficient, but the sampling isn't the computational bottleneck of the experiments anyway

    print("Sampling points on two interlaced 2-spheres in R5")
    
    # center first sphere = [1,0,0,0,0], plan x1,x2,x3
    center_1 = np.array([1,0,0,0,0])
    # center second sphere = [0,0,0,0,0], plan x1, x4, x5
    center_2 = np.zeros(5)

    # Sampling points on the first sphere
    M_1 = floor(N/2)
    X_1 = []
    while len(X_1)<M_1:
        # points uniformly sampled on a sphere centered in center_1
        candidate = np.concatenate((np.random.randn(3),np.zeros(2)))
        candidate =  candidate/np.linalg.norm(candidate)
        candidate = candidate + center_1
        # only keep those at distance less than epsilon from center_2
        if np.linalg.norm(candidate-center_2)<epsilon:
            X_1.append(candidate)
    
    # Sampling points on the second sphere
    M_2 = N - M_1
    X_2 = []
    while len(X_2)<M_2:
        # points uniformly sampled on a sphere centered in center_2
        candidate = np.concatenate((np.random.randn(1),np.zeros(2),np.random.randn(2)))
        candidate =  candidate/np.linalg.norm(candidate)
        candidate = candidate + center_2
        # only keep those at distance less than epsilon from center_1
        if np.linalg.norm(candidate-center_1)<epsilon:
            X_2.append(candidate)


    X_total = np.concatenate((np.array(X_1),np.array(X_2)), axis = 0)
    print(f"Number of sampled points = {np.shape(X_total)[0]}")
    return X_total


if __name__ == "__main__":
    print(points_on_interlaced_2_spheres_in_R5(10,0.1))
    print(points_on_interlaced_2_circles_in_R3(10,0.1))
    print(points_on_parallel_hyperplanes_in_Rd(20,3))