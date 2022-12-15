import numpy                 as np
from math                    import floor, ceil,sqrt
import itertools


from basic_functions import test_criticality, circumcenter_of_simplex

def points_on_J_T_parabolas_example(N : int, epsilon : float):
    """
    Returns 4*floor(N/4) points with step epsilon on the two parabolas (2*N on each) described in the article Relative Neighborhood Graphs and Their Relatives
    
    """
    N = floor(N/4)
    x_1 = []
    for i in range(-N,N,1):
        t_1= i*epsilon
        x_1.append(np.array([1-(t_1**2)/4,0,t_1]))
    x_1 = np.array(x_1)

    x_2 = []
    for i in range(-N,N,1):
        t_2 =  i*epsilon
        x_2.append(np.array([(t_2**2)/4-1,t_2,0]))
    x_2 = np.array(x_2)
    X_total = np.concatenate((x_1,x_2))
    print(f"Number of sampled points = {np.shape(X_total)[0]}")
    return X_total



def regular_simplex(d : int):
    """
    Returns a regular d-simplex in Rd with vertices e1, ... , ed, -t*(e1 +... +ed), where t is chosen so that the simplex is regular

    The simplex is a numpy array of shape (d+1,d)
    """
    last_point = -np.ones((1,d))*(sqrt(float(d)+1)-1)/float(d)
    points = np.concatenate((np.eye(d),last_point), axis = 0 )
    return points


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

def points_on_generalized_critical_configuration(N : int, d : int, theta = np.pi/2):
    """
    Samples points on a generalization of the two interlaced circles

    There are d+1 circular arc in ambient dimension 2*d+1

    Each circular arc is centered on the center of one of the faces of a standard d simplex, and passes through the opposite vertex of the simplex

    They each belong to different planes (that each contain the vector from the center of the corresponding face to the opposite vertex, and a different vector of R^(2d+1) canonical basis)

    Each circular arc describes an angle of theta

    floor(N/(d+1)) points are sampled regularly on each arc
    
    """
    N = floor(N/(d+1))
    print("Sampling points on the generalization of the construction from Relative Neighborhood Graphs and Their Relatives")
    simplex = regular_simplex(d)
    X = []
    for index, point in enumerate(simplex):
        opposite_face = simplex[[i for i in range(len(simplex)) if i!=index]]
        center = np.concatenate((circumcenter_of_simplex(opposite_face),np.zeros(d+1)))
        point_2 = np.concatenate((point,np.zeros(d+1)))
        radius = np.linalg.norm(point_2-center)
        e_1 = point_2-center
        e_1 = e_1/np.linalg.norm(e_1)
        e_2 = np.concatenate((np.zeros(d+index), np.ones(1), np.zeros(2*d+1-(d+index)-1)))
        

        for i in range(N):
            angle = float(i)/N*theta - theta/2
            point_on_circle = (np.cos(angle)*e_1 + np.sin(angle)*e_2)*radius + center
            X.append(point_on_circle)
        
    X = np.array(X)
    print(f"Number of sampled points = {np.shape(X)[0]}")
    
    return X





def points_on_parallel_cubes_in_RD(N : int, d :  int, D : int, length = 10.0) -> np.ndarray:
    """
    Samples N points on two parallel cubes of dimension d and side "length" in R^D, facing each other and at distance 1
    
    d must be <= D-1 ; one of the cubes has coordinate x_D = 0, the other coordinate x_D = 1

    The points are sampled uniformly (floor(N/2) on one cube and N-floor(N/2) on the other)

    Returns an (N,D) numpy array
    """
    print(f"Sampling points on two parallel cubes of dimension {d} in R^{D}")

    # Sampling points on the first cube, in the hyperplane {x_d = 0}
    M_1 = floor(N/2)
    X_1 = np.concatenate((np.random.rand(M_1,d)*length,np.zeros((M_1,D-d))), axis = 1)
    
    M_2 = N - M_1
    X_2 = np.concatenate((np.random.rand(M_2,d)*length,np.zeros((M_2,D-d-1)),np.ones((M_2,1))), axis = 1)
    
    X_total = np.concatenate((np.array(X_1),np.array(X_2)), axis = 0)
    print(f"Number of sampled points = {np.shape(X_total)[0]}")
    
    return X_total

def points_on_regular_grid_on_parallel_cubes_in_RD(epsilon : float, d: int, D : int, length = 10.0) -> np.ndarray:
    """
    Samples points on two parallel cubes of dimension d and side "length" in R^D, facing each other and at distance 1
    
    d must be <= D-1 ; one of the cubes has coordinate x_D = 0, the other coordinate x_D = 1

    The points are sampled on a regular grid of step size epsilon

    The number N of points sampled is 2*ceil(length/epsilon)^d

    Returns an (N,D) numpy array
    """
    print(f"Sampling points on a regular grid of step size {epsilon} on two parallel cubes of dimension {d} in R^{D}")


    M = ceil(length/epsilon)
    # Sampling points on the first cube, in the hyperplane {x_d = 0}

    grid_positions = list(itertools.product(*[range(0,M+1)]*d))
    N = len(grid_positions)
    
    positions = np.array(grid_positions)*epsilon
    X_1 = np.concatenate((positions,np.zeros((N,D-d))), axis = 1)
    X_2 = np.concatenate((positions,np.zeros((N,D-d-1)),np.ones((N,1))), axis = 1)
    
    X_total = np.concatenate((np.array(X_1),np.array(X_2)), axis = 0)
    print(f"Number of sampled points = {np.shape(X_total)[0]}")
    
    return X_total

def points_on_two_segments_in_R3(N : int, x_0 = [2,-2], v_1 = 0.1, v_2 = 0.4, length = 5.0) -> np.ndarray:
    """
    Samples N points on two segments of specified length in R^3 at distance 1 (floor(N/2) on one segment and N-floor(N/2) on the other segment)
    
    One segment is of the shape Conv({(0,0,0),(length,0,0)}) ; the other starts in x_0 (in R^3) and has direction (v_1,v_2,0) (this situation is generic up to rotation and symmetry)
    
    The points are sampled uniformly 

    Returns an (N,d) numpy array
    """
    print(f"Sampling points on two segments in R3")

    # Sampling points on the first segment
    M_1 = floor(N/2)
    X_1 = np.concatenate((np.random.rand(M_1,1)*length,np.zeros((M_1,2))), axis = 1)
    
    # Sampling points on the second segment
    M_2 = N - M_1
    X_2 = np.concatenate((np.random.rand(M_2,1)*length,np.zeros((M_2,1))), axis = 1)
    v = np.array([v_1,v_2])/np.linalg.norm(np.array([v_1,v_2]))
    rotation_matrix = np.array([[v[0], -v[1]],[v[1], v[0]]])
    X_2 = np.matmul(rotation_matrix, X_2.T).T
    X_2 = np.concatenate((X_2+ np.array(x_0),np.ones((M_2,1))), axis = 1)
    
    
    X_total = np.concatenate((np.array(X_1),np.array(X_2)), axis = 0)
    print(f"Number of sampled points = {np.shape(X_total)[0]}")
    
    return X_total





def points_on_lissajous_curve_in_R2(N : int) -> np.ndarray:
    """
    Returns a numpy array of shape (N,2) of points sampled (uniformly ?) on Lissajou's curve
    """
    print("Sampling points on Lissajou's curve in R2")
    
    t = 2*np.pi*np.random.rand(N)
    X = np.transpose(np.concatenate((np.sin(2*t),np.sin(3*t))).reshape(2,N))
    print(f"Number of sampled points = {np.shape(X)[0]}")
    return X

def points_on_torus_in_R3(N : int, r : float, R : float) -> np.ndarray:
    """
    Return a numpy array of shape (N,3) of points NON-UNIFORMLY sampled on the torus in R3 with transverse circle of radius r, 
    and distance from the center of the torus to any center of the transverse circle equal to R.
    """
    print("Sampling points on the torus in R3")

    u = 2*np.pi*np.random.rand(N)
    v = 2*np.pi*np.random.rand(N)
    X = np.transpose(np.vstack([(R+r*np.cos(v))*np.cos(u),(R+r*np.cos(v))*np.sin(u),r*np.sin(v)]))
    print(f"Number of sampled points = {np.shape(X)[0]}")
    return X


if __name__ == "__main__":
    print(points_on_interlaced_2_spheres_in_R5(10,0.1))
    print(points_on_interlaced_2_circles_in_R3(10,0.1))
    print(points_on_parallel_cubes_in_RD(20,2,3))
    print(points_on_lissajous_curve_in_R2(10))
    print(points_on_torus_in_R3(10,1,4))
    print(points_on_regular_grid_on_parallel_cubes_in_RD(1, 2, 3, length = 3) )
    print(points_on_two_segments_in_R3(10))
    points_on_generalized_critical_configuration(10,3)



