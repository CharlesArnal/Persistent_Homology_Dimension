import numpy                 as np

from points_clouds_generation import points_on_J_T_parabolas_example, points_on_generalized_critical_configuration, points_on_interlaced_2_spheres_in_R5,\
 points_on_interlaced_2_circles_in_R3, points_on_parallel_cubes_in_RD, points_on_regular_grid_on_parallel_cubes_in_RD, points_on_two_segments_in_R3, \
    points_on_torus_in_R3, points_on_lissajous_curve_in_R2
import matplotlib.pyplot     as plt

from explore_critical_simplices import plot_crit_points_3D_points_cloud, plot_crit_circumspheres_3D_points_cloud
from persistence_computations import  plot_diagrams, get_number_of_intervals, get_mean_and_std_of_intervals_lengths
from basic_functions import compute_alpha_complex, dgm_square_to_correct_values



np.random.seed(0)



precision = 1e-15


positive_or_negative = "positive"

experience_num = 8

def get_number_and_mean_length_intervals(X, max_degree):
    Acomplex = compute_alpha_complex(X,precision)
    print("Computing simplex tree")
    st = Acomplex.create_simplex_tree()
    print("Computing persistence diagram")
    st.compute_persistence()


    print("Computing stats")
    number_of_intervals = {}
    mean_length_of_intervals = {}

    for degree in range(max_degree+1):
        number_of_intervals[degree ] = get_number_of_intervals(st, degree, min_len_intervals = 1e-14, min_birth_value = 0.2)
        mean_length_of_intervals[degree], _ = get_mean_and_std_of_intervals_lengths(st, degree, min_len_intervals = 1e-10, min_birth_value = 0.2)
    
    return number_of_intervals, mean_length_of_intervals



if experience_num == 1:
    number_of_points = [10,30,60,100,200,400,1000]
    max_degree = 2
    number_of_intervals = {}
    mean_lengths_of_intervals = {}
    # needed for plotting purposes
    adapted_number_of_points = {}
    real_number_of_points = []
    for degree in range(max_degree+1):
        number_of_intervals[degree] = []
        mean_lengths_of_intervals[degree] = []
        adapted_number_of_points[degree] = []
    for N in number_of_points:
        print(f"N = {N}")
        X = points_on_interlaced_2_circles_in_R3(N, 0.5)
        num_int, mean_len =  get_number_and_mean_length_intervals(X, max_degree)
        real_number_of_points.append(len(X))
        for degree in range(max_degree+1):
            # to avoid taking the log of zero
            if num_int[degree] > 0:
                adapted_number_of_points[degree].append(len(X))
                number_of_intervals[degree].append(num_int[degree])
                mean_lengths_of_intervals[degree].append(mean_len[degree])

    #fig = plt.figure( num = f"Number and mean lengths of intervals for the interlaced circles in R3")
    for degree in range(max_degree+1):
        """"
        plt.subplot(2*max_degree+2,1,2*degree+1)
        plt.plot(adapted_number_of_points[degree], number_of_intervals[degree])
        plt.xscale('log')
        plt.yscale('log')
        plt.subplot(2*max_degree+2,1,2*degree+2)
        plt.plot(adapted_number_of_points[degree],mean_lengths_of_intervals[degree])
        plt.xscale('log')
        plt.yscale('log')
        """
        print(f"For degree {degree} and N = {real_number_of_points}, log(number of intervals)/log(N) : {np.log(np.array(number_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
        print(f"For degree {degree} and N = {real_number_of_points}, log(mean length of intervals)/log(N) : {np.log(np.array(mean_lengths_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
    
        
    #plt.show()       

"""
For degree 0 and N = [10, 30, 60, 100, 200, 400, 1000], log(number of intervals)/log(N) : []
For degree 0 and N = [10, 30, 60, 100, 200, 400, 1000], log(mean length of intervals)/log(N) : []
For degree 1 and N = [10, 30, 60, 100, 200, 400, 1000], log(number of intervals)/log(N) : [1.2787536  1.43561263 1.48158222 1.49226366 1.50845412 1.51773733  1.50549416]
For degree 1 and N = [10, 30, 60, 100, 200, 400, 1000], log(mean length of intervals)/log(N) : [-2.5188373  -2.49647941 -2.41831246 -2.37696369 -2.35488381 -2.3378469  -2.29558413]
For degree 2 and N = [10, 30, 60, 100, 200, 400, 1000], log(number of intervals)/log(N) : [1.04139269 1.39507697 1.45654066 1.47663817 1.49983797 1.51232272  1.50228105]
For degree 2 and N = [10, 30, 60, 100, 200, 400, 1000], log(mean length of intervals)/log(N) : [-2.63636335 -2.6280188  -2.45951809 -2.42015518 -2.38238373 -2.35856904  -2.30980961]
"""


if experience_num == 2:
    number_of_points = [10,30,60,100,150,200,300,400]
    max_degree = 4
    number_of_intervals = {}
    mean_lengths_of_intervals = {}
    # needed for plotting purposes
    adapted_number_of_points = {}
    for degree in range(max_degree+1):
        number_of_intervals[degree] = []
        mean_lengths_of_intervals[degree] = []
        adapted_number_of_points[degree] = []
    for N in number_of_points:
        print(f"N = {N}")
        X =  points_on_interlaced_2_spheres_in_R5(N, 0.5)
        num_int, mean_len =  get_number_and_mean_length_intervals(X, max_degree)
        for degree in range(max_degree+1):
            # to avoid taking the log of zero
            if num_int[degree] > 0:
                adapted_number_of_points[degree].append(N)
                number_of_intervals[degree].append(num_int[degree])
                mean_lengths_of_intervals[degree].append(mean_len[degree])

    for degree in range(max_degree+1):

        print(f"For degree {degree} and N = {number_of_points}, log(number of intervals)/log(N) : {np.log(np.array(number_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
        print(f"For degree {degree} and N = {number_of_points}, log(mean length of intervals)/log(N) : {np.log(np.array(mean_lengths_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
        
   
"""
For degree 0 and N = [10, 30, 60, 100, 150, 200, 300, 400], log(number of intervals)/log(N) : []
For degree 0 and N = [10, 30, 60, 100, 150, 200, 300, 400], log(mean length of intervals)/log(N) : []
For degree 1 and N = [10, 30, 60, 100, 150, 200, 300, 400], log(number of intervals)/log(N) : [1.32221929 1.56360228 1.64514244 1.67493004 1.67055615 1.66677966
 1.68114998 1.65799575]
For degree 1 and N = [10, 30, 60, 100, 150, 200, 300, 400], log(mean length of intervals)/log(N) : [-2.17179628 -2.0066756  -1.90321642 -1.82280074 -1.7442186  -1.75620573
 -1.69172342 -1.65594013]
For degree 2 and N = [10, 30, 60, 100, 150, 200, 300, 400], log(number of intervals)/log(N) : [1.38021124 1.79620495 1.83447622 1.86735992 1.85403449 1.83857999
 1.84977631 1.82688024]
For degree 2 and N = [10, 30, 60, 100, 150, 200, 300, 400], log(mean length of intervals)/log(N) : [-2.53576605 -2.24540641 -1.99967408 -1.90499472 -1.79922796 -1.79750442
 -1.7460747  -1.70479741]
For degree 3 and N = [10, 30, 60, 100, 150, 200, 300, 400], log(number of intervals)/log(N) : [1.         1.71808447 1.75066096 1.82437511 1.81913862 1.80269145
 1.82654154 1.81381142]
For degree 3 and N = [10, 30, 60, 100, 150, 200, 300, 400], log(mean length of intervals)/log(N) : [-2.79792295 -2.50180897 -2.11128213 -2.00610805 -1.86073301 -1.84076261
 -1.80897997 -1.75797876]
For degree 4 and N = [10, 30, 60, 100, 150, 200, 300, 400], log(number of intervals)/log(N) : [1.3197224  1.39308805 1.54618485 1.56723197 1.56039142 1.6113072
 1.61830763]
For degree 4 and N = [10, 30, 60, 100, 150, 200, 300, 400], log(mean length of intervals)/log(N) : [-2.62985338 -2.2338662  -2.10698668 -1.92459814 -1.88307015 -1.87083234
 -1.8088941 ]
"""





if experience_num == 3:
    number_of_points = [10,30,60,100,200,400,1000, 1500]
    max_degree = 2
    number_of_intervals = {}
    mean_lengths_of_intervals = {}
    # needed for plotting purposes
    adapted_number_of_points = {}
    for degree in range(max_degree+1):
        number_of_intervals[degree] = []
        mean_lengths_of_intervals[degree] = []
        adapted_number_of_points[degree] = []
    for N in number_of_points:
        print(f"N = {N}")
        X = points_on_J_T_parabolas_example(N, 0.05)
        num_int, mean_len =  get_number_and_mean_length_intervals(X, max_degree)
        for degree in range(max_degree+1):
            # to avoid taking the log of zero
            if num_int[degree] > 0:
                adapted_number_of_points[degree].append(N)
                number_of_intervals[degree].append(num_int[degree])
                mean_lengths_of_intervals[degree].append(mean_len[degree])

    for degree in range(max_degree+1):

        print(f"For degree {degree} and N = {number_of_points}, log(number of intervals)/log(N) : {np.log(np.array(number_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
        print(f"For degree {degree} and N = {number_of_points}, log(mean length of intervals)/log(N) : {np.log(np.array(mean_lengths_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
     

"""
For degree 0 and N = [10, 30, 60, 100, 200, 400, 1000, 1500], log(number of intervals)/log(N) : []
For degree 0 and N = [10, 30, 60, 100, 200, 400, 1000, 1500], log(mean length of intervals)/log(N) : []
For degree 1 and N = [10, 30, 60, 100, 200, 400, 1000, 1500], log(number of intervals)/log(N) : [1.17609126 1.55033624 1.53300324 1.36295582 1.18464846 1.04759753
 0.90863721 0.85825986]
For degree 1 and N = [10, 30, 60, 100, 200, 400, 1000, 1500], log(mean length of intervals)/log(N) : [-4.10836593 -2.83180094 -2.50141458 -2.22394673 -1.93300108 -1.70937391
 -1.48263116 -1.40043   ]
For degree 2 and N = [10, 30, 60, 100, 200, 400, 1000, 1500], log(number of intervals)/log(N) : [0.95424251 1.50826257 1.50990832 1.34242268 1.16680155 1.03181532
 0.89494845 0.84533004]
For degree 2 and N = [10, 30, 60, 100, 200, 400, 1000, 1500], log(mean length of intervals)/log(N) : [-4.10861884 -2.83296774 -2.50635397 -2.22833822 -1.93681806 -1.71274931
 -1.48555881 -1.40319534]
"""



if experience_num == 4:
    number_of_points = [10,30,60,100,200,400,1000, 1500, 2000, 3000]
    max_degree = 2
    number_of_intervals = {}
    mean_lengths_of_intervals = {}
    # needed for plotting purposes
    adapted_number_of_points = {}
    for degree in range(max_degree+1):
        number_of_intervals[degree] = []
        mean_lengths_of_intervals[degree] = []
        adapted_number_of_points[degree] = []
    for N in number_of_points:
        print(f"N = {N}")
        X = points_on_torus_in_R3(N, 1, 5)
        num_int, mean_len =  get_number_and_mean_length_intervals(X, max_degree)
        for degree in range(max_degree+1):
            # to avoid taking the log of zero
            if num_int[degree] > 0:
                adapted_number_of_points[degree].append(N)
                number_of_intervals[degree].append(num_int[degree])
                mean_lengths_of_intervals[degree].append(mean_len[degree])

    for degree in range(max_degree+1):

        print(f"For degree {degree} and N = {number_of_points}, log(number of intervals)/log(N) : {np.log(np.array(number_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
        print(f"For degree {degree} and N = {number_of_points}, log(mean length of intervals)/log(N) : {np.log(np.array(mean_lengths_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
     

"""
For degree 0 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000], log(number of intervals)/log(N) : []
For degree 0 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000], log(mean length of intervals)/log(N) : []
For degree 1 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000], log(number of intervals)/log(N) : [0.30103    0.92187952 1.0232785  1.05360498 1.03907168 1.0166646
 0.9819869  0.95096804 0.91510089 0.84808712]
For degree 1 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000], log(mean length of intervals)/log(N) : [-1.18234389 -0.60113062 -0.59478106 -0.55796989 -0.54000095 -0.49438944
 -0.46904751 -0.45411388 -0.4643175  -0.52956279]
For degree 2 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000], log(number of intervals)/log(N) : [0.20379505 0.7861761  0.88542601 0.88544109 0.88431089 0.87815909
 0.88856643 0.89406565 0.89612727]
For degree 2 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000], log(mean length of intervals)/log(N) : [-1.4457377  -1.15126406 -1.01639019 -0.90827126 -0.86619857 -0.94465532
 -0.95533868 -0.95835318 -0.94114086]
"""




if experience_num == 5:
    number_of_points = [10,30,60,100,200,400,1000, 1500, 2000, 3000, 4000]
    max_degree = 3
    number_of_intervals = {}
    mean_lengths_of_intervals = {}
    # needed for plotting purposes
    adapted_number_of_points = {}
    for degree in range(max_degree+1):
        number_of_intervals[degree] = []
        mean_lengths_of_intervals[degree] = []
        adapted_number_of_points[degree] = []
    for N in number_of_points:
        print(f"N = {N}")
        X = points_on_parallel_cubes_in_RD(N,2,3,length = 2)
        num_int, mean_len =  get_number_and_mean_length_intervals(X, max_degree)
        for degree in range(max_degree+1):
            # to avoid taking the log of zero
            if num_int[degree] > 0:
                adapted_number_of_points[degree].append(N)
                number_of_intervals[degree].append(num_int[degree])
                mean_lengths_of_intervals[degree].append(mean_len[degree])

    for degree in range(max_degree+1):

        print(f"For degree {degree} and N = {number_of_points}, log(number of intervals)/log(N) : {np.log(np.array(number_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
        print(f"For degree {degree} and N = {number_of_points}, log(mean length of intervals)/log(N) : {np.log(np.array(mean_lengths_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")


"""
For degree 0 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(number of intervals)/log(N) : []
For degree 0 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(mean length of intervals)/log(N) : []
For degree 1 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(number of intervals)/log(N) : [1.         0.95792633 1.0792547  1.05190186 1.00466046 0.98517371
 0.99088853 0.99173306 0.99720775 0.99958297 0.99771788]
For degree 1 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(mean length of intervals)/log(N) : [-1.61866814 -1.21559224 -1.07336117 -1.06302283 -1.06952234 -1.13130785
 -1.12185975 -1.11247765 -1.11662292 -1.11135213 -1.10752116]
For degree 2 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(number of intervals)/log(N) : [0.30103    0.52680255 0.93510484 0.92254902 0.95788407 0.94747349
 0.98394101 0.98609865 0.98903004 0.99403092 0.99131507]
For degree 2 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(mean length of intervals)/log(N) : [-2.30318124 -2.37938876 -1.39687955 -1.26139722 -1.21427997 -1.15204891
 -1.14831493 -1.13641477 -1.13418344 -1.12238956 -1.11248888]
"""



if experience_num == 6:
    number_of_points = [30,100,200,400,1000,  2000, 3000, 4000]
    max_degree = 3
    number_of_intervals = {}
    mean_lengths_of_intervals = {}
    # needed for plotting purposes
    adapted_number_of_points = {}
    for degree in range(max_degree+1):
        number_of_intervals[degree] = []
        mean_lengths_of_intervals[degree] = []
        adapted_number_of_points[degree] = []
    for N in number_of_points:
        print(f"N = {N}")
        length = 2
        X = points_on_parallel_cubes_in_RD(N,3,6,length = length)
        num_int, mean_len =  get_number_and_mean_length_intervals(X, max_degree)
        for degree in range(max_degree+1):
            # to avoid taking the log of zero
            if num_int[degree] > 0:
                adapted_number_of_points[degree].append(N)
                number_of_intervals[degree].append(num_int[degree])
                mean_lengths_of_intervals[degree].append(mean_len[degree])

    for degree in range(max_degree+1):
        epsilon = 2*(np.array(adapted_number_of_points[degree], float)/2)**(-1.0/3.0)
        print(f"epsilon {epsilon}")

        print(f"For degree {degree} and N = {number_of_points}, log(number of intervals)/log(N) ")
        my_string = " & ".join([ f'{result:.2f}' for result in np.log(np.array(number_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))])
        print(my_string)
        print(f"For degree {degree} and N = {number_of_points}, log(mean length of intervals)/log(N) : ")
        my_string = " & ".join([ f'{result:.2f}' for result in np.log(np.array(mean_lengths_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))])
        print(my_string)
        print(f"For degree {degree} and N = {number_of_points}, number of intervals/N : ")
        my_string = " & ".join([ f'{result:.2f}' for result in np.array(number_of_intervals[degree], float)/np.array(adapted_number_of_points[degree], float)])
        print(my_string)
        print(f"For degree {degree} and N = {number_of_points}, mean length of intervals/epsilon^2*200 : ")
        my_string = " & ".join([ f'{result:.2f}' for result in np.array(mean_lengths_of_intervals[degree], float)*epsilon**(-2.0)*200])
        print(my_string)
        print("----\n")

"""
For degree 0 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(number of intervals)/log(N) : []
For degree 0 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(mean length of intervals)/log(N) : []
For degree 1 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(number of intervals)/log(N) : [1.07918125 1.30620213 1.25001541 1.22512455 1.22902138 1.18529884
 1.11081282 1.09215699 1.08235054 1.07509206 1.07353482]
For degree 1 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(mean length of intervals)/log(N) : [-1.5516813  -1.07534616 -0.98977609 -0.89376746 -0.84838711 -0.83538036
 -0.89691577 -0.89767584 -0.89117548 -0.88321346 -0.8766894 ]
For degree 2 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(number of intervals)/log(N) : [0.69897    1.253288   1.22051581 1.22357902 1.25295775 1.23549951
 1.20137267 1.18917207 1.18026346 1.15975788 1.16076043]
For degree 2 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(mean length of intervals)/log(N) : [-2.63318514 -1.46125856 -1.3127155  -1.18087298 -1.09283801 -1.01230632
 -0.93358439 -0.91649513 -0.93699585 -0.93530442 -0.93295264]
For degree 3 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(number of intervals)/log(N) : [0.         0.81518019 0.90700038 0.99561304 1.03362316 1.06683718
 1.06885201 1.0765732  1.08090031 1.07190234 1.07993028]
For degree 3 and N = [10, 30, 60, 100, 200, 400, 1000, 1500, 2000, 3000, 4000], log(mean length of intervals)/log(N) : [-4.2157405  -1.98227359 -1.45998899 -1.3716474  -1.25033011 -1.17735082
 -1.08475902 -1.05366603 -1.04438778 -1.01060827 -0.99851388]
"""

if experience_num == 7:
    number_of_points = [10,30,60,100,200,400,1000, 1500, 2000, 3000, 4000]
    max_degree = 2
    number_of_intervals = {}
    mean_lengths_of_intervals = {}
    # needed for plotting purposes
    adapted_number_of_points = {}
    real_number_of_points = []
    for degree in range(max_degree+1):
        number_of_intervals[degree] = []
        mean_lengths_of_intervals[degree] = []
        adapted_number_of_points[degree] = []
    for N in number_of_points:
        print(f"N = {N}")
        epsilon = np.sqrt(8.0/float(N))
        X = points_on_regular_grid_on_parallel_cubes_in_RD(epsilon,2,3,length = 2)
        num_int, mean_len =  get_number_and_mean_length_intervals(X, max_degree)
        real_number_of_points.append(len(X))
        for degree in range(max_degree+1):
            # to avoid taking the log of zero
            if num_int[degree] > 0:
                adapted_number_of_points[degree].append(len(X))
                number_of_intervals[degree].append(num_int[degree])
                mean_lengths_of_intervals[degree].append(mean_len[degree])

    for degree in range(max_degree+1):

        print(f"For degree {degree} and N = {real_number_of_points}, log(number of intervals)/log(N) : {np.log(np.array(number_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
        print(f"For degree {degree} and N = {real_number_of_points}, log(mean length of intervals)/log(N) : {np.log(np.array(mean_lengths_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")


"""
    For degree 0 and N = [32, 50, 98, 162, 242, 512, 1152, 1682, 2178, 3200, 4232], log(number of intervals)/log(N) : []
For degree 0 and N = [32, 50, 98, 162, 242, 512, 1152, 1682, 2178, 3200, 4232], log(mean length of intervals)/log(N) : []
For degree 1 and N = [32, 50, 98, 162, 242, 512, 1152, 1682, 2178, 3200, 4232], log(number of intervals)/log(N) : [1.00887882 1.02896933 0.84432464 0.86131571 0.87220733 0.88826149
 0.90142436 0.90652109 0.90969929 0.91404034 0.91693603]
For degree 1 and N = [32, 50, 98, 162, 242, 512, 1152, 1682, 2178, 3200, 4232], log(mean length of intervals)/log(N) : [-0.49689918 -0.62127218 -0.74874635 -0.77275163 -0.84078626 -0.85011291
 -0.88188075 -0.89144274 -0.89885609 -0.90620388 -0.91030068]
For degree 2 and N = [32, 50, 98, 162, 242, 512, 1152, 1682, 2178, 3200, 4232], log(number of intervals)/log(N) : [0.633985   0.70873528 0.78158001 0.8174554  0.83899115 0.86819791
 0.88959593 0.89723254 0.90181179 0.90784395 0.9117285 ]
For degree 2 and N = [32, 50, 98, 162, 242, 512, 1152, 1682, 2178, 3200, 4232], log(mean length of intervals)/log(N) : [-0.5769286  -0.73495188 -0.76161323 -0.78003984 -0.84429089 -0.8516846
 -0.88244369 -0.89179985 -0.89911526 -0.90636864 -0.9104202 ]
"""


if experience_num == 8:
    number_of_points = [20,50,100, 200]
    max_degree = 4
    number_of_intervals = {}
    mean_lengths_of_intervals = {}
    # needed for plotting purposes
    adapted_number_of_points = {}
    real_number_of_points = []
    for degree in range(max_degree+1):
        number_of_intervals[degree] = []
        mean_lengths_of_intervals[degree] = []
        adapted_number_of_points[degree] = []
    for N in number_of_points:
        print(f"N = {N}")
        X = points_on_generalized_critical_configuration(N,2)
        num_int, mean_len =  get_number_and_mean_length_intervals(X, max_degree)
        real_number_of_points.append(len(X))
        for degree in range(max_degree+1):
            # to avoid taking the log of zero
            if num_int[degree] > 0:
                adapted_number_of_points[degree].append(len(X))
                number_of_intervals[degree].append(num_int[degree])
                mean_lengths_of_intervals[degree].append(mean_len[degree])

    for degree in range(max_degree+1):

        print(f"For degree {degree} and N = {number_of_points}, log(number of intervals)/log(N) ")
        my_string = " & ".join([ f'{result:.2f}' for result in np.log(np.array(number_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))])
        print(my_string)
        print(f"For degree {degree} and N = {number_of_points}, log(mean length of intervals)/log(N) : ")
        my_string = " & ".join([ f'{result:.2f}' for result in np.log(np.array(mean_lengths_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))])
        print(my_string)




"""
For degree 0 and N = [9, 18, 39, 60, 78, 99], log(number of intervals)/log(N) : []
For degree 0 and N = [9, 18, 39, 60, 78, 99], log(mean length of intervals)/log(N) : []
For degree 1 and N = [9, 18, 39, 60, 78, 99], log(number of intervals)/log(N) : [1.46497352 1.6134392  1.67710721 1.67589983 1.67599503 1.66823797]
For degree 1 and N = [9, 18, 39, 60, 78, 99], log(mean length of intervals)/log(N) : [-1.77020395 -1.86099788 -1.92629626 -1.93627226 -1.94187581 -1.93752509]
For degree 2 and N = [9, 18, 39, 60, 78, 99], log(number of intervals)/log(N) : [1.65553681 1.96044909 2.08700619 2.09358803 2.09419587 2.0892793 ]
For degree 2 and N = [9, 18, 39, 60, 78, 99], log(mean length of intervals)/log(N) : [-2.17169241 -2.21251024 -2.27181128 -2.27486723 -2.27523445 -2.26991577]
For degree 3 and N = [9, 18, 39, 60, 78, 99], log(number of intervals)/log(N) : [1.51655163 1.99893052 2.17924027 2.19318829 2.19544926 2.19454184]
For degree 3 and N = [9, 18, 39, 60, 78, 99], log(mean length of intervals)/log(N) : [-2.40433556 -2.38465994 -2.39242267 -2.36835255 -2.35520857 -2.34278962]
For degree 4 and N = [9, 18, 39, 60, 78, 99], log(number of intervals)/log(N) : [0.94639463 1.66770297 1.95077165 1.99699646 2.01743008 2.02767278]
For degree 4 and N = [9, 18, 39, 60, 78, 99], log(mean length of intervals)/log(N) : [-2.46287882 -2.39743742 -2.39718022 -2.37078936 -2.35976474 -2.34671655]

"""


if experience_num == 9:
    number_of_points = [10,20,30,40,50, 60, 80, 100]
    max_degree = 6
    number_of_intervals = {}
    mean_lengths_of_intervals = {}
    # needed for plotting purposes
    adapted_number_of_points = {}
    real_number_of_points = []
    for degree in range(max_degree+1):
        number_of_intervals[degree] = []
        mean_lengths_of_intervals[degree] = []
        adapted_number_of_points[degree] = []
    for N in number_of_points:
        print(f"N = {N}")
        X = points_on_generalized_critical_configuration(N,3)
        num_int, mean_len =  get_number_and_mean_length_intervals(X, max_degree)
        real_number_of_points.append(len(X))
        for degree in range(max_degree+1):
            # to avoid taking the log of zero
            if num_int[degree] > 0:
                adapted_number_of_points[degree].append(len(X))
                number_of_intervals[degree].append(num_int[degree])
                mean_lengths_of_intervals[degree].append(mean_len[degree])

    for degree in range(max_degree+1):

        print(f"For degree {degree} and N = {real_number_of_points}, log(number of intervals)/log(N) : {np.log(np.array(number_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")
        print(f"For degree {degree} and N = {real_number_of_points}, log(mean length of intervals)/log(N) : {np.log(np.array(mean_lengths_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))}")


"""
For degree 0 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(number of intervals)/log(N) : []
For degree 0 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(mean length of intervals)/log(N) : []
For degree 1 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(number of intervals)/log(N) : [1.46497352 1.6134392  1.67205096 1.67710721 1.68439229 1.67589983
 1.67599503 1.66823797]
For degree 1 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(mean length of intervals)/log(N) : [-1.77020395 -1.86099788 -1.91800856 -1.92629626 -1.93975642 -1.93627226
 -1.94187581 -1.93752509]
For degree 2 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(number of intervals)/log(N) : [1.65553681 1.96044909 2.06403508 2.08700619 2.09036921 2.09358803
 2.09419587 2.0892793 ]
For degree 2 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(mean length of intervals)/log(N) : [-2.17169241 -2.21251024 -2.25984506 -2.27181128 -2.2731813  -2.27486723
 -2.27523445 -2.26991577]
For degree 3 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(number of intervals)/log(N) : [1.51655163 1.99893052 2.14782837 2.17924027 2.18656585 2.19318829
 2.19544926 2.19454184]
For degree 3 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(mean length of intervals)/log(N) : [-2.40433556 -2.38465994 -2.39771803 -2.39242267 -2.37891891 -2.36835255
 -2.35520857 -2.34278962]
For degree 4 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(number of intervals)/log(N) : [0.94639463 1.66770297 1.89791942 1.95077165 1.97703371 1.99699646
 2.01743008 2.02767278]
For degree 4 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(mean length of intervals)/log(N) : [-2.46287882 -2.39743742 -2.40847747 -2.39718022 -2.38424395 -2.37078936
 -2.35976474 -2.34671655]
For degree 5 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(number of intervals)/log(N) : []
For degree 5 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(mean length of intervals)/log(N) : []
For degree 6 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(number of intervals)/log(N) : []
For degree 6 and N = [9, 18, 30, 39, 48, 60, 78, 99], log(mean length of intervals)/log(N) : []


# the last two degrees are wrong, must redo
"""





if experience_num == 10:
    number_of_points = [30,100,200,400,1000,  1500]
    max_degree = 3
    number_of_intervals = {}
    mean_lengths_of_intervals = {}
    # needed for plotting purposes
    adapted_number_of_points = {}
    for degree in range(max_degree+1):
        number_of_intervals[degree] = []
        mean_lengths_of_intervals[degree] = []
        adapted_number_of_points[degree] = []
    for N in number_of_points:
        print(f"N = {N}")
        length = 2
        X = points_on_parallel_cubes_in_RD(N,5,7,length = length)
        num_int, mean_len =  get_number_and_mean_length_intervals(X, max_degree)
        for degree in range(max_degree+1):
            # to avoid taking the log of zero
            if num_int[degree] > 0:
                adapted_number_of_points[degree].append(N)
                number_of_intervals[degree].append(num_int[degree])
                mean_lengths_of_intervals[degree].append(mean_len[degree])

    for degree in range(max_degree+1):
        epsilon = 2*(np.array(adapted_number_of_points[degree], float)/2)**(-1.0/5.0)
        print(f"epsilon {epsilon}")

        print(f"For degree {degree} and N = {number_of_points}, log(number of intervals)/log(N) ")
        my_string = " & ".join([ f'{result:.2f}' for result in np.log(np.array(number_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))])
        print(my_string)
        print(f"For degree {degree} and N = {number_of_points}, log(mean length of intervals)/log(N) : ")
        my_string = " & ".join([ f'{result:.2f}' for result in np.log(np.array(mean_lengths_of_intervals[degree], float))/np.log(np.array(adapted_number_of_points[degree], float))])
        print(my_string)
        print(f"For degree {degree} and N = {number_of_points}, number of intervals/N : ")
        my_string = " & ".join([ f'{result:.2f}' for result in np.array(number_of_intervals[degree], float)/np.array(adapted_number_of_points[degree], float)])
        print(my_string)
        print(f"For degree {degree} and N = {number_of_points}, mean length of intervals/epsilon^2*200 : ")
        my_string = " & ".join([ f'{result:.2f}' for result in np.array(mean_lengths_of_intervals[degree], float)*epsilon**(-2)*200])
        print(my_string)
        print("----\n")