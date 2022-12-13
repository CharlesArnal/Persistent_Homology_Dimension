# Persistent_Homology_Dimension
A few experiments on persistent homology dimension

points_clouds_generation.py contains a few functions that can generate points clouds in various geometric configurations as numpy arrays

compute_total_persistence.py defines a function that computes the total persistence, for the Wasserstein distance of a given order, of a points cloud (for the alpha-filtration)
It can also plot the persistence diagrams of the points cloud for each homology degree

plot_crit_points.py plots in 3D the centers of the (smallest) circumsphere of the critical simplices of a given dimension of the alpha filtration of a points cloud in R³
You can choose whether to plot positive or negative simplices
