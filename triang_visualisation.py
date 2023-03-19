import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.spatial import Delaunay
from mpl_toolkits import mplot3d


def plot_tri_simple(ax, points, tri):
    for tr in tri.simplices:
        pts = points[tr, :]
        ax.plot3D(pts[[0,1],0], pts[[0,1],1], pts[[0,1],2], color='g', lw='0.1')
        ax.plot3D(pts[[0,2],0], pts[[0,2],1], pts[[0,2],2], color='g', lw='0.1')
        ax.plot3D(pts[[0,3],0], pts[[0,3],1], pts[[0,3],2], color='g', lw='0.1')
        ax.plot3D(pts[[1,2],0], pts[[1,2],1], pts[[1,2],2], color='g', lw='0.1')
        ax.plot3D(pts[[1,3],0], pts[[1,3],1], pts[[1,3],2], color='g', lw='0.1')
        ax.plot3D(pts[[2,3],0], pts[[2,3],1], pts[[2,3],2], color='g', lw='0.1')

    ax.scatter(points[:,0], points[:,1], points[:,2], color='b')




#np.random.seed(0)

def f1(t):
    return np.array([np.ones(len(t))-0.2*(t-0.5)**2,t-0.5 ,np.zeros(len(t))]).T
    #return np.array([t, np.zeros(len(t)),np.zeros(len(t))]).T

def f2(t):
    return np.array([0.1*(t-0.5)**2 ,np.zeros(len(t)), t-0.5]).T
    #return np.array([np.zeros(len(t)),np.zeros(len(t)),t]).T

points1 = f1(np.random.uniform(0,1,8))
points2 = f2(np.random.uniform(0,1,8))


x = 2.0 * np.random.rand(20) - 1.0
y = 2.0 * np.random.rand(20) - 1.0
z = 2.0 * np.random.rand(20) - 1.0
points_0 = np.vstack([x, y, z]).T

print(np.shape(points_0))

points = np.concatenate((points1,points2))

print(np.shape(points))


tri = Delaunay(points)

fig = plt.figure()
ax = plt.axes(projection='3d')
plot_tri_simple(ax, points, tri)
plt.show()

"""

points = np.array([[0, 0,3], [0, 1.1, 2], [1, 0,-1], [1, 1,5]])

tri = Delaunay(points)
print(tri.simplices)

plt.triplot(points[:,0], points[:,1], points[:,2], tri.simplices)

#plt.plot(points[:,0], points[:,1], 'o')

#plt.show()
"""
"""
x=[0,2,1]
y=[0,0,1]

my_triang = tri.Triangulation(x, y)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(my_triang,[0,0,0])
plt.show()

#plot_trisurf(triangulation, ...)
"""