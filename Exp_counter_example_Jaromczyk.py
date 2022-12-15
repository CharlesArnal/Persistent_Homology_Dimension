import numpy                 as np
import itertools
import matplotlib.pyplot as plt
from points_clouds_generation import points_on_interlaced_2_circles_in_R3, points_on_J_T_parabolas_example



# An example


#  -----------
# Two spheres interlaced, points regularly spaced
# map of the distance between points as a function of the angles


N = 100



u, v = np.mgrid[0:N:1, 0:N:1]

z = np.zeros((N,N))
for i in range(N):
    for j in range(N):
        teta1 = float(i)/N*(np.pi/2) - np.pi/4
        teta2 = float(j)/N*(np.pi/2) - np.pi/4
        x_1 = np.array([1-np.cos(teta1),0,np.sin(teta1)])
        x_2 = np.array([np.cos(teta2),np.sin(teta2),0])
        z[i,j] =  np.linalg.norm(x_1-x_2)
        # z[i,j] = min(np.linalg.norm(x_1-x_2),1.001)


fig = plt.figure( num = "Exp 1" )
ax = plt.axes(projection='3d')
ax.plot_surface(u, v, z, alpha=0.3)
plt.show()


# 3D plot of the dataset


x_1 = []
for i in range(N):
    teta1 = float(i)/N*(np.pi/2) - np.pi/4
    x_1.append(np.array([1-np.cos(teta1),0,np.sin(teta1)]))
x_1 = np.array(x_1)

x_2 = []
for i in range(N):
    teta2 = float(i)/N*(np.pi/2) - np.pi/4
    x_2.append(np.array([np.cos(teta2),np.sin(teta2),0]))
x_2 = np.array(x_2)


fig = plt.figure( num = "Exp 2" )
ax = plt.axes(projection='3d')

ax.scatter3D(x_1[:,0],x_1[:,1],x_1[:,2])
ax.scatter3D(x_2[:,0],x_2[:,1],x_2[:,2])
plt.show()

# --------
# The exact example from Jaromczyk and Toussaint - parabolas and non-regular spacing
# map of the distance between points as a function of position on the parabolas

N = 100

epsilon = 0.05
u, v = np.mgrid[-N:N:1, -N:N:1]

z = np.zeros((2*N,2*N))
for i in range(-N,N,1):
    for j in range(-N,N,1):
        t_1=  i*epsilon
        t_2 = j*epsilon
        x_1 = np.array([1-(t_1**2)/4,0,t_1])
        x_2 = np.array([(t_2**2)/4-1,t_2,0])
        z[i+N,j+N] =  np.linalg.norm(x_1-x_2)
        # z[i,j] = min(np.linalg.norm(x_1-x_2),1.001)

print()
fig = plt.figure( num = "Exp 3" )
ax = plt.axes(projection='3d')
ax.plot_surface(u, v, z, alpha=0.3)
plt.show()

X = points_on_J_T_parabolas_example(N, epsilon)


fig = plt.figure( num = "Exp 4" )
ax = plt.axes(projection='3d')

ax.scatter3D(X[:,0],X[:,1],X[:,2])
plt.show()