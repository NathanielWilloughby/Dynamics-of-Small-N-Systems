import numpy as np
import matplotlib.pyplot as plt

xs = np.loadtxt('xpositions.txt')
ys = np.loadtxt('ypositions.txt')
zs = np.loadtxt('zpositions.txt')
vxs = np.loadtxt('xvelocities.txt')
vys = np.loadtxt('yvelocities.txt')
vzs = np.loadtxt('zvelocities.txt')
origDists = np.loadtxt('origin_distances.txt')
mwithin = np.loadtxt('mass_within.txt')

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
ax.scatter(xs, ys, zs)
ax.set_xlabel('Position along x-axis')
ax.set_ylabel('Position along y-axis')
ax.set_zlabel('Position along z-axis')
plt.show()

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
ax.quiver(xs, ys, zs, vxs, vys, vzs, length=0.1, normalize=True)
plt.show()

rs = np.linspace(np.min(origDists), np.max(origDists), 1000)
plumMassWithin = np.zeros(shape=np.shape(rs))
for i in range(np.shape(rs)[0]):
    plumMassWithin[i] = (rs[i]**3)/((1+(rs[i]**2))**1.5) # Assuming M = 1 and R = 1

plt.scatter(origDists, mwithin)
plt.plot(rs, plumMassWithin)
plt.xlabel('Distance from origin')
plt.ylabel('Mass enclosed')
plt.show()