#! /bin/bash/python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

filename = 'beta_check2.txt'

check = np.loadtxt(filename)

fig = plt.figure()
ax = plt.axes(projection='3d')

x = check[:,0]
y = check[:,1]
z = check[:,2]

xmin = min(x)
xmax = max(x)
step = y[1] - y[0]
num_steps = int((xmax-xmin)/step)+1

newx = np.linspace(xmin, xmax, num_steps)
newy = newx

X, Y = np.meshgrid(newx, newy)

nsize = int(np.sqrt(len(z)))
Z = np.zeros((nsize, nsize))
k = 0

for i in range(nsize):
    for j in range(nsize):
        Z[i,j] = z[k]
        k += 1

# beta1, beta2, nmse
ax.plot_surface(X, Y, Z, cmap='cool')
plt.show()

