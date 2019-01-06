#! /bin/bash/python

import numpy as np
import matplotlib.pyplot as plt


filename = 'beta_check1.txt'

check = np.loadtxt(filename)

x = check[:,0]
y = check[:,1]
z = check[:,2]

xmin = min(x)
xmax = max(x)
step = y[1] - y[0]
num_steps = int((xmax-xmin)/step)

num = 91
out = np.zeros((num))
k = 0
for i in range(num):
    out[k] = z[k]
    k += 1

plt.plot(out)
plt.show()
