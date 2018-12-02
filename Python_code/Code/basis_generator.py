#! /bin/bash/python

## creates pretty pictures of shapelet basis functions
import numpy as np
import matplotlib.pyplot as plt

class Basis:
    def __init__(self, beta=None, extent=100, coords=None, nmax=3):
        self.beta = beta
        self.extent = extent
        self.coords = coords
        self.nmax = nmax

    def create_coordinates(self):
        self.coords = np.zeros((self.extent,1))
        for i in range(self.extent):
            self.coords[i, 0] = i - self.extent/2

    def calcShapeletBasis(self):
        self.coords = self.coords / self.beta
        sq_coords = np.power(self.coords, 2)
        gauss = np.exp(-0.5 * sq_coords)
        norm_const1 = m.sqrt(self.beta * m.sqrt(m.pi))
        shapelet_basis = np.zeros((self.nmax, self.coords))

        for n in range(self.nmax + 1):
            norm_const2 = m.sqrt(m.pow(2, n) * m.factorial(n))
            norm_constant = 1.0/(norm_const1 * norm_const2)
            hermite_z = self.calcHermite(n)
            shapelet_basis[n, :] = norm_constant * np.multiply(hermite_z, gauss)

        return shapelet_basis

    def calcHermite(self, n):
        hermites = np.loadtxt('hermite_coeffs.txt')
        k = 0
        hermite_poly = 0.0

        while k <= n:
            hermite_poly = hermite_poly+hermites[n, k]*np.power(self.coords, (n-k))
            k += 1
        return hermite_poly

    def makeShapelet(self, basis1, basis2):
        td_basis = np.zeros((self.extent, self.extent))
        for i in range(self.extent):
            for j in range(self.extent):
                td_basis[i, j] = basis1[i]*basis2[j]

        return td_basis

ghp = Basis()
ghp.create_coordinates()
basis_fn = ghp.calcShapeletBasis()
nzerozero = ghp.makeShapelet(basis_fn[0, :], basis_fn[0, :])
nzeroone = ghp.makeShapelet(basis_fn[0, :], basis_fn[1, :])
noneone = ghp.makeShapelet(basis_fn[1, :], basis_fn[1, :])
nzerotwo = ghp.makeShapelet(basis_fn[0, :], basis_fn[2, :])
nonetwo = ghp.makeShapelet(basis_fn[1, :], basis_fn[2, :])
ntwotwo = ghp.makeShapelet(basis_fn[2, :], basis_fn[2, :])

plt.figure(0)
plt.imshow(nzerozero)
plt.colorbar()
plt.figure(1)
plt.imshow(nzerozero)
plt.colorbar()
plt.figure(2)
plt.imshow(nzerozero)
plt.colorbar()
plt.figure(3)
plt.imshow(nzerozero)
plt.colorbar()
plt.figure(4)
plt.imshow(nzerozero)
plt.colorbar()
plt.figure(5)
plt.imshow(nzerozero)
plt.colorbar()
