#! /bin/bash/python

# simple script to create the shapelet basis for pretty plots.

import numpy as np
import math as m
import matplotlib.pyplot as plt

xrange = 1000
yrange = 1000

class Basis:
    def __init__(self, coords=None, moments=None, major=1.0, minor=1.0, ubasis=None, vbasis=None,
                 ucoords=None, vcoords=None, nmax=5):
        self.coords = coords
        self.moments = moments
        self.major = major
        self.minor = minor
        self.ubasis = ubasis
        self.vbasis = vbasis
        self.ucoords = ucoords
        self.vcoords = vcoords
        self.nmax = nmax

    def calcBasis(self):
        self.calcCoords()
        self.ucoords = self.coords[:, 0]/self.major
        self.vcoords = self.coords[:, 1]/self.minor
        self.ubasis = self.calcShapeletBasis(self.ucoords, self.major)
        self.vbasis = self.calcShapeletBasis(self.vcoords, self.minor)
        self.plotBasis()

    def calcCoords(self):
        self.coords = np.zeros((xrange*yrange, 2))
        k = 0

        for i in range(xrange):
            for j in range(yrange):
                self.coords[k, 0] = (i-xrange/2.0)*10.0/xrange
                self.coords[k, 1] = (j-yrange/2.0)*10/yrange
                k += 1

    def calcShapeletBasis(self, zcoords, zbeta):
        total_coords = zcoords.shape[0]
        sq_coords = np.power(zcoords, 2)
        gauss = np.exp(-0.5 * sq_coords)
        norm_const1 = m.sqrt(zbeta * m.sqrt(m.pi))
        shapelet_basis = np.zeros((self.nmax+1, total_coords))

        for n in range(self.nmax + 1):
            norm_const2 = m.sqrt(m.pow(2, n) * m.factorial(n))
            norm_constant = 1.0/(norm_const1 * norm_const2)
            hermite_z = self.calcHermite(n, zcoords)
            shapelet_basis[n, :] = norm_constant * np.multiply(hermite_z, gauss)

        return shapelet_basis

    def calcHermite(self, n, zcoords):
        hermites = np.loadtxt('hermite_coeffs.txt')
        k = 0
        hermite_poly = 0.0

        while k <= n:
            hermite_poly = hermite_poly+hermites[n, k]*np.power(zcoords, (n-k))
            k += 1
        return hermite_poly

    def twoDArray(self, n1, n2):
        new_array = np.zeros((xrange, yrange))
        k = 0

        for i in range(xrange):
            for j in range(yrange):
                new_array[i,j] = self.ubasis[n1, k]*self.vbasis[n2, k]
                k += 1

        return new_array

    def oneDArray(self, n1, n2):
        new_array = np.zeros((xrange*yrange, 3))
        k = 0

        for i in range(xrange):
            for j in range(yrange):
                new_array[k, 0] = self.coords[k, 0]
                new_array[k, 1] = self.coords[k, 1]
                new_array[k, 2] = self.ubasis[n1, k]*self.vbasis[n2, k]
                k += 1

        return new_array

    def plotBasis(self):
        basisTD = self.twoDArray(0, 0)
        pltMax = 0.57
        pltMin = -1*pltMax

        plt.figure(0)
        plt.imshow(basisTD, vmin=pltMin, vmax=pltMax, cmap='hot')
        plt.colorbar()
        plt.axis('off')
        plt.show()

beta = Basis()
beta.calcBasis()