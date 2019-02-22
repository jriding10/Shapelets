import numpy as np
import math as m

class Shapelets:
    def __init__(self, colRA=None, rowDec=None, moments=None, major=3.0, minor=1.0, min_moments=None,
                position_angle=0.0, ubasis=None, vbasis=None, ucoords=1.0, vcoords=1.0):
        self.colRA = colRA
        self.rowDec = rowDec
        self.moments = moments
        self.major = major
        self.minor = minor
        self.position_angle = position_angle
        self.ubasis = ubasis
        self.vbasis = vbasis
        self.ucoords = ucoords
        self.vcoords = vcoords
        self.min_moments = min_moments

    def calcMoments(self, source, coords, n_max):
        if abs(self.position_angle) < m.pi/np.sqrt(2):
            beta_u = self.major
            beta_v = self.minor
        else:
            beta_u = self.minor
            beta_v = self.major

        self.ucoords = coords[:, 0]/beta_u
        self.vcoords = coords[:, 1]/beta_v
        self.ubasis = self.calcShapeletBasis(self.ucoords, beta_u, n_max)
        self.vbasis = self.calcShapeletBasis(self.vcoords, beta_v, n_max)
        length_basis = self.ubasis.shape[1]
        basis_function = np.zeros((length_basis, 1))

        total_moments = int(0.5*(n_max + 1)*(n_max + 2))
        self.moments = np.zeros((total_moments, 3))

        k = 0
        for n1 in range(n_max + 1):
            for n2 in range(n_max + 1):
                if (n1 + n2) <= n_max:
                    basis_function[:, 0] = np.multiply(self.ubasis[n1, :], self.vbasis[n2, :])
                    a_moment = (self.calcLeastSquares(basis_function, source))
                    self.moments[k, 0] = n1
                    self.moments[k, 1] = n2
                    self.moments[k, 2] = a_moment
                    k += 1
        self.moments = self.moments[(np.argsort(abs(self.moments[:,2])))]
        self.moments = np.flipud(self.moments)

    def calcShapeletBasis(self, zcoords, zbeta, nmax):
        total_coords = zcoords.shape[0]
        sq_coords = np.power(zcoords, 2)
        gauss = np.exp(-0.5 * sq_coords)
        norm_const1 = m.sqrt(zbeta * m.sqrt(m.pi))
        shapelet_basis = np.zeros((nmax+1, total_coords))

        for n in range(nmax + 1):
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

    def calcLeastSquares(self, basis_function, source):
        num_coords = basis_function.shape[0]
        basis = np.zeros((num_coords,1))
        if len(basis_function.shape) == 1:
            basis[:,0] = basis_function[:]
        else:
            basis[:,0] = basis_function[:,0]

        flat_source = source.flatten()
        transpose_basis = np.transpose(basis_function)
        denominator = np.dot(transpose_basis, basis_function)
        numerator = np.dot(transpose_basis, flat_source)
        a_moment = numerator/denominator
        return a_moment

    def minimiseMoments(self, numToInc):
        min_moments = np.zeros((numToInc, 3))
        for i in range(numToInc):
            min_moments[i, :] = self.moments[i, :]
        return min_moments

    def calcShapeletModel(self, coords, moments):
        num_coords = coords.shape[0]
        flat_model = np.zeros((num_coords, 1))
        basis_function = np.zeros((num_coords, 1))
        total_moments = moments.shape[0]

        for k in range(total_moments):
            n1 = int(moments[k, 0])
            n2 = int(moments[k, 1])
            a_moment = moments[k, 2]
            basis_function[:, 0] = np.multiply(self.ubasis[n1, :], self.vbasis[n2, :])
            flat_model[:, 0] = np.add(flat_model[:, 0], a_moment*basis_function[:, 0])

        return flat_model

    def calcResidual(self, source, model):
        irange = model.shape[0]
        resids = np.zeros((irange,1))
        for i in range(irange):
            resids[i] = source[i] - model[i]
        return resids

    def makeMomentArray(self):
        max_moments = int(max(max(self.moments[:,0]), max(self.moments[:,1])) + 1)
        moment_array = np.zeros((max_moments, max_moments))

        for i in range(self.moments.shape[0]):
            x = max_moments - int(self.moments[i,0]) - 1
            y = int(self.moments[i,1])
            moment_array[x,y] = self.moments[i,2]
        return moment_array