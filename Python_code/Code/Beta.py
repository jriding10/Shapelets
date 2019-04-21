import numpy as np
import math as m

class Beta:
    def __init__(self, beta1=None, beta2=None, currentMSE=50, currentSSIM=0, maxBeta=1000, nmax=5,
                 smoothedSource=None, minBeta=1.0, stepSize=1.0, checkBeta=None):
        self.beta1 = beta1
        self.beta2 = beta2
        self.currentMSE = currentMSE
        self.currentSSIM = currentSSIM
        self.maxBeta = maxBeta
        self.minBeta = minBeta
        self.stepSize = stepSize
        self.nmax = nmax
        self.smoothedSource = smoothedSource
        self.checkBeta = checkBeta

    def calcAllBetas(self, source, extent, minResolution):
        self.minBeta = minResolution
        maxDim = max(extent[0], extent[1])
        self.maxBeta = 0.9*maxDim*m.pow(self.nmax, -0.52)*minResolution
        number_of_steps = int(self.maxBeta/self.minBeta)
        self.checkBeta = np.zeros((number_of_steps**2, 3))
        self.stepSize = minResolution
        self.beta1 = self.minBeta - self.stepSize
        self.reduceResolution(source)
        k = 0

        for i in range(number_of_steps):
            self.beta1 = self.beta1 + self.stepSize
            self.beta2 = self.minBeta - self.stepSize
            print("Gone through %.0f iterations of %.0f" % (i, number_of_steps))
            for j in range(number_of_steps):
                self.beta2 = self.beta2 + self.stepSize
                model = self.calcModel()
                self.calcSSIM(model)
                self.checkBeta[k, 0] = self.beta1
                self.checkBeta[k, 1] = self.beta2
                self.checkBeta[k, 2] = self.currentSSIM
                k += 1

        np.savetxt("beta_check3.txt", self.checkBeta)
        print("Finished calculating betas")
        major, minor = self.findBestBetas()
        return major, minor

    def findBestBetas(self):
        bestBeta = max(self.checkBeta[:,2])
        k = 0
        while self.checkBeta[k, 2] < bestBeta:
            k += 1
        return self.checkBeta[k, 0], self.checkBeta[k, 1]

    def reduceResolution(self, source, extent):
        max_dim = max(extent[0], extent[1])
        temp_source = np.zeros((int(round(extent[0]/2)), int(round(extent[1]/2))))
        while max_dim > 101:
            l = -2
            rows = int(round(source.shape[0]/2))
            cols = int(round(source.shape[1]/2))
            for i in range(rows-1):
                l += 2
                m = -2
                for j in range(cols-1):
                    m += 2
                    temp_source[i, j] = (source[l, m] + source[l+1, m] + source[l, m+1] +
                                        source[l+1, m+1])/4

            new_row_extent = temp_source.shape[0]
            new_column_extent = temp_source.shape[1]
            max_dim = max(new_row_extent, new_column_extent)
            if max_dim < 100:
                self.smoothedSource = temp_source
            else:
                source = temp_source
                temp_source = np.zeros((int(round(source.shape[0]/2)), int(round(source.shape[1]/2))))
        del temp_source

    def checkBounds(self, beta):
        if beta > self.maxBeta:
            beta = self.maxBeta
        if beta <= 0.0:
            beta = self.maxBeta
        return beta

    def calcModel(self, coords):
        beta1, beta2 = self.checkPositionAngle()
        ucoords = coords[:, 0] / beta1
        vcoords = coords[:, 1] / beta2
        ubasis = shapelet.calcShapeletBasis(ucoords, beta1)
        vbasis = shapelet.calcShapeletBasis(vcoords, beta2)
        length_basis = ubasis.shape[1]
        basis_function = np.zeros((length_basis, 1))
        simple_model = np.zeros((length_basis, 1))

        for n1 in range(self.nmax + 1):
            for n2 in range(self.nmax + 1):
                if (n1 + n2) <= self.nmax:
                    basis_function[:, 0] = np.multiply(ubasis[n1, :], vbasis[n2, :])
                    a_moment = (builder.calcLeastSquares(basis_function))
                    simple_model[:, 0] = simple_model[:, 0] + a_moment*basis_function[:, 0]

        return simple_model

    def checkPositionAngle(self):
        if abs(builder.position_angle) < m.pi / np.sqrt(2):
            beta1 = self.beta1
            beta2 = self.beta2
        else:
            beta1 = self.beta2
            beta2 = self.beta1

        return beta1, beta2

    def calcNMSE(self, model):
        check = self.smoothedSource.shape()
        if len(check) > 0:
            source = self.smoothedSource.flatten()
            resids = source - model
            sum_source = np.sum(source)
            sqerr = np.square(resids)
            self.currentMSE = np.sum(sqerr) / (sum_source * sum_source)
        else:
            self.currentMSE = 1.0

    def calcSSIM(self, model):
        check = self.smoothedSource.shape()
        if len(check) == 0:
            self.currentSSIM = 0.0
            return

        source = self.smoothedSource.flatten()
        source_model = np.multiply(source, model)
        source_mean = np.mean(source)
        source_var = np.var(source)
        model_mean = np.mean(model)
        model_var = np.var(model)
        dynamic_range = max(source) - min(source)
        c1 = 0.01*dynamic_range
        c2 = 0.03*dynamic_range

        cov_xy = np.mean(source_model) - source_mean*model_mean
        luminance = (2*source_mean*model_mean + c1) / (source_mean**2 + model_mean**2 + c1)
        contrast = (2*np.sqrt(source_var*model_var) + c2) / (source_var + model_var + c2)
        structure = (cov_xy + c2/2) / (np.sqrt(source_var*model_var) + c2/2)
        self.currentSSIM = luminance*contrast*structure



