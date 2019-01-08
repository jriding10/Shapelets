import numpy as np
import math as m


class TheSource:
    def __init__(self, extended_source=None, max_basis_index=25, normalised_mse=0.0, avg_ssim=1.0, sum_source=0.0):
        self.extended_source = extended_source
        self.max_basis_index = max_basis_index
        self.normalised_mse = normalised_mse
        self.avg_ssim = avg_ssim
        self.sum_source = sum_source

    def sumOfSource(self):
        self.sum_source = np.sum(np.sum(self.extended_source))

    def calcWeightedCentre(self, xextent, yextent, coords):
        k = 0
        row_centre = 0.0
        col_centre = 0.0
        for i in range(xextent):
            for j in range(yextent):
                row_centre += self.extended_source[i, j]*coords[k, 0]
                col_centre += self.extended_source[i, j]*coords[k, 1]
                k += 1

        row_centre = row_centre/self.sum_source
        col_centre = col_centre/self.sum_source
        weighted_centre = ([row_centre, col_centre])
        return weighted_centre

    def calcFluxMatrix(self, xextent, yextent, coords):
        # Matrix terms
        row_x_row = 0.
        col_x_col = 0.
        row_x_col = 0.

        k = 0
        for i in range(xextent):
            for j in range(yextent):
                row_x_row = row_x_row + self.extended_source[i, j]*coords[k, 0]*coords[k, 0]
                col_x_col = col_x_col + self.extended_source[i, j]*coords[k, 1]*coords[k, 1]
                row_x_col = row_x_col + self.extended_source[i, j]*coords[k, 0]*coords[k, 1]
                k += 1

        row_x_row = row_x_row/self.sum_source
        col_x_col = col_x_col/self.sum_source
        row_x_col = row_x_col/self.sum_source
        covar_matrix = ([row_x_row, row_x_col, col_x_col])
        return covar_matrix

    def calcEigenvalueDecomp(self, covar_matrix):
        eigenmatrix = np.array([[covar_matrix[0], covar_matrix[1]],
                                [covar_matrix[1], covar_matrix[2]]])

        (eigenval, eigenvect) = np.linalg.eig(eigenmatrix)
        return (eigenval, eigenvect)

    def calcPositionAngle(self, eigenvect):
        defaultMajor, defaultMinor, pa = self.calcDefaultValues()
        position_angle = -1*m.atan2(eigenvect[1, 1], eigenvect[0, 1])-m.pi
        return defaultMajor, defaultMinor, position_angle

    def calcDefaultValues(self, angScale, extent):
        row_axis = 0.9*extent[0]*m.pow(self.max_basis_index, -0.52)
        col_axis = 0.9*extent[1]*m.pow(self.max_basis_index, -0.52)

        if row_axis > col_axis:
            defaultMajor = row_axis*abs(angScale[0])
            defaultMinor = col_axis*abs(angScale[1])
            position_angle = m.pi/2
        else:
            defaultMajor = col_axis*abs(angScale[1])
            defaultMinor = row_axis*abs(angScale[0])
            position_angle = 0.0

        return defaultMajor/2.0, defaultMinor/2.0, position_angle

    def calcTransformCoordinates(self, pa, coords):
        transform_matrix = [[m.cos(pa), -m.sin(pa)], [m.sin(pa), m.cos(pa)]]
        coord_list = np.dot(transform_matrix, np.transpose(coords))
        coords = np.transpose(coord_list)
        return coords

    def calcSSIM(self, source, model):
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
        self.avg_ssim = luminance*contrast*structure
        return self.avg_ssim

    def calcNMSE(self, source, resids):
        sum_source = np.sum(source)
        sqerr = np.square(resids)
        self.normalised_mse = np.sum(sqerr) / (sum_source * sum_source)
        return self.normalised_mse
