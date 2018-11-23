#! /bin/bash python
# Shapelets.py

##########################################################################
## A remake of Shapelets.py based upon what I now know. Decomposes an   ##
## image into moments using gauss hermite polynomials as a basis -      ##
## shapelets to astronomers -> f(x,y) = A.g(a, b) where f(x, y) is the  ##
## original function (image), g(a, b) is the new basis (shapelets) and  ##
## A is the transform (or moments).                                     ##
##                                                                      ##
##                                                                      ##
##   Version History                                                    ##
##                                                                      ##
##   v1.1  2018/08/18  JLR     Rewrite of original code                 ##
##                                                                      ##
##                                                                      ##
##########################################################################

import os
import math as m
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

main_window_width = 1000
main_window_height = 1000
myDirectory = 'Users/jriding/Documents/Projects/Guider/TestData/'


class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'Shapelet-ifer'
        self.left = 0
        self.top = 0
        self.width = main_window_width
        self.height = main_window_height
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.main_gui = MainWindow(self)
        self.setCentralWidget(self.noise_app)

        self.show()

class MainWindow(QWidget):

    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout()

        # initialise tab
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tabs.resize(main_window_width, main_window_height)

        # initialise buttons, labels, figures etc
        self.button_file_select = QPushButton('Select Fits File')
        self.label_nbasis = QLabel('Order of basis')
        self.entry_nbasis = QLineEdit(self)
        self.label_name = QLabel('Name of Object')
        self.entry_name = QLineEdit(self)
        self.button_resize = QPushButton('Resize Image')
        self.button_undo = QPushButton('Undo')
        self.button_calculate = QPushButton('Calculate')
        self.canvas = DataDisplay(self.tab1)
        self.canvas2 = DataDisplay(self.tab2)
        self.button_save = QPushButton('Save')
        self.label_ssim = QLabel('SSIM')
        self.data_ssim = QLabel()
        self.label_nmse = QLabel('NMSE')
        self.data_nmse = QLabel()
        self.button_moments = QPushButton('Moment Map')

        #create tabs
        self.tabs.addTab(self.tab1, "Data Selection")
        self.tabs.addTab(self.tab2, "Results")

        self.createTab1()
        self.createTab2()
        self.show()

        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

    def createTab1(self):
        layout = QVBoxLayout(self)
        self.file_select.clicked.connect(lambda: self.findFitsFile())
        layout.addWidget(self.file_select)
        layout.addWidget(self.canvas)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        hbox1 = QHBoxLayout(self)
        self.button_undo.clicked.connect(lambda: self.canvas.displayData(fits_data.source_data))
        self.button_resize.clicked.connect(lambda: coords.selectROI())
        self.button_calculate.clicked.connect(lambda: radio_source.calcShapelet())
        hbox1.addWidget(self.button_undo)
        hbox1.addWidget(self.button_resize)
        hbox1.addWidget(self.button_calculate)
        layout.addWidget(hbox1)
        self.tab1.setLayout(layout)


    @pyqtSlot()
    def findFitsFile(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filetypes = 'Fits files (*.fits);;All files (*)'
        fits_data.filename, _ = QFileDialog.getOpenFileName(self, 'Data', myDirectory, filetypes, options=options)
        fits_data.getFitsFile()
        canvas.displayData(fits_data.source_data)

    def on_press(self, event):
        coords.row0 = int(event.xdata)
        coords.column0 = int(event.ydata)

    def on_release(self, event):
        coords.row1 = int(event.xdata)
        coords.column1 = int(event.ydata)
        coords.row_extent = abs(coords.row1 - coords.row0)
        coords.column_extent = abs(coords.column1 - coords.column0)
        self.canvas.displayRectangle()

    def getUserInputs(self):
        fits_data.source_name = self.entry_name.text()
        radio_source.max_basis_index = self.entry_nbasis.text()

    def updateDisplay(self):


class DataDisplay(FigureCanvas):

    def __init__(self, parent=None):
        fig = Figure(figsize=im_size, dpi=10)
        self.axes = fig.add_subplot(111)
        self.displayDefault()

        FigureCanvas.__init__(self, fig)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.setParent(parent)

    def displayDefault(self):
        dataToDisplay = np.zeros((int(main_window_width), int(main_window_height/2)))
        self.axes.imshow(dataToDisplay)

    def displayData(self, dataToDisplay):
        self.axes.cla()
        if len(dataToDisplay) != 0:
            self.axes.imshow(dataToDisplay)
            self.draw()
        else:
            self.displayDefault()

    def displayRectangle(self):
        delta_row = abs(coords.row0 - coords.row1)
        delta_column = abs(coords.roi_column0 - coords.roi_column1)
        rect = plt.Rectangle((coords.row0, coords.column0), delta_row, delta_column,
                             linewidth=2, edgecolor='y', facecolor='none')
        self.axes.add_patch(rect)
        self.draw()


class FitsData:
    def __init__(self, filename=None, RA=None, Dec=None, xAngScale=None, yAngScale=None, rowDecpxl=None,
                 colRApxl=None, xExtent=None, yExtent=None, source_data=None, source_name='temp'):
        self.filename = filename
        self.RA = RA
        self.Dec = Dec
        self.xAngScale = xAngScale
        self.yAngScale = yAngScale
        self.rowDecpxl = rowDecpxl
        self.colRApxl = colRApxl
        self.xExtent = xExtent
        self.yExtent = yExtent
        self.source_data = source_data
        self.source_name = source_name

    def getFitsFile(self):
        fitsFile = pyfits.open(self.filename)
        fitsHeader = fitsFile[0].header
        fitsData = fitsFile[0].data

        # RA == y axis/column indices, Dec == x axis/row indices
        self.RA = np.radians(fitsHeader['CRVAL1'])
        self.Dec = np.radians(fitsHeader['CRVAL2'])
        self.xAngScale = np.radians(fitsHeader['CDELT2'])
        self.yAngScale = np.radians(fitsHeader['CDELT2'])
        self.colDecpxl = fitsHeader['CRPIX2']
        self.colRApxl = fitsHeader['CRPIX1']
        self.xExtent = fitsHeader['NAXIS2']
        self.yExtent = fitsHeader['NAXIS1']
        fitsFile.close()

        self.checkSourceData(fitsData)
        self.copyData()

    def checkSourceData(self, fitsdata):
        size_of_data = fitsdata.shape
        dim_of_data = len(size_of_data)

        if dim_of_data == 2:
            self.source_data = fitsdata
        elif dim_of_data == 4:
            self.source_data = fitsdata[0, 0, :, :]
        else:
            print('Fits data is of unknown dimensions: %.0f' % dim_of_data)
            sys.exit(app.exit())

    def copyData(self):
        coords.colRA_pxl = self.colRApxl
        coords.rowDec_pxl = self.rowDecpxl
        coords.xextent = self.xExtent
        coords.yextent = self.yExtent
        builder.colRA = self.RA
        builder.rowDec = self.Dec
        radio_source.extended_source = self.source_data

class TheSource:
    def __init__(self, extended_source=None, max_basis_index=5, normalised_mse=0.0, avg_ssim=1.0):
        self.extended_source = extended_source
        self.max_basis_index = max_basis_index
        self.normalised_mse = normalised_mse
        self.avg_ssim = avg_ssim

    def calcShapelet(self):
        row_centre_pxl = m.floor(coords.row_extent/2)
        col_centre_pxl = m.floor(coords.column_extent/2)
        coords.changeRADecRef(row_centre_pxl, col_centre_pxl)
        coords.calcCoordinateList()
        self.calcCovarianceFit()
        self.calcTransformCoordinates()
    #        beta.calcMajorMinor()
        builder.calcMoments()
        builder.saveShapelet()
        self.setupModelSource()
        builder.calcShapeletModel()
        builder.calcResidual()
        self.calcSSIM()
        self.calcNMSE()
        ex.main_gui.updateDisplay()

    def calcCovarianceFit(self):
        sum_over_source = np.sum(np.sum(self.extended_source))

        weighted_centre = self.calcWeightedCentre(sum_over_source)
        self.pxl_coords[:, 0] -= weighted_centre[0]
        self.pxl_coords[:, 1] -= weighted_centre[1]
        builder.colRA += weighted_centre[1]
        builder.rowDec -= weighted_centre[0]

        covar_matrix = self.calcFluxMatrix(sum_over_source)
        (eigenval, eigenvect) = self.calcEigenvalueDecomp(covar_matrix)
        self.calcPositionAngle(eigenvect)

    def calcWeightedCentre(self, factor):
        k = 0
        row_centre = 0.0
        col_centre = 0.0
        for i in range(coords.row_extent):
            for j in range(coords.column_extent):
                row_centre += self.extended_source[i, j]*coords.pxl_coords[k, 0]
                col_centre += self.extended_source[i, j]*coords.pxl_coords[k, 1]
                k += 1

        row_centre = row_centre/factor
        col_centre = col_centre/factor
        weighted_centre = ([row_centre, col_centre])
        return weighted_centre

    def calcFluxMatrix(self, factor):
        # Matrix terms
        row_x_row = 0.
        col_x_col = 0.
        row_x_col = 0.

        k = 0
        for i in range(coords.row_extent):
            for j in range(coords.column_extent):
                row_x_row = row_x_row + self.extended_source[i, j]*coords.pxl_coords[k, 0]*coords.pxl_coords[k, 0]
                col_x_col = col_x_col + self.extended_source[i, j]*coords.pxl_coords[k, 1]*coords.pxl_coords[k, 1]
                row_x_col = row_x_col + self.extended_source[i, j]*coords.pxl_coords[k, 0]*coords.pxl_coords[k, 1]
                k += 1

        row_x_row = row_x_row/factor
        col_x_col = col_x_col/factor
        row_x_col = row_x_col/factor
        covar_matrix = ([row_x_row, row_x_col, col_x_col])
        return covar_matrix

    def calcEigenvalueDecomp(self, covar_matrix):
        eigenmatrix = np.array([[covar_matrix[0], covar_matrix[1]],
                                [covar_matrix[1], covar_matrix[2]]])

        (eigenval, eigenvect) = np.linalg.eig(eigenmatrix)
        return (eigenval, eigenvect)

    def calcPositionAngle(self, eigenvect):
        defaultMajor, defaultMinor = self.calcDefaultValues()
        builder.position_angle = -1*m.atan2(eigenvect[1, 1], eigenvect[0, 1])
        builder.major = defaultMajor
        builder.minor = defaultMinor

    def calcDefaultValues(self):
        row_axis = 0.9*coords.row_extent*m.pow(self.max_basis_index, -0.52)
        col_axis = 0.9*coords.column_extent*m.pow(self.max_basis_index, -0.52)

        if row_axis > col_axis:
            defaultMajor = row_axis*abs(fits_data.xAngScale)
            defaultMinor = col_axis*abs(fits_data.yAngScale)
            builder.position_angle = m.pi/2
        else:
            defaultMajor = col_axis*abs(fits_data.yAngScale)
            defaultMinor = row_axis*abs(fits_data.xAngScale)
            builder.position_angle = 0.0

        return defaultMajor, defaultMinor

    def calcTransformCoordinates(self):
        PA = builder.position_angle
        transform_matrix = [[m.cos(PA), -m.sin(PA)], [m.sin(PA), m.cos(PA)]]
        coord_list = np.dot(transform_matrix, np.transpose(coords.pxl_coords))
        coords.pxl_coords = np.transpose(coord_list)

    def setupModelSource(self):
        model_source.__dict__ = radio_source.__dict__.copy()
        residual.__dict__ = radio_source.__dict__.copy()

    def calcSSIM(self):
        cov_xy = 1

        source_mean = np.mean(source.extended_source.flatten())
        source_var = np.var(source.extended_source.flatten())
        model_mean = np.mean(self.extended_source.flatten())
        model_var = np.var(self.extended_source.flatten())

        numerator_ssim = 4 * source_mean * model_mean * cov_xy
        denominator_ssim = (source_mean ** 2 + model_mean ** 2) * (source_var ** 2 + model_var ** 2)
        self.avg_ssim = numerator_ssim / denominator_ssim

        ssim_txt = '%3.4f' % self.avg_ssim
        app.ssim_value.configure(text=ssim_txt)

    def calcNMSE(self):
        sum_source = np.sum(radio_source.extended_source.flatten())
        sqerr = np.square(residual.extended_source.flatten())

        self.normalised_mse = np.sum(sqerr) / (sum_source * sum_source)
        nmse_txt = '%3.4f' % self.normalised_mse
        app.nmse_value.configure(text=nmse_txt)


class Coordinates:
    def __init__(self, colRA_pxl=None,  rowDec_pxl=None, row_extent=None, column_extent=None, pxl_coords=None,
                 roi_row0=None, roi_column0=None, roi_row1=None, roi_column1=None):
        self.colRA_pxl = colRA_pxl
        self.rowDec_pxl = rowDec_pxl
        self.row_extent = row_extent
        self.column_extent = column_extent
        self.pxl_coords = pxl_coords
        self.roi_row0 = roi_row0
        self.roi_column0 = roi_column0
        self.roi_row1 = roi_row1
        self.roi_column1 = roi_column1

    def changeRADecRef(self, new_row, new_column):
        new_row = int(new_row)
        new_column = int(new_column)
        deltaRowPixel = abs(self.rowDec_pxl - new_row)
        deltaColumnPixel = abs(self.colRA_pxl - new_column)
        newRA = builder.rowRA - fits_data.yAngScale*deltaColumnPixel
        newDec = builder.rowDec - fits_data.xAngScale*deltaRowPixel

        self.colRA_pxl = new_column
        self.rowDec_pxl = new_row
        builder.rowRA = newRA
        builder.colDec = newDec

    def selectROI(self):
        self.changeRADecRef(coords.roi_row0, coords.roi_column0)
        new_source = np.zeros((self.row_extent, self.column_extent))

        for i in range(self.row_extent):
            for j in range(self.column_extent):
                rows = int(i + self.roi_row0)
                cols = int(j + self.roi_column0)
                new_source[i, j] = radio_source.extended_source[rows, cols]

        self.colRA_pxl = 0
        self.rowDec_pxl = 0
        radio_source.extended_source = new_source
        ex.main_gui.getUserInputs()
        ex.main_gui.canvas.displayFitsData(new_source)

    def calcCoordinateList(self):
        total_coordinates = self.row_extent * self.column_extent
        coordinates = np.zeros((total_coordinates, 2))
        row_start = int(-1 * self.row_extent/2)
        column_start = int(-1 * self.column_extent/2)

        k = 0
        for i in range(self.row_extent):
            for j in range(self.column_extent):
                coordinates[k, 0] = (row_start + i) * abs(fits_data.xAngScale)
                coordinates[k, 1] = (column_start + j) * abs(fits_data.yAngScale)
                k += 1
        self.pxl_coords = coordinates

    def expandArray(self, flat_array):
        num_rows = self.row_extent
        num_cols = self.column_extent
        new_array = np.zeros((num_rows, num_cols))
        k = 0

        for i in range(num_rows+1):
            for j in range(num_cols+1):
                new_array[i,j] = flat_array[k]
                k += 1

        return new_array


class Shapelets:
    def __init__(self, colRA=None, rowDec=None, moments=None, major=3.0, minor=1.0,
                position_angle=0.0, ubasis=None, vbasis=None, ucoords=None, vcoords=None):
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

    def calcMoments(self):
        if abs(self.position_angle) < m.pi/np.sqrt(2):
            beta_u = self.major
            beta_v = self.minor
        else:
            beta_u = self.minor
            beta_v = self.major

        self.ucoords = coords.pxl_coords[:, 0]/beta_u
        self.vcoords = coords.pxl_coords[:, 1]/beta_v

        n_max = radio_source.max_basis_index
        self.ubasis = self.calcShapeletBasis(self.ucoords, beta_u)
        self.vbasis = self.calcShapeletBasis(self.vcoords, beta_v)

        total_moments = int(0.5*(n_max + 1)*(n_max + 2))
        self.moments = np.zeros((total_moments, 3))

        k = 0
        for n1 in range(n_max + 1):
            for n2 in range(n_max + 1):
                if (n1 + n2) <= n_max:
                    basis_function = np.multiply(self.ubasis[n1, :], self.vbasis[n2, :])
                    a_moment = (self.calcLeastSquares(basis_function))
                    self.moments[k, 0] = n1
                    self.moments[k, 1] = n2
                    self.moments[k, 2] = a_moment
                    k += 1

        self.moments = self.moments[(np.argsort(abs(self.moments[:,2])))]
        self.moments = np.flipud(self.moments)

    def calcShapeletBasis(self, zcoords, zbeta):
        total_coords = zcoords.shape[0]
        sq_coords = np.power(zcoords, 2)
        gauss = np.exp(-0.5 * sq_coords)
        norm_const1 = m.sqrt(zbeta * m.sqrt(m.pi))
        shapelet_basis = np.zeros((radio_source.max_basis_index+1, total_coords))

        for n in range(radio_source.max_basis_index + 1):
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

    def calcLeastSquares(self, basis_function):
        flat_source = radio_source.extended_source.flatten()
        transpose_basis = np.transpose(basis_function)
        denominator = np.dot(transpose_basis, basis_function)
        numerator = np.dot(transpose_basis, flat_source)
        a_moment = numerator/denominator
        return a_moment

    def minimiseMoments(self, value):
        new_moments = np.zeros((value, 3))
        for i in range(value+1):
            new_moments[i, :] = self.moments[i, :]
        self.moments = new_moments

    def saveShapelet(self):
        file_name = '../Models/' + fits_data.source_name + '_moments.txt'
        np.savetxt(file_name, self.moments)

        beta1 = np.degrees(self.major)
        beta2 = np.degrees(self.minor)
        shapelet_parameters = ([np.degrees(self.colRA), np.degrees(self.rowDec), beta1,
                                beta2, np.degrees(self.position_angle)])

        file_name = '../Models/' + fits_data.source_name + '_parameters.txt'
        np.savetxt(file_name, shapelet_parameters)

    def calcShapeletModel(self):
        flat_model = np.zeros((coords.row_extent*coords.column_extent, 1))
        total_moments = self.moments.shape[0]

        for k in range(total_moments):
            n1 = int(self.moments[k, 0])
            n2 = int(self.moments[k, 1])
            a_moment = self.moments[k, 2]
            basis_function = np.multiply(self.ubasis[n1, :], self.vbasis[n2, :])
            flat_model = np.add(flat_model, a_moment*basis_function)

        model_source.extended_source = flat_model
        model_source.extended_source = coords.expandArray(flat_model)

    def calcResidual(self):
        residual.extended_source = radio_source.extended_source - model_source.extended_source


class Beta:
    def __init__(self, previousMajor=3.0, previousMinor=1.0,
                 previousMSE=100, currentMSE=50, maxBeta=1000,
                 minBeta=1.0, stepSize=1.0, bestMajor=0, bestMinor=0):
        self.previousMajor = previousMajor
        self.previousMinor = previousMinor
        self.previousMSE = previousMSE
        self.currentMSE = currentMSE
        self.maxBeta = maxBeta
        self.minBeta = minBeta
        self.stepSize = stepSize
        self.bestMajor = bestMajor
        self.bestMinor = bestMinor

    def calcMajorMinor(self):
        maxDim = max(coords.row_extent, coords.column_extent)
        minResolution = min(fits_data.xAngScale, fits_data.yAngScale)
        self.minBeta = minResolution
        self.maxBeta = 0.9*maxDim*m.pow(radio_source.max_basis_index, -0.52)*self.minBeta
        self.previousMajor = builder.major
        self.previousMinor = builder.minor
        solution = 0

        while solution != 2:
            self.currentMSE = 5
            self.previousMSE = 10
            if self.bestMajor == 0:
                self.calcMajor()

            self.checkMajorSolution()
            self.currentMSE = 5
            self.previousMSE = 10
            if self.bestMinor == 0:
                self.calcMinor()

            self.checkMinorSolution()
            solution = self.bestMajor + self.bestMinor


    def calcMajor(self):
        builder.major = self.previousMajor + 5 * self.stepSize
        if builder.major > self.maxBeta:
            builder.major = self.maxBeta
        while beta.currentMSE < beta.previousMSE:
            self.calcMajorShapelet()
            self.currentMSE = radio_source.normalised_mse


    def calcMajorShapelet(self):
        self.previousMSE = self.currentMSE
        builder.major = builder.major - self.stepSize
        if builder.major < 0.0:
            builder.major = self.maxBeta
        builder.calcMoments()
        builder.calcShapeletModel()
        builder.calcResidual()
        radio_source.calcNMSE()


    def checkMajorSolution(self):
        currentMajor = builder.major + self.stepSize
        if round((currentMajor * 10e4)) == round(self.previousMajor * 10e4):
            builder.major = self.previousMajor
            self.bestMajor = 1
        else:
            self.previousMajor = builder.major + self.stepSize
            builder.major = self.previousMajor
            self.bestMinor = 0
            self.bestMajor = 0


    def calcMinor(self):
        builder.minor = self.previousMinor + 5 * self.stepSize
        if builder.minor > self.maxBeta:
            builder.minor = self.maxBeta
        while self.currentMSE < self.previousMSE:
            self.calcMinorShapelet()
            self.currentMSE = radio_source.normalised_mse


    def calcMinorShapelet(self):
        self.previousMSE = self.currentMSE
        builder.minor = builder.minor - self.stepSize
        if builder.minor < 0.0:
            builder.minor = self.maxBeta
        builder.calcMoments()
        builder.calcShapeletModel()
        builder.calcResidual()
        radio_source.calcNMSE()


    def checkMinorSolution(self):
        currentMinor = builder.minor + self.stepSize
        if round((currentMinor * 10e4)) == round(self.previousMajor * 10e4):
            builder.minor = self.previousMinor
            beta.bestMinor = 1
        else:
            self.previousMinor = builder.minor + self.stepSize
            builder.minor = self.previousMinor
            self.bestMinor = 0
            self.bestMajor = 0


if __name__ == '__main__':
    app = QApplication(sys.argv)
    fits_data = FitsData()
    radio_source = TheSource()
    model_source = TheSource()
    residual = TheSource()
    coords = Coordinates()
    builder = Shapelets()
    beta = Beta()
    ex = App()
    sys.exit(app.exec_())
