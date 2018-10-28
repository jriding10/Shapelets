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
import matplotlib
import astropy.io.fits as pyfits
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk
import matplotlib.pyplot as plt


class ShapeletGUI(ttk.Frame):
    def __init__(self, master):
        ttk.Frame.__init__(self, master)
        self.grid()

        self.file_label = Label(self, fg='black', justify='center', font='Helvetica 12')
        self.file_label.configure(relief='groove', bg='white')
        self.get_file = Button(self, text='Select', command=self.getFilename, width=12)
        self.basis_label = Label(self, text='Number of basis:', fg='black', font='Helvetica 12')
        self.basis_entry = Entry(self, width=3, fg='red')
        self.basis_entry.bind('<Return>', self.getNumberBasis)
        self.source_label = Label(self, text='Source name:', fg='black', font='Helvetica 12')
        self.source_entry = Entry(self, width=10, fg='red')
        self.source_entry.bind('<Return>', self.getSourceName)
        self.coeff_label = Button(self, text='Display Moments', command=self.momentPlot, width=12)
        self.quitter = Button(self, text='Quit', command=self.quitter, width=12)
        self.load_file = Button(self, text='Load', command=dataset.getFitsFile, width=12)
        self.resize_im = Button(self, text='Resize', command=im_dims.resizeImageOfSource, width=12)
        self.undo_resize = Button(self, text='Undo', command=self.undoResize, width=12)
        self.process_data = Button(self, text='Go', command=self.createShapelet)
        self.calc_done = Label(self, justify='center', fg='black', font='Helvetica 12')
        self.ssim_label = Label(self, text='SSIM', justify='center', font='Helvetica 12')
        self.ssim_value = Label(self, justify='left', font='Helvetica 12', bg='white', width=12)
        self.nmse_label = Label(self, text='NMSE', justify='center', font='Helvetica 12')
        self.nmse_value = Label(self, justify='left', font='Helvetica 12', bg='white', width=12)

        self.fig = plt.figure(1)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(column=1, row=3, sticky=(N, W, E, S))

        self.fig.canvas.mpl_connect('button_press_event', self.xy)
        self.fig.canvas.mpl_connect('button_release_event', self.selectROI)

        self.get_file.grid(column=0, row=0)
        self.file_label.grid(column=1, row=0, sticky=(E, W))
        self.load_file.grid(column=2, row=0)
        self.basis_label.grid(column=0, row=1, sticky='E')
        self.basis_entry.grid(column=1, row=1, sticky='W')
        self.source_label.grid(column=0, row=2, sticky='E')
        self.source_entry.grid(column=1, row=2, sticky='W')
        self.coeff_label.grid(column=2, row=5)
        self.quitter.grid(column=2, row=6)
        self.resize_im.grid(column=2, row=1)
        self.undo_resize.grid(column=2, row=2, sticky='W')
        self.process_data.grid(column=1, row=4, sticky=(E, W))
        self.calc_done.grid(column=0, row=4)
        self.ssim_label.grid(column=0, row=5, sticky='E')
        self.ssim_value.grid(column=1, row=5, sticky='W')
        self.nmse_label.grid(column=0, row=6, sticky='E')
        self.nmse_value.grid(column=1, row=6, sticky='W')

    def quitter(self):
        if messagebox.askyesno('Verify', 'Do you really want to quit?'):
            plt.close('all')
            root.destroy()

    def getFilename(self):
        dir_name = os.getcwd()
        window_title = 'Select fits file'
        file_types = (('fits files', '*.fits'), ('all files', '*.*'))
        dataset.filename = filedialog.askopenfilename(initialdir = dir_name,
                                          title = window_title,
                                          filetypes = file_types)
        self.file_label.configure(text=dataset.filename)
        self.update_idletasks()

    def displayFitsData(self):
        self.ax.clear()
        plt.imshow(source.extended_source, cmap='hot')
        self.fig.canvas.draw()

    def undoResize(self):
        self.copyData()
        self.displayFitsData()

    def getNumberBasis(self, event):
        n = self.basis_entry.get()
        source.max_basis_index = int(n)
        self.basis_entry.config(fg='black')
        self.update_idletasks()

    def getSourceName(self, event):
        dataset.source_name = self.source_entry.get()
        self.source_entry.config(fg='black')
        self.update_idletasks()

    def xy(self, event):
        lastx, lasty = event.xdata, event.ydata
        region.x0 = int(lastx)
        region.y0 = int(lasty)

    def selectROI(self, event):
        lastx = region.x0
        lasty = region.y0
        x, y = event.xdata, event.ydata
        region.x1 = int(x)
        region.y1 = int(y)
        deltax = x - lastx
        deltay = y - lasty
        rectangle = plt.Rectangle((lastx,lasty), deltax, deltay, linewidth=2,
                                  edgecolor='y', facecolor='none')
        self.ax.add_patch(rectangle)
        self.fig.canvas.draw()

    def createShapelet(self):
        self.calc_done.configure(text='Thinking', bg='red')
        self.update_idletasks()
        im_dims.calcCoordinateList()
        source.calcCovarianceFit()
        shapelet.calcMoments()
        shapelet.saveShapelet()

        self.setupModelSource()
        shapelet.calcShapeletModel()
        shapelet.calcResidual()
        model.calcNMSE()
        model.calcSSIM()
        self.calc_done.configure(text='Done', bg='green')
        self.update_idletasks()
        self.displayResults()

    def setupModelSource(self):
        model.__dict__ = source.__dict__.copy()
        model_dims.__dict__= im_dims.__dict__.copy()
        residual.__dict__ = source.__dict__.copy()
        residual_dims.__dict__= im_dims.__dict__.copy()

    def displayResults(self):
        results = Toplevel()
        results.wm_title('Shapelet Model')

        scale_max, scale_min = source.getMinMax()

        data = plt.figure(2)
        canvas = FigureCanvasTkAgg(data, master=results)
        canvas.get_tk_widget().grid(column=0, row=0, sticky=(N,W,E,S))

        plt.subplot(221)
        plt.imshow(source.extended_source, vmin=scale_min, vmax=scale_max, cmap='hot')

        plt.subplot(222)
        plt.imshow(model.extended_source, vmin=scale_min, vmax=scale_max, cmap='hot')

        plt.subplot(223)
        plt.imshow(residual.extended_source, vmin=scale_min, vmax=scale_max, cmap='hot')

        data.canvas.draw()

    def momentPlot(self):
        coeffs = Toplevel()
        coeffs.wm_title('Coefficient Power Distribution')

        data = plt.figure(3)
        canvas = FigureCanvasTkAgg(data, master=coeffs)
        canvas.get_tk_widget().grid(column=0, row=0, sticky=(N,W,E,S))

        moment_plot = self.makeMomentArray()
        plt.subplot(111)
        plt.imshow(moment_plot, cmap='hot')

    def makeMomentArray(self):
        totalMoments = shapelet.moments.shape[0]
        moment_array = np.zeros((source.max_basis_index+1), (source.max_basis_index+1))

        for i in range(totalMoments):
            j = int(shapelet.moments[i, 0])
            k = int(shapelet.moments[i, 1])
            moment_array[j, k] = shapelet.moments[i, 2]

        return moment_array


class Fits_data:
    def __init__(self, filename=None, RA=0.0, Dec=0.0, xAngScale=1.0, yAngScale=1.0, xDecpxl=0,
                 yRApxl=0, xExtent=1, yExtent=1, source_data=[], source_name='temp'):
        self.filename = filename
        self.RA = RA
        self.Dec = Dec
        self.xAngScale = xAngScale
        self.yAngScale = yAngScale
        self.xDecpxl = xDecpxl
        self.yRApxl = yRApxl
        self.xExtent = xExtent
        self.yExtent = yExtent
        self.source_data = source_data
        self.source_name = source_name

    def getFitsFile(self):
        fitsFile = pyfits.open(self.filename)
        fitsHeader = fitsFile[0].header
        table_data = fitsFile[0].data

        # RA == y axis/column indices, Dec == x axis/row indices
        self.RA = np.radians(fitsHeader['CRVAL1'])
        self.Dec = np.radians(fitsHeader['CRVAL2'])
        self.xAngScale = np.radians(fitsHeader['CDELT2'])
        self.yAngScale = np.radians(fitsHeader['CDELT2'])
        self.xDecpxl = fitsHeader['CRPIX2']
        self.yRApxl = fitsHeader['CRPIX1']
        self.xExtent = fitsHeader['NAXIS2']
        self.yExtent = fitsHeader['NAXIS1']
        fitsFile.close()

        self.checkSourceData(table_data)
        self.copyData()
        app.displayFitsData()

    def checkSourceData(self, image_table):
        size_of_data = image_table.shape
        dim_of_data = len(size_of_data)

        if dim_of_data == 2:
            self.source_data = image_table
        elif (dim_of_data == 4):
            self.source_data = image_table[0, 0, :, :]
        else:
            sys.exit()

    def copyData(self):
        im_dims.yRA_pxl = self.yRApxl
        im_dims.xDec_pxl = self.xDecpxl
        im_dims.xextent = self.xExtent
        im_dims.yextent = self.yExtent
        shapelet.yRA = self.RA
        shapelet.xDec = self.Dec
        source.extended_source = self.source_data


class Rectangle:
    def __init__(self, x0=0, x1=1, y0=0, y1=1):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1


class ImageDimensions:

    def __init__(self, yRA_pxl=0,  xDec_pxl=0, xextent=100, yextent=100,):
        self.yRA_pxl = yRA_pxl
        self.xDec_pxl = xDec_pxl
        self.xextent = xextent
        self.yextent = yextent

    def resizeImageOfSource(self):
        self.xextent = region.x1 - region.x0
        self.yextent = region.y1 - region.y0
        new_source = np.zeros((self.yextent, self.xextent))

        for i in range(self.yextent):
            for j in range(self.xextent):
                rows = int(i + region.y0)
                cols = int(j + region.x0)
                new_source[i, j] = source.extended_source[rows, cols]

        self.changeRADecRef(region.x0, region.y0)
        self.yRA_pxl = 0
        self.xDec_pxl = 0
        source.extended_source = new_source
        app.displayFitsData()

    def changeRADecRef(self, new_xpxl, new_ypxl):
        if (self.yRA_pxl != new_ypxl):
            shapelet.yRA = shapelet.yRA - (self.yRA_pxl - new_ypxl) * dataset.yAngScale
            self.yRApxl = new_ypxl
        if (self.xDec_pxl != new_xpxl):
            shapelet.xDec = shapelet.xDec - (self.xDec_pxl - new_xpxl) * dataset.xAngScale
            self.xDec_pxl = new_xpxl

    def calcCoordinateList(self):
        xcentre_pxl = m.floor(self.xextent / 2)
        ycentre_pxl = m.floor(self.yextent / 2)
        self.changeRADecRef(xcentre_pxl, ycentre_pxl)
        total_coordinates = self.xextent * self.yextent
        coordinates = np.zeros((total_coordinates, 2))
        xstart = -1 * xcentre_pxl
        ystart = -1 * ycentre_pxl

        k = 0
        for i in range(self.xextent):
            for j in range(self.yextent):
                coordinates[k, 0] = (xstart + i) * abs(dataset.xAngScale)
                coordinates[k, 1] = (ystart + j) * abs(dataset.yAngScale)
                k += 1
        source.pxl_coords_list = coordinates

    def expandArray(self, flat_array):
        twoDarray = np.zeros((self.yextent, self.xextent))
        twoDarray = np.reshape(flat_array, (self.yextent, self.xextent))

        return twoDarray


class SourceImage:
    def __init__(self, extended_source=None, pxl_coords_list=None,
                 max_basis_index=5, normalised_mse=0.0, avg_ssim=1.0):

        self.extended_source = extended_source
        self.pxl_coords_list = pxl_coords_list
        self.max_basis_index = max_basis_index
        self.normalised_mse = normalised_mse
        self.avg_ssim = avg_ssim

    def calcCovarianceFit(self):
        sum_of_source = np.sum(np.sum(self.extended_source))

        weighted_centre = self.calcWeightedCentre(sum_of_source)
        self.pxl_coords_list[:, 0] -= weighted_centre[0]
        self.pxl_coords_list[:, 1] -= weighted_centre[1]
#        im_dims.changeRADecRef(weighted_centre[0], weighted_centre[1])
# Not in original code but maybe should be. move to new centre

        covar_matrix = self.calcFluxMatrix(sum_of_source)
        (eigenval, eigenvect) = self.calcEigenvalueDecomp(covar_matrix)

        self.calcPositionAngle(eigenvect)
        self.calcTransformCoordinates()
#        beta.calcMajorMinor()

    def calcWeightedCentre(self, sum_of_source):
        k = 0
        xcentre = 0.0
        ycentre = 0.0
        for i in range(im_dims.yextent):
            for j in range(im_dims.xextent):
                xcentre += self.extended_source[i, j] * self.pxl_coords_list[k, 0]
                ycentre += self.extended_source[i, j] * self.pxl_coords_list[k, 1]
                k += 1

        xcentre = int(xcentre / sum_of_source)
        ycentre = int(ycentre / sum_of_source)
        weighted_centre = ([xcentre, ycentre])
        return weighted_centre

    def calcFluxMatrix(self, sum_of_source):
        # Matrix terms
        xx_value = 0.
        yy_value = 0.
        xy_value = 0.

        k = 0
        for i in range(im_dims.yextent):
            for j in range(im_dims.xextent):
                xx_value = xx_value + self.extended_source[i, j] * \
                           self.pxl_coords_list[k, 0] * self.pxl_coords_list[k, 0]
                yy_value = xx_value + self.extended_source[i, j] * \
                           self.pxl_coords_list[k, 1] * self.pxl_coords_list[k, 1]
                xy_value = xx_value + self.extended_source[i, j] * \
                           self.pxl_coords_list[k, 0] * self.pxl_coords_list[k, 1]
                k += 1

        xx_value = xx_value / sum_of_source
        yy_value = yy_value / sum_of_source
        xy_value = xy_value / sum_of_source
        covar_matrix = ([xx_value, xy_value, yy_value])
        return covar_matrix

    def calcEigenvalueDecomp(self, covar_matrix):
        eigenmatrix = np.array([[covar_matrix[0], covar_matrix[1]],
                                [covar_matrix[1], covar_matrix[2]]])

        (eigenval, eigenvect) = np.linalg.eig(eigenmatrix)
        return (eigenval, eigenvect)

    def calcPositionAngle(self, eigenvect):
        defaultMajor, defaultMinor = self.calcDefaultValues()
        shapelet.position_angle = -1*m.atan2(eigenvect[1, 1], eigenvect[0, 1])
        shapelet.major = defaultMajor
        shapelet.minor = defaultMinor

    def calcDefaultValues(self):
        x_axis = 0.9 * im_dims.xextent * m.pow(self.max_basis_index, -0.52)
        y_axis = 0.9 * im_dims.yextent * m.pow(self.max_basis_index, -0.52)

        if x_axis > y_axis:
            defaultMajor = x_axis * abs(dataset.xAngScale)
            defaultMinor = y_axis * abs(dataset.yAngScale)
        else:
            defaultMajor = y_axis * abs(dataset.yAngScale)
            defaultMinor = x_axis * abs(dataset.xAngScale)

        return defaultMajor, defaultMinor

    def calcTransformCoordinates(self):
        PA = shapelet.position_angle
        transform_matrix = [[m.cos(PA), -m.sin(PA)], [m.sin(PA), m.cos(PA)]]
        coord_list = np.dot(transform_matrix, np.transpose(self.pxl_coords_list))
        self.pxl_coords_list = np.transpose(coord_list)

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
        sum_source = np.sum(source.extended_source.flatten())
        sqerr = np.square(residual.extended_source.flatten())

        self.normalised_mse = np.sum(sqerr) / (sum_source * sum_source)
        nmse_txt = '%3.4f' % self.normalised_mse
        app.nmse_value.configure(text=nmse_txt)

    def getMinMax(self):
        source_max = np.max(self.extended_source.flatten())
        source_min = np.min(self.extended_source.flatten())

        return source_max, source_min


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
        maxDim = max(im_dims.xextent, im_dims.yextent)
        minResolution = min(dataset.xAngScale, dataset.yAngScale)
        self.minBeta = minResolution
        self.maxBeta = 0.9 * maxDim * m.pow(source.max_basis_index, -0.52) * self.minBeta
        self.previousMajor = shapelet.major
        self.previousMinor = shapelet.minor
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
        shapelet.major = self.previousMajor + 5 * self.stepSize
        if shapelet.major > self.maxBeta:
            shapelet.major = self.maxBeta
        while beta.currentMSE < beta.previousMSE:
            self.calcMajorShapelet()
            self.currentMSE = source.normalised_mse


    def calcMajorShapelet(self):
        self.previousMSE = self.currentMSE
        shapelet.major = shapelet.major - self.stepSize
        if shapelet.major < 0.0:
            shapelet.major = self.maxBeta
        shapelet.calcMoments()
        shapelet.calcShapeletModel()
        shapelet.calcResidual()
        source.calcNMSE()


    def checkMajorSolution(self):
        currentMajor = shapelet.major + self.stepSize
        if round((currentMajor * 10e4)) == round(self.previousMajor * 10e4):
            shapelet.major = self.previousMajor
            self.bestMajor = 1
        else:
            self.previousMajor = shapelet.major + self.stepSize
            shapelet.major = self.previousMajor
            self.bestMinor = 0
            self.bestMajor = 0


    def calcMinor(self):
        shapelet.minor = self.previousMinor + 5 * self.stepSize
        if shapelet.minor > self.maxBeta:
            shapelet.minor = self.maxBeta
        while self.currentMSE < self.previousMSE:
            self.calcMinorShapelet()
            self.currentMSE = source.normalised_mse


    def calcMinorShapelet(self):
        self.previousMSE = self.currentMSE
        shapelet.minor = shapelet.minor - self.stepSize
        if shapelet.minor < 0.0:
            shapelet.minor = self.maxBeta
        shapelet.calcMoments()
        shapelet.calcShapeletModel()
        shapelet.calcResidual()
        source.calcNMSE()


    def checkMinorSolution(self):
        currentMinor = shapelet.minor + self.stepSize
        if round((currentMinor * 10e4)) == round(self.previousMajor * 10e4):
            shapelet.minor = self.previousMinor
            beta.bestMinor = 1
        else:
            self.previousMinor = shapelet.minor + self.stepSize
            shapelet.minor = self.previousMinor
            self.bestMinor = 0
            self.bestMajor = 0


class Shapelet:
    def __init__(self, yRA=0.0, xDec=0.0, moments=[], major=3.0, minor=1.0,
                position_angle=90 ):

        self.yRA = yRA
        self.xDec = xDec
        self.moments = moments
        self.major = major
        self.minor = minor        
        self.position_angle = position_angle

    def calcMoments(self):
        u_coords = source.pxl_coords_list[:,0]/self.major
        v_coords = source.pxl_coords_list[:,1]/self.minor
        
        n_max = source.max_basis_index
        u_basis = self.calcShapeletBasis(u_coords, self.major)
        v_basis = self.calcShapeletBasis(v_coords, self.minor)

        total_moments = int(0.5*(n_max+1)*(n_max+2))
        self.moments = np.zeros((total_moments, 3))

        k=0
        for n1 in range(n_max+1):
            for n2 in range(n_max+1):
                if (n1+n2) <= n_max:
                    basis_function = np.multiply(u_basis[n1, :], v_basis[n2, :])
                    a_moment = (self.calcLeastSquares(basis_function))
                    self.moments[k,0] = n1
                    self.moments[k,1] = n2 
                    self.moments[k,2] = a_moment
                    k += 1
                    
    def calcShapeletBasis(self, coords, beta):
        total_coords = coords.shape[0]
        sq_coords = np.power(coords, 2)
        gauss = np.exp(-0.5*sq_coords)
        norm_const1 = m.sqrt(beta*m.sqrt(m.pi))
        shapelet_basis = np.zeros((source.max_basis_index+1, total_coords))
                       
        for n in range(source.max_basis_index+1):
            norm_const2 = m.sqrt(m.pow(2, n)*m.factorial(n))
            norm_constant = 1./(norm_const1*norm_const2)
            hermite_z = self.calcHermite(n, coords)
            shapelet_basis[n, :] = norm_constant*np.multiply(hermite_z, gauss)
            
        return shapelet_basis

    def calcHermite(self, n, coords):
        hermites = np.loadtxt('hermite_coeffs.txt')

        k=0
        hermite_poly=0.0
        while (k <= n):
            hermite_poly = hermite_poly+hermites[n, k]*np.power(coords, (n-k))
            k += 1
        return hermite_poly

    def calcLeastSquares(self, basis_function):
        flat_source = source.extended_source.flatten()
        transpose_basis = np.transpose(basis_function)
        demoninator = np.dot(transpose_basis, basis_function)
        numerator = np.dot(transpose_basis, flat_source)
        a_moment = numerator/demoninator
        return a_moment

    def saveShapelet(self):
        file_name = '../Models/' + dataset.source_name + '_moments.txt'
        np.savetxt(file_name, self.moments)

        beta1 = np.degrees(self.major)
        beta2 = np.degrees(self.minor)
        shapelet_parameters = ([np.degrees(self.yRA), np.degrees(self.xDec), beta1,
                                beta2, np.degrees(self.position_angle)])

        file_name = '../Models/' + dataset.source_name + '_parameters.txt'
        np.savetxt(file_name, shapelet_parameters)

    def calcShapeletModel(self):
        flat_model = np.zeros((model_dims.xextent * model_dims.yextent, 1))

        total_moments = self.moments.shape[0]

        u_coords = model.pxl_coords_list[:,0]/self.major
        v_coords = model.pxl_coords_list[:,1]/self.minor
        u_basis = self.calcShapeletBasis(u_coords, self.major)
        v_basis = self.calcShapeletBasis(v_coords, self.minor)
        
        for k in range(total_moments):
            n1 = int(self.moments[k, 0])
            n2 = int(self.moments[k, 1])
            a_moment = self.moments[k, 2]
            basis_function = np.multiply(u_basis[n1, :], v_basis[n2, :])
            flat_model = np.add(flat_model, a_moment * basis_function)

        model.extended_source = flat_model
#        model.extended_source = model_dims.expandArray(flat_model)

    def calcResidual(self):
        for i in range(im_dims.yextent):
            for j in range(im_dims.xextent):
                residual.extended_source[i, j] = source.extended_source[i, j] - \
                                             model.extended_source[i, j]

    def getShapelet(self):
        file_name = '../Models/' + dataset.source_name + '_moments.txt'
        self.moments = np.loadtxt(file_name)

        file_name = '../Models/' + dataset.source_name + '_parameters.txt'
        shapelet_parameters = np.loadtxt(file_name)

        self.yRA = shapelet_parameters[0]
        self.xDec = shapelet_parameters[1]
        self.major = shapelet_parameters[2]
        self.minor = shapelet_parameters[3]
        self.position_angle = shapelet_parameters[4]


root = Tk()
shapelet = Shapelet()
source = SourceImage()
model = SourceImage()
residual = SourceImage()
dataset = Fits_data()
region = Rectangle()
im_dims = ImageDimensions()
model_dims = ImageDimensions()
residual_dims = ImageDimensions()
beta = Beta()

app = ShapeletGUI(root)

root.mainloop()

