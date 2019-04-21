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

# Custom classes
import Beta
import Shapes
import Fits
import Coords
import Source

# Standard libraries
import sys
import math as m
import numpy as np
import matplotlib.pyplot as plt

# GUI libraries
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

main_window_width = 1000
main_window_height = 1000
im_size = (600, 600)
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
        self.setCentralWidget(self.main_gui)

        self.show()


class MainWindow(QWidget):

    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout()

        # initialise tab
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()
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
        self.canvas3 = DisplayResults(self.tab3)
        self.button_save = QPushButton('Save')
        self.button_save2 = QPushButton('Save')
        self.label_ssim = QLabel('SSIM')
        self.data_ssim = QLabel()
        self.label_nmse = QLabel('NMSE')
        self.data_nmse = QLabel()
        self.button_moments = QPushButton('Moment Map')
        self.label_major = QLabel('Major (arcmin)')
        self.data_major = QLabel()
        self.label_minor = QLabel('Minor (arcmin)')
        self.data_minor = QLabel()
        self.label_pa = QLabel('Position Angle (degrees)')
        self.data_pa = QLabel()
        self.label_min_moments = QLabel('How many moments?')
        self.entry_min_moments = QLineEdit()
        self.label_default = QLabel('Default is all')

        #create tabs
        self.tabs.addTab(self.tab1, "Data Selection")
        self.tabs.addTab(self.tab2, "Numerical Results")
        self.tabs.addTab(self.tab3, "Results")

        self.createTab1()
        self.createTab2()
        self.createTab3()
        self.show()

        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

    def createTab1(self):
        layout = QVBoxLayout(self)
        self.button_file_select.clicked.connect(lambda: self.findFitsFile())
        layout.addWidget(self.button_file_select)
        layout.addWidget(self.canvas)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        hbox1 = QWidget()
        box1 = QGridLayout(self)
        self.button_undo.clicked.connect(lambda: self.undoROI())
        self.button_resize.clicked.connect(lambda: self.createROI())
        self.button_calculate.clicked.connect(lambda: calcShapelet())
        box1.addWidget(self.button_undo, 0, 0)
        box1.addWidget(self.button_resize, 0, 1)
        box1.addWidget(self.button_calculate, 0, 3)
        box1.addWidget(self.label_name, 1, 0)
        box1.addWidget(self.entry_name, 1, 1)
        box1.addWidget(self.label_nbasis, 1, 2)
        box1.addWidget(self.entry_nbasis, 1, 3)
        box1.addWidget(self.label_min_moments, 2, 0)
        box1.addWidget(self.entry_min_moments, 2, 1)
        box1.addWidget(self.label_default, 2, 2)

        hbox1.setLayout(box1)
        layout.addWidget(hbox1)
        self.tab1.setLayout(layout)

    def createTab2(self):
        layout = QVBoxLayout(self)
        hbox1 = QWidget()
        box1 = QHBoxLayout(self)
        box1.addWidget(self.label_major)
        box1.addWidget(self.data_major)
        box1.addWidget(self.label_minor)
        box1.addWidget(self.data_minor)
        box1.addWidget(self.label_pa)
        box1.addWidget(self.data_pa)
        hbox1.setLayout(box1)
        layout.addWidget(hbox1)
        layout.addWidget(self.canvas2)
        self.button_save.clicked.connect(lambda: shapelet.saveMoments())
        layout.addWidget(self.button_save)
        self.tab2.setLayout(layout)

    def createTab3(self):
        layout = QVBoxLayout(self)
        hbox1 = QWidget()
        box1 = QHBoxLayout(self)
        box1.addWidget(self.label_nmse)
        box1.addWidget(self.data_nmse)
        box1.addWidget(self.label_ssim)
        box1.addWidget(self.data_ssim)
        hbox1.setLayout(box1)
        layout.addWidget(hbox1)
        layout.addWidget(self.canvas3)
        self.button_save2.clicked.connect(lambda: saveResults())
        layout.addWidget(self.button_save2)
        self.tab3.setLayout(layout)

    @pyqtSlot()
    def findFitsFile(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filetypes = 'Fits files (*.fits);;All files (*)'
        fits_data.filename, _ = QFileDialog.getOpenFileName(self, 'Data', myDirectory, filetypes, options=options)
        fits_data.getFitsFile()
        data.source = fits_data.source_data
        self.canvas.displayData(fits_data.source_data)

    def on_press(self, event):
        coords.row0 = int(event.ydata)
        coords.column0 = int(event.xdata)

    def on_release(self, event):
        coords.row1 = int(event.ydata)
        coords.column1 = int(event.xdata)
        coords.row_extent = abs(coords.row1 - coords.row0)
        coords.column_extent = abs(coords.column1 - coords.column0)
        self.canvas.displayRectangle()

    def createROI(self):
        old_ref_pxl = [fits_data.rowDecpxl, fits_data.colRApxl]
        old_ref_coords = [fits_data.RA, fits_data.Dec]
        data.angScale = [fits_data.xAngScale, fits_data.yAngScale]
        data.ref_pxl = [0, 0]
        data.ref_coords = coords.selectROI(data.source, old_ref_pxl, old_ref_coords, data.angScale)
        data.source = coords.data
        self.canvas.displayData(data.source)

    def undoROI(self):
        radio_source.extended_source = fits_data.source_data
        data.source = fits_data.source_data
        data.ref_coords = [fits_data.RA, fits_data.Dec]
        data.ref_pxl = [fits_data.rowDecpxl, fits_data.colRApxl]
        self.canvas.displayData(radio_source.extended_source)

    def getUserInputs(self):
        name = self.entry_name.text()
        if len(name) != 0:
            data.name = name
        else:
            data.name = 'test'

        basis = self.entry_nbasis.text()
        if len(basis) != 0:
            radio_source.max_basis_index = int(basis)
            data.nmax = int(basis)
        else:
            radio_source.max_basis_index = 25

        n = radio_source.max_basis_index
        min_moments = self.entry_min_moments.text()
        if len(min_moments) != 0:
            shapelet.min_moments = int(min_moments)
        else:
            shapelet.min_moments = int(0.5*(n+1)*(n+2))

    def updateDisplay(self):
        qApp.processEvents()

    def displayResults(self):
        moment_array = shapelet.makeMomentArray()
        self.canvas2.displayData(moment_array)
        self.canvas3.displayData()


class DisplayResults(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure(figsize=im_size, dpi=10)
        self.axes1 = fig.add_subplot(221)
        self.axes2 = fig.add_subplot(222)
        self.axes3 = fig.add_subplot(223)
        self.displayDefault()

        FigureCanvas.__init__(self, fig)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.setParent(parent)

    def displayDefault(self):
        data = np.zeros((int(main_window_width/4), int(main_window_height/4)))
        model = np.zeros((int(main_window_width/4), int(main_window_height/4)))
        resids = np.zeros((int(main_window_width/4), int(main_window_height/4)))
        self.axes1.imshow(data)
        self.axes2.imshow(model)
        self.axes3.imshow(resids)

    def displayData(self):
        data = radio_source.extended_source
        model = model_source.extended_source
        resids = residual.extended_source

        self.axes1.imshow(data, cmap='hot')
        self.axes2.imshow(model, cmap='hot')
        self.axes3.imshow(resids, cmap='hot')


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
        delta_column = abs(coords.column0 - coords.column1)
        rect = plt.Rectangle((coords.column0, coords.row0), delta_column, delta_row,
                             linewidth=2, edgecolor='y', facecolor='none')
        self.axes.add_patch(rect)
        self.draw()


class Common:
    def __init__(self):
        self.coords_list = None
        self.source = None
        self.model = None
        self.resids = None
        self.extent = None
        self.angScale = None
        self.ref_pxl = None
        self.ref_coords = None
        self.nmax = 25
        self.name = 'test'


def calcShapelet():
    ex.main_gui.getUserInputs()
    setupSource()
    setupCoords()
    calcBetas()
    shapelet.calcMoments(data.source, data.coords_list, data.nmax)
    shapelet.saveShapelet(data.name)
    setupModelSource()
    calcModel()
    calcPerformance()
    ex.main_gui.displayResults()
    ex.main_gui.updateDisplay()

def setupSource():
    data.angScale = [fits_data.xAngScale, fits_data.yAngScale]
    old_ref_pxl = [fits_data.rowDecpxl, fits_data.colRApxl]
    old_ref_coords = [fits_data.RA, fits_data.Dec]
    if coords.row1 == None:
        data.extent = [fits_data.xExtent, fits_data.yExtent]
        data.ref_coords = coords.useWholeImage(data.extent, fits_data.source_data, old_ref_pxl, old_ref_coords, data.angScale)

    data.source = coords.data
    radio_source.extended_source = coords.data
    row_centre_pxl = m.floor(coords.row_extent/2)
    col_centre_pxl = m.floor(coords.column_extent/2)
    data.ref_pxl = [row_centre_pxl, col_centre_pxl]
    old_ref_pxl = [0, 0]
    old_ref_coords = data.ref_coords
    data.ref_coords = coords.changeRADecRef(old_ref_pxl, data.ref_pxl, old_ref_coords, data.angScale)
    data.coords_list = coords.calcCoordinateList(data.angScale)

def setupCoords():
    moveCentreToFluxWeighted()
    rotateCoordinates()

def calcBetas():
    minResolution = min(data.angScale[0], data.angScale[1])
    beta1, beta2 = beta.calcAllBetas(data.source, data.extent, minResolution)
    if beta1 > beta2:
        shapelet.major = beta1
        shapelet.minor = beta2
    else:
        shapelet.major = beta2
        shapelet.minor = beta1
    major = np.degrees(shapelet.major) * 60
    minor = np.degrees(shapelet.minor) * 60
    major_text = '%3.4f' % major
    minor_text = '%3.4f' % minor
    ex.main_gui.data_major.setText(major_text)
    ex.main_gui.data_minor.setText(minor_text)

def setupModelSource():
    model_source.__dict__ = radio_source.__dict__.copy()
    residual.__dict__ = radio_source.__dict__.copy()

def calcModel():
    flat_model = shapelet.calcShapeletModel(data.coords_list, shapelet.moments)
    data.model = coords.expandArray(flat_model)
    data.resids = shapelet.calcResidual(data.source, data.model)

def calcPerformance():
    flat_data = data.source.flatten()
    flat_model = data.model.flatten()
    flat_resids = data.resids.flatten()
    nmse = radio_source.calcNMSE(flat_data, flat_resids)
    ssim = radio_source.calcSSIM(flat_data, flat_model)

    nmse_txt = '%E' % nmse
    ex.main_gui.data_nmse.setText(nmse_txt)
    ssim_txt = '%3.4f' % ssim
    ex.main_gui.data_ssim.setText(ssim_txt)

def saveResults():
    plt.imsave(data.name + '_data.png', radio_source.extended_source)
    plt.imsave(data.name + '_model.png', model_source.extended_source)
    plt.imsave(data.name + '_resids.png', residual.extended_source)

def moveCentreToFluxWeighted():
    weighted_centre = radio_source.calcWeightedCentre(data.coords_list)
    data.coords_list[:, 0] -= weighted_centre[0]
    data.coords_list[:, 1] -= weighted_centre[1]
    new_row_pxl = data.ref_pxl - weighted_centre[0]
    new_col_pxl = data.ref_pxl - weighted_centre[1]
    new_ref_pxl = [new_row_pxl, new_col_pxl]
    data.ref_coords = coords.changeRADecRef(data.ref_pxl, new_ref_pxl, data.ref_coords, data.angScale)
    data.ref_pxl = new_ref_pxl
    shapelet.colRA = data.ref_coords[0]
    shapelet.rowDec = data.ref_coords[1]

def rotateCoordinates():
    covar_matrix = radio_source.calcFluxMatrix(data.extent, data.coords_list)
    eigenvalue, eigenvector = radio_source.calcEigenvalueDecomp(covar_matrix)
    position_angle = radio_source.calcPositionAngle(eigenvector)
    shapelet.position_angle = position_angle
    data.coords_list = coords.calcTransformCoordinates(position_angle, data.coords_list)
    pa_text = '%3.5f' % np.degrees(position_angle)
    ex.main_gui.data_pa.setText(pa_text)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    data = Common()
    fits_data = Fits.FitsData()
    radio_source = Source.ImageData()
    model_source = Source.ImageData()
    residual = Source.ImageData()
    coords = Coords.Coordinates()
    shapelet = Shapes.Shapelets()
    beta = Beta.Beta()
    ex = App()
    sys.exit(app.exec_())
