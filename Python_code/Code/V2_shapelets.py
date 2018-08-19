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

import math as m
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import pyfits

## Load gui modules - all of PyQt4 for python 2.7 and selected parts of
## PyQt5 for python 3

try:
    from PyQt4.QtGui import *
    from PyQt4.QtCore import *
except:
    try:
        from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, \    
            QHBoxLayout, QVBoxLayout, QLabel, QLineEdit, QTextEdit, \
            QComboBox, QCheckBox
        from PyQt5.QtCore import pyqtSlot
        from PyQt5.QtGui import QImage, QPixmap
    except:
        print("Please install PyQt4 or PyQt5.")
        raise UserWarning

##########################################################################
## The Shapelet Class: A shapelet is defined by its major axis, minor   ##
## axis and its moments. To give meaning to an astronomical source, a   ##
## (RA, Dec) and position angle are also specified.                     ##
##                                                                      ##
## The final two parameters in this class measure the goodness-of-fit to## 
## the original data used to create the moments. These are the          ##    
## normalised mean square error and the structial similarly index       ## 
## measure. The nMSE is common in astronomy but the SSIM is better      ##
## suited at accurately measuring how well an extracted image/model     ##
## "looks like" the data.                                               ##
##########################################################################

class Shapelet:
    def __init__(self, centreRA=0.0, centreDec=0.0, moments=[], \
    majorAxis=3.0, minorAxis=1.0, positionAngle=90, maxBasisIndex=5, \ 
    normalisedMSE=0.0, averageSSIM=1.0) 
        self.centreRA = centreRA
        self.centreDec = centreDec
        self.maxBasisIndex = maxBasisIndex
        self.majorAxis = majorAxis
        self.minorAxis = minorAxis
        self.positionAngle = positionAngle
        self.moments = moments
        self.normalisedMSE = normalisedMSE
        self.averageSSIM = averageSSIM
    
    def calcNumberMoments(self):
        totalNumMoments = 0.5*(self.maxBasisIndex + 1)*(self.maxBasisIndex + 2)
        return totalNumMoments
        
    def getMoments(self, filename, sourceName): 
        hermite_coeffs_file = 'hermite_coeffs.txt'    
        if (needCalcMoments):
            hermitePolynomialCoeffs = np.loadtxt( filename )
            self.moments = decomposeSource( sourceName.imageOfSource )            
        else:
            self.moments = np.loadtxt( filename )
         
##########################################################################
## The Image Class: An image can be the original data or the reproduced ##
## model. What defines it is a matrix of flux (Jy.beam^-1) and its      ##
## dimensions - angular scale in both x and y, extend in both x and y,  ##
## and a matrix/pixel/coordinate set that maps (RA, Dec) to (x, y).     ##
## Pixels are used by the fits file in the header information.          ##
##                                                                      ##
##          matrix(row, column) = (pxl2, pxl1) = (Dec, RA).             ##
##                                                                      ##
##########################################################################       
    
class ExtendedSource:
    def __init__(self, imageOfSource=[], RA=0.0, Dec=0.0, yRApxl=0, \
    xDecpxl=0, xAngularScale=1.0, yAngularScale=1.0, xExtent=100, \
    yExtent=100, yRApxl=0, xDecpxl=0)
        self.imageOfSource = imageOfSource
        self.RA = RA
        self.Dec = Dec 
        self.yRApxl = yRApxl
        self.xDecpxl = xDecpxl
        self.xAngularScale = xAngularScale      
        self.yAngularScale = yAngularScale   
    
    def getFitsFile(self, filename):
        fitsFile = pyfits.open( filename )
        fitsHeader = fitsFile[0].header
        self.imageOfSource = fitsFile[1].data
        
        # RA == y axis/column indices, Dec == x axis/row indices
        self.RA = np.radians(fitsHeader['CRVAL1'])
        self.Dec = np.radians(fitsHeader['CRVAL2'])
        self.xAngularScale = np.radians(fitsHeader['CDELT2'])
        self.yAngularScale = np.radians(fitsHeader['CDELT2'])
        self.xDecpxl = fitsHeader['CRPIX2']
        self.yRApxl = fitsHeader['CRPIX1']
        self.xExtent = fitsHeader['NAXIS2']
        self.yExtent = fits_info['NAXIS1']
        fitsFile.close()


    def resizeImageOfSource(self, newXstart, newYstart, newXend, newYend):
        changeRADecRef( newXpxl, newYpxl)      
        self.xExtent = newXend-newXstart
        self.yExtend = newYend-newYstart
        self.yRApxl = 0
        self.xDecpxl = 0
       
        # change reference so that pxl (0,0) occurs at image centre
        centreXpxl = floor(self.xExtent/2)
        centreYpxl = floor(self.yExtent/2)
        changeRADecRef( centreXpxl, centreYpxl )
        self.yRApxl = 0
        self.xDecpxl = 0       
        
        k = 0
        for i in range( newXstart, newXend+1 ):
            for j in range( newYstart, newYend+1 ):
                newImageOfSource[k, :] = self.imageOfSource[i, j]
                k = k+1          
        self.imageOfSource = newImageOfSource
       
    def changeRADecRef(self, newXpxl, newYpxl ):
        if (self.yRApxl ne newYpxl):
            self.RA = self.RA - (self.yRApxl - newYpxl)*self.yAngularScale
            self.yRApxl = newXpxl
        if (self.xDecpxl ne newXpxl):
            self.Dec = self.Dec - (self.xDecpxl - newXpxl)*self.xAngularScale
            self.xDecpxl = newXpxl

   def coordinateList(self):
       for i in range( self. 
   def covarianceFit(self, shapeletName)
           
    S = sum(sum(data))
    nside = m.sqrt(coords.shape[0])
    nside = int(nside)
    npix = nside*nside
    offsets = (0.0,0.0)
    addmat = np.zeros((npix,2))       
       
       
class ShapeletGUI:

# Load fits file
sourceName = ExtendedSource()
getFitsFile(filename, sourceName )

# Resize extended source
sourceName.resizeImageOfSource(newXstart, newYstart, newXend, newYend)

# Decompose into moments
shapeletName = Shapelet()
shapeletName.maxBasisIndex = 25
shapeletName.calcNumberMoments()
pxlCoordinateList = sourceName.coordinateList()

sourceName.covarianceFit( shapeletName )

