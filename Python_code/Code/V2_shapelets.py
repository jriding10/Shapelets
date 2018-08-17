#! /bin/bash python
# Shapelets.py
# A remake of Shapelets.py based upon what I now know. Decomposes an image
# into moments using gauss hermite polynomials as a basis - shapelets to 
# astronomers -> f(x,y) = A.g(a, b) where f(x, y) is the original function
# (image), g(a, b) is the new basis (shapelets) and A is the transform 
# (or moments).
#
#
# Version History 
#
# v1.1  2018/08/18  JLR     Rewrite of original code  
#
#
##########################################################################

import math as m
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from PyQt4.QtGui import *
from PyQt4.QtCore import *

class Shapelet:
    RA = 0.0
    Dec = 0.0
    maxBasisIndex = 5
    majorAxis = 3.0
    minorAxis = 1.0
    positionAngle = 90
    moments = []
    normalisedMSE = 0.0
    averageSSIM = 1.0
    
    def calcNumberMoments(self):
    
        totalNumMoments = 0.5*(self.maxBasisIndex + 1)*(self.maxBasisIndex + 2)
        return totalNumMoments
        
    def getMoments(self, filename, extendedSource):
        
        if (needCalcMoments):
            hermitePolynomialCoeffs = np.loadtxt( filename )
            sourceMoments = decomposeSource( extendedSource )            
        else:
            sourceMoments = np.loadtxt( filename )
         return sourceMoments
         
     
    
class Image:
    extendedSource = []
    xAngularScale = 1.0
    yAngularScale = 1.0
    xExtent = 100.0
    yExtent = 100.0
    RA = 0.0
    Dec = 0.0
    yRApxl = 0.0
    xDecpxl = 0.0
    
    def getExtendedSource(self, filename ):
        
        fitsFile = fits.open( filename )
        fitsHeader = fitsFile[0].header
        extendedSource = fitsFile[1].data
        
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

        return extendedSource

   def resizeExtendedSource(self, newXstart, newYstart, newXend, newYend):
       
       changeRADecRef()
       self.RA = newYstart*yAngularScale + self.RA
       self.Dec = newXstart*xAngularScale + self.Dec
       
       self.xExtent = newXend-newXstart
       self.yExtend = newYend-newYstart
       
       k = 0
       for i in range( newXstart, newXend+1 ):
           for j in range( newYstart, newYend+1 ):
               newExtendedSource[k, :] = self.extendedSource[i, j]
               k = k+1
           
       self.extendedSource = newExtendedSource
       return
           
       
       
       
       
       
class ShapeletGUI:


