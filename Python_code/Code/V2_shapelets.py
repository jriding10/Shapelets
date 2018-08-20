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
    def __init__(self, centreRA=0.0, centreDec=0.0, extendedSource=[], \
    yRApxl=0, xDecpxl=0, pxlCoordList = [], xAngularScale=1.0, \ 
    yAngularScale=1.0, xExtent=100, yExtent=100, yRApxl=0, xDecpxl=0, \
    moments=[], majorAxis=3.0, minorAxis=1.0, positionAngle=90, \
    maxBasisIndex=5, normalisedMSE=0.0, averageSSIM=1.0): 
    
        self.centreRA = centreRA
        self.centreDec = centreDec
        self.extendedSource = extendedSource
        self.yRApxl = yRApxl
        self.xDecpxl = xDecpxl
        self.pxlCoordinateList = pxlCoordList
        self.xAngularScale = xAngularScale      
        self.yAngularScale = yAngularScale
        self.xExtent = xExtent
        self.yExtent = yExtent
        self.moments = moments
        self.majorAxis = majorAxis
        self.minorAxis = minorAxis        
        self.positionAngle = positionAngle        
        self.maxBasisIndex = maxBasisIndex
        self.normalisedMSE = normalisedMSE
        self.averageSSIM = averageSSIM

#------------------------------------------------------------------------#    
    def calcNumberMoments(self):
        totalNumMoments = 0.5*(self.maxBasisIndex + 1)* \
        (self.maxBasisIndex + 2)
        return totalNumMoments

#------------------------------------------------------------------------#        
    def calcMoments(self):
        uCoords = self.pxlCoordList[:,0]/self.major
        vCoords = self.pxlCoordList[:,1]/self.minor
             
        uShapeletBasis = calcShapeletBasis(uCoords, self.major)
        vShapeletBasis = calcShapeletBasis(vCoords, self.minor)
        flatSourceArray = flattenExtendedSourceArray()
        
        self.moments=np.zeros((self.maxBasisIndex+1,self.maxBasisIndex+1) 
        
        k = 0
        for n1 in range( self.maxBasisIndex ):
            for n2 in range( self.maxBasisIndex ):
                if (n1+n2) <= self.maxBasisIndex:
                    basisFunction = np.multiple((uShapeletBasis[n1,:],\
                    vShapeletBasis[n2, :]), axis=1)
                    self.moments[n1,n2]=calcLeastSquares(basisFunction, \
                    flatSourceArray)
                
#------------------------------------------------------------------------#
    def calcShapeletBasis( self, zCoords, beta ):
        zGauss = np.exp(0.5*np.array(zCoords)**2)
        normConst1 = m.sqrt(m.pow(beta,2)*m.sqrt(m.pi))
        
        zHermiteBasis = calcHermiteFunction( zCoords )
                       
        for n in range( self.maxBasisIndex ):
            normConst2 = m.sqrt(m.pow(2, n)*m.factorial(n))
            normConstant = 1./(normConst1*normConst2)
            shapeletBasis = normConstant*zGauss*hermite[n,:]
            
        return shapeletBasis  
            
#------------------------------------------------------------------------#
   def calcHermiteFunction( zCoords ):
       hermitePolynomials = np.loadtxt('hermite_coeffs.txt')
       for n in range( self.maxBasisIndex ):
            normConstant = 1./normConst2
            k = 0
            hermite = 0.
            while (k <= n ):
                hermite[n, :] += hermitePolynomials[n, k]* \
                m.pow((np.array(zCoords)),(n-k))
                k += 1
            
        return hermite           

#------------------------------------------------------------------------#
    def calcLeastSquares(self, basisFunction, flatSourceArray ):
        transposeBasis = np.transpose(basisFunction)
        demoninator = np.dot(transposeBasis, basisFunction)
        numerator = np.dot(transposeBasis, flatSourceArray)
        aMoment = numerator/demoninator
        return aMoment 
                
#------------------------------------------------------------------------#    
    def getFitsFile(self, filename):
        fitsFile = pyfits.open( filename )
        fitsHeader = fitsFile[0].header
        self.extendedSource = fitsFile[1].data
        
        # RA == y axis/column indices, Dec == x axis/row indices
        RA = np.radians(fitsHeader['CRVAL1'])
        Dec = np.radians(fitsHeader['CRVAL2'])
        self.xAngularScale = np.radians(fitsHeader['CDELT2'])
        self.yAngularScale = np.radians(fitsHeader['CDELT2'])
        self.xDecpxl = fitsHeader['CRPIX2']
        self.yRApxl = fitsHeader['CRPIX1']
        self.xExtent = fitsHeader['NAXIS2']
        self.yExtent = fits_info['NAXIS1']
        fitsFile.close()
        
        # Reference the extended source to the pxl (0,0)   
        changeRADecRef( 0, 0 ) 
                

#------------------------------------------------------------------------#
    def resizeImageOfSource(self, newXstart, newYstart, newXend, newYend): 
        self.xExtent = newXend-newXstart
        self.yExtend = newYend-newYstart    

        for i in range( self.xExtent ):
            for j in range( self.yExtent ):
                newSource[i, j] = self.extendedSource[i+newXstart, \
                j+newYstart]         
         
        changeRADecRef( newXstart, newYstart )      
        self.yRApxl = 0
        self.xDecpxl = 0
        
#------------------------------------------------------------------------#       
    def changeRADecRef( self, newXpxl, newYpxl ):       
        if (self.yRApxl ne newYpxl):
            self.RA = self.RA-(self.yRApxl-newYpxl)*self.yAngularScale
            self.yRApxl = newYpxl
        if (self.xDecpxl ne newXpxl):
            self.Dec = self.Dec-(self.xDecpxl-newXpxl)*self.xAngularScale
            self.xDecpxl = newXpxl

#------------------------------------------------------------------------#
   def calcCoordinateList(self, xCentrePxl, yCentrePxl):
       changeRADecRef( xCentrePxl, yCentrePxl )
       xStart = -1*xCentrePxl
       yStart = -1*yCentrePxl
       k = 0
       for i in range( self.xExtent ):
           for j in range( self.yExtent ):
                coordinates[k, 0] = xStart + i
                coordinates[k, 1] = yStart + j
                k += 1
       self.pxlCoordList = coordinates

#------------------------------------------------------------------------#
   def flattenExtendedSourceArray( self ):
       k = 0
       for i in range( self.xExtent ):
           for j in range( self.yExtent ):
               flatSourceArray[k] = self.extendedSource[i, j]
               k += 1
       return flatSourceArray

#------------------------------------------------------------------------#
   def calcCovarianceFit(self):
       totalFluxOfSource = np.sum( np.sum( self.extendedSource ))
       numberOfPxls = self.xExtent*self.yExtent
       
       fluxWeightedCentre=calcWeightedCentre(totalFluxOfSource)
       changeRADecRef( fluxWeightedCentre[0], fluxWeightedCentre[1] )
       
       fluxMatrix = calcFluxMatrix( totalFluxOfSource )
       (eigenval, eigenvect) = calcEigenValueDecomp( fluxMatrix )
       
       self.major = np.sqrt(eigenval[1])
       self.minor = np.sqrt(eigenval[0])
       self.positionAngle = m.atan2(eigenvect[1,1],eigenvect[0,1])
       
#------------------------------------------------------------------------#       
    def calcWeightedCentre(self, totalFluxOfSource ):
        k = 0
        xCentre = 0.
        yCentre = 0.
        for i in range( xExtent ):
            for j in range( yExtent ):
                xCentre += self.extendedSource[i,j]*self.pxlCoordList[k,0]
                yCentre += self.extendedSource[i,j]*self.pxlCoordList[k,1]
                k += 1
       
        xCentre = int(xCentre/totalFluxOfSource)
        yCentre = int(yCentre/totalFluxOfSource)
        return fluxWeightedCentre = [xCentre, yCentre]
        
#------------------------------------------------------------------------#
    def calcFluxMatrix(self, totalFluxOfSource):
       # Matrix terms
       xxFlux = 0.
       yyFlux = 0.
       xyFlux = 0.
       
       k = 0
       for i in range( self.xExtent ):
           for j in range( self.yExtent ):              
               xxFlux = xxFlux + self.extendedSource[i,j]* \
               self.pxlCoordList[k, 0]*self.pxlCoordList[k,0]
               yyFlux = xxFlux + self.extendedSource[i,j]* \
               self.pxlCoordList[k, 1]*self.pxlCoordList[k,1]       
               xyFlux = xxFlux + self.extendedSource[i,j]* \
               self.pxlCoordList[k, 0]*self.pxlCoordList[k,1] 
               k += 1
        
       xxFlux = xxFlux/totalFluxOfSource
       yyFlux = yyFlux/totalFluxOfSource
       xyFlux = xyFlux/totalFluxOfSource
       return fluxMatrix = [xxFlux, xyFlux, yyFlux]
       
#------------------------------------------------------------------------#
    def calcEigenValueDecomp(self, fluxMatrix ):
        trueFluxMatrix = np.array([[fluxMatrix[0],fluxMatrix[1]], \
        [fluxMatrix[1],fluxMatrix[2]]]) 
        
        (eigenval, eigenvect) = np.linalg.eig(trueFluxMatrix)
        return (eigenval, eigenvect)
       
#------------------------------------------------------------------------#
    def calcTransformCoordinates( self ):
        PA = self.positionAngle
        transformMatrix=[[m.cos(PA), -m.sin(PA)],[m.sin(PA), m.cos(PA)]]
        coordList=np.dot(transformMatrix,np.transpose(self. pxlCoordList))
        self.pxlCoordList = np.transpose(coordList)
          
        
##########################################################################       
#class ShapeletGUI:



