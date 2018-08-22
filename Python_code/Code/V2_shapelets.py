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
import astropy.io.fits as pyfits
import shelve
from scipy.special import eval_hermite

## Load gui modules - all of PyQt4 for python 2.7 and selected parts of
## PyQt5 for python 3

try:
    from PyQt4.QtGui import *
    from PyQt4.QtCore import *
except:
    try:
        from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton,    
            QHBoxLayout, QVBoxLayout, QLabel, QLineEdit, QTextEdit,
            QComboBox, QCheckBox)
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
## I have left the reference to yRA and xDec (instead of just RA and    ##
## Dec) since I keep mixing the two around.                             ##
##########################################################################

class Shapelet:
    def __init__(self, yRA=0.0, xDec=0.0, moments=[], major=3.0, minor=1.0,
                positionAngle=90):  

        self.yRA = yRA
        self.xDec = xDec
        self.moments = moments
        self.major = major
        self.minor = minor        
        self.positionAngle = positionAngle        
        
#------------------------------------------------------------------------#
    def saveShapelet(self, shapelet):
        filename = '../Models/' + shapelet + '_moments.txt'
        np.savetxt(filename, self.moments)
        
        toBeSaved = ([self.yRA, self.xDec, self.major, 
                    self.minor, self.positionAngle])
                    
        filename = '../Models/' + shapelet + '_parameters.txt'
        np.savetxt(filename, toBeSaved)
        
#------------------------------------------------------------------------#  
    def getShapelet(self, shapelet):
        filename = '../Models/' + shapelet + '_moments.txt'
        self.moments = np.loadtxt(filename)
        
        filename = '../Models/' + shapelet + '_parameters.txt'
        toBeDivided = np.loadtxt(filename)

        self.yRA = toBeDivided[0]
        self.xDec = toBeDivided[1]
        self.majorAxis = toBeDivided[2]
        self.minorAxis = toBeDivided[3]
        self.positionAngle = toBeDivided[4]
        
#------------------------------------------------------------------------# 
    def calcMoments( self, sourceName ):
        uCoords = sourceName.pxlCoordList[:,0]/self.major
        vCoords = sourceName.pxlCoordList[:,1]/self.minor
        
        n_max = sourceName.maxBasisIndex     
        uShapeletBasis = self.calcShapeletBasis(uCoords, self.major, n_max)
        vShapeletBasis = self.calcShapeletBasis(vCoords, self.minor, n_max)
        flatSourceArray = sourceName.flattenExtendedSourceArray()

        totalMoments = int(0.5*(n_max+1)*(n_max+2))
        self.moments = np.zeros((totalMoments, 3)) 
        
        k=0
        for n1 in range(n_max+1 ):
            for n2 in range( n_max+1 ):
                if (n1+n2) <= n_max:
                    basisFunction=uShapeletBasis[n1,:]*vShapeletBasis[n2,:]
                    aMoment = (self.calcLeastSquares(basisFunction,
                                          flatSourceArray))
                    self.moments[k,0] = n1
                    self.moments[k,1] = n2 
                    self.moments[k,2] = aMoment
                    k += 1
                    
#------------------------------------------------------------------------#
    def calcShapeletBasis( self, zCoords, beta, n_max ):
        totalCoords = zCoords.shape[0]
        zGauss = np.exp(0.5*np.power(zCoords,2))
        normConst1 = m.sqrt(beta*m.sqrt(m.pi))
        shapeletBasis = np.zeros(( n_max+1, totalCoords )) 
                       
        for n in range( n_max+1 ):
            normConst2 = m.sqrt(m.pow(2, n)*m.factorial(n))
            normConstant = 1./(normConst1*normConst2)
            hermite_z = eval_hermite( n, zCoords )
            shapeletBasis[n, :] = normConstant*np.multiply( 
                                  hermite_z, zGauss)
            
        return shapeletBasis 

#------------------------------------------------------------------------#
    def calcLeastSquares(self, basisFunction, flatSourceArray ):
        transposeBasis = np.transpose(basisFunction)
        demoninator = np.dot(transposeBasis, basisFunction)
        numerator = np.dot(transposeBasis, flatSourceArray)
        aMoment = numerator/demoninator
        return aMoment 
        
#------------------------------------------------------------------------#
    def calcShapeletModel( self, modelName ):
        xExt = modelName.xExtent   
        yExt = modelName.yExtent
        flatModelArray = np.zeros((xExt*yExt, 1))

        totalMoments = self.moments.shape[0]
        n_max = int(-3./2 + (m.sqrt(1.+8.*totalMoments))/2)
        modelName.maxBasisIndex = n_max 
               
        uCoords = modelName.pxlCoordList[:,0]/self.major
        vCoords = modelName.pxlCoordList[:,1]/self.minor
        uShapeletBasis = self.calcShapeletBasis(uCoords, self.major, n_max)
        vShapeletBasis = self.calcShapeletBasis(vCoords, self.minor, n_max)        
        
        modelName.maxBasisIndex = n_max
            
        for k in range( totalMoments ):
            n1 = int(self.moments[k, 0])
            n2 = int(self.moments[k, 1])
            aMoment = self.moments[k, 2]
            basisFunction = uShapeletBasis[n1,:]*vShapeletBasis[n2,:]
            basisFunction = np.reshape(xExt*yExt,1)
            flatModelArray += aMoment*basisFunction
#            counter[k,0] = n1
#            counter[k,1] = n2
#            counter[k,2] = 
            
        return flatModelArray     
               
##########################################################################
## The sourceImage Class: This object contains all the pertinent info   ##
## used to calculate and compare the original dataset with the model.   ##
##                                                                      ##
##                                                                      ##
## The final two parameters in this class measure the goodness-of-fit   ## 
## to the original data used to create the moments. These are the       ##    
## normalised mean square error and the structial similarly index       ## 
## measure. The nMSE is common in astronomy but the SSIM is better      ##
## suited at accurately measuring how well an extracted image/model     ##
## "looks like" the data.                                               ##
##########################################################################
class SourceImage:
    def __init__(self, RA=0.0, Dec=0.0, extendedSource=[], yRApxl=0, 
                 xDecpxl=0, pxlCoordList=[], xAngularScale=1.0, 
                 yAngularScale=1.0, xExtent=100, yExtent=100, 
                 maxBasisIndex=5, normalisedMSE=0.0, averageSSIM=1.0, 
                 majorAxis=30, minorAxis=30, posAng=0. ): 
    
        self.RA = RA
        self.Dec = Dec
        self.extendedSource = extendedSource
        self.yRApxl = yRApxl
        self.xDecpxl = xDecpxl
        self.pxlCoordList = pxlCoordList
        self.xAngularScale = xAngularScale      
        self.yAngularScale = yAngularScale
        self.xExtent = xExtent
        self.yExtent = yExtent        
        self.maxBasisIndex = maxBasisIndex
        self.normalisedMSE = normalisedMSE
        self.averageSSIM = averageSSIM
        self.majorAxis = majorAxis
        self.minorAxis = minorAxis
        self.posAng = posAng
                
#------------------------------------------------------------------------#    
    def getFitsFile(self, filename):
        fitsFile = pyfits.open( filename )
        fitsHeader = fitsFile[0].header
        source = fitsFile[0].data
                
        # RA == y axis/column indices, Dec == x axis/row indices
        self.RA = np.radians(fitsHeader['CRVAL1'])
        self.Dec = np.radians(fitsHeader['CRVAL2'])
        self.xAngularScale = np.radians(fitsHeader['CDELT2'])
        self.yAngularScale = np.radians(fitsHeader['CDELT2'])
        self.xDecpxl = fitsHeader['CRPIX2']
        self.yRApxl = fitsHeader['CRPIX1']
        self.xExtent = fitsHeader['NAXIS2']
        self.yExtent = fitsHeader['NAXIS1']
        fitsFile.close()
        
        # Reference the extended source to the pxl (0,0)   
        self.changeRADecRef( 0, 0 ) 
        self.checkSourceData( source ) 
               
#------------------------------------------------------------------------#
    def checkSourceData(self, source):
        sizeOfData = source.shape
        dimOfData = len(sizeOfData)
        
        if (dimOfData == 2):
            self.extendedSource = source
        elif (dimOfData < 2):
            sys.exit()
        else:
             self.extendedSource = source[0,0,:,:]

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
        if (self.yRApxl != newYpxl):
            self.RA = self.RA-(self.yRApxl-newYpxl)*self.yAngularScale
            self.yRApxl = newYpxl
        if (self.xDecpxl != newXpxl):
            self.Dec = self.Dec-(self.xDecpxl-newXpxl)*self.xAngularScale
            self.xDecpxl = newXpxl

#------------------------------------------------------------------------#
    def calcCoordinateList(self, xCentrePxl, yCentrePxl):
        self.changeRADecRef( xCentrePxl, yCentrePxl )
        totalCoordinates = self.xExtent*self.yExtent
        coordinates = np.zeros((totalCoordinates, 2)) 
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
        totalCoords = self.xExtent*self.yExtent
        flatSourceArray = np.zeros((totalCoords, 1))
        
        flatSourceArray = np.flatten(self.extendedSource)

        return flatSourceArray

#------------------------------------------------------------------------#
    def expandArray( self, flatArray ):
        twoDarray = np.zeros((self.xExtent, self.yExtent))
        twoDarray = np.reshape( flatArray, (self.xExtent, self.yExtent)) 
        
        return twoDarray
        
#------------------------------------------------------------------------#
    def calcCovarianceFit( self ):
        totalFluxOfSource = np.sum( np.sum( self.extendedSource ))
        numberOfPxls = self.xExtent*self.yExtent
       
        fluxWeightedCentre = self.calcWeightedCentre(totalFluxOfSource)
        self.changeRADecRef( fluxWeightedCentre[0], fluxWeightedCentre[1] )
        
        fluxMatrix = self.calcFluxMatrix( totalFluxOfSource )
        (eigenval, eigenvect) = self.calcEigenValueDecomp( fluxMatrix )
        
        self.calcMajorMinor( eigenval )
        self.posAng = m.atan2(eigenvect[1,1],eigenvect[0,1])
        self.calcTransformCoordinates( self.posAng )

#------------------------------------------------------------------------#
   def calcMajorMinor( self, eigenval ):
   
        if eigenval[1] >= 0:
            self.majorAxis = (0.9*self.xExtent*
                                 m.pow(self.maxBasisIndex, -0.52))
        else:
            self.majorAxis = np.sqrt(eigenval[1])
            
        if eigenval[0] >= 0: 
            self.minorAxis = self.majorAxis
        else:
            self.minor = np.sqrt(eigenval[0])
       
#------------------------------------------------------------------------#       
    def calcWeightedCentre(self, totalFluxOfSource ):
        k = 0
        xCentre = 0.
        yCentre = 0.
        for i in range( self.xExtent ):
            for j in range( self.yExtent ):
                xCentre += self.extendedSource[i,j]*self.pxlCoordList[k,0]
                yCentre += self.extendedSource[i,j]*self.pxlCoordList[k,1]
                k += 1
       
        xCentre = int(xCentre/totalFluxOfSource)
        yCentre = int(yCentre/totalFluxOfSource)
        fluxWeightedCentre = ([xCentre, yCentre])
        return fluxWeightedCentre
        
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
       fluxMatrix = ([xxFlux, xyFlux, yyFlux])
       return fluxMatrix
       
#------------------------------------------------------------------------#
    def calcEigenValueDecomp(self, fluxMatrix ):
        trueFluxMatrix = np.array([[fluxMatrix[0],fluxMatrix[1]], \
        [fluxMatrix[1],fluxMatrix[2]]]) 
        
        (eigenval, eigenvect) = np.linalg.eig(trueFluxMatrix)
        return (eigenval, eigenvect)
       
#------------------------------------------------------------------------#
    def calcTransformCoordinates( self, PA ):
        transformMatrix=[[m.cos(PA), -m.sin(PA)],[m.sin(PA), m.cos(PA)]]
        coordList=np.dot(transformMatrix,np.transpose(self. pxlCoordList))
        self.pxlCoordList = np.transpose(coordList)
        
#------------------------------------------------------------------------#
    def getImageScale( self, filename, shapeletName ):
        # Get user parameters to draw new image, otherwise draw from fits
        # file
       
#        if (useFits == 1):
#            self.getFitsFile(filename)        
#        else:
#            self.xAngularScale = 1.0      
#            self.yAngularScale = 1.0
#            self.xExtent = 100
#            self.yExtent = 100
        self.getFitsFile(filename)
        self.RA = shapeletName.yRA
        self.Dec = shapeletName.xDec
           
        xCentrePxl = m.floor(self.xExtent/2)
        yCentrePxl = m.floor(self.yExtent/2)
        self.calcCoordinateList(xCentrePxl, yCentrePxl) 
        self.calcTransformCoordinates( shapeletName.positionAngle )
       
#------------------------------------------------------------------------#
    def setShapeletFromImage( self, shapeletName ):
        shapeletName.yRA = self.RA
        shapeletName.xDec = self.Dec
        shapeletName.positionAngle = self.posAng
        shapeletName.major = self.majorAxis*self.xAngularScale
        shapeletName.minor = self.minorAxis*self.yAngularScale
        
#------------------------------------------------------------------------#
    def setImageFromShapelet( self, shapeletName ):
        self.RA = shapeletName.yRA
        self.Dec = shapeletName.xDec
        self.posAng = shapeletName.positionAngle
        self.majorAxis = shapeletName.major/self.xAngularScale
        self.minorAxis = shapeletName.minor/self.yAngularScale


#    def calcModelPerformance( self, flatResidualArray ):      
        
        
        
        
        
##########################################################################       
#class ShapeletGUI:
########################## MAIN PROGRAM ##################################

filename = '../Fits_files/FnxA.fits'

source = 'FornaxA'
shapelet = 'ForA'
model = 'FnxA'

sourceName = 'FornaxA'
shapeletName = 'ForA'
modelName = 'FnxA'

shapeletName = Shapelet()
sourceName = SourceImage()
modelName = SourceImage()
residualName = SourceImage()

sourceName.getFitsFile(filename)

xCentrePxl = m.floor(sourceName.xExtent/2)
yCentrePxl = m.floor(sourceName.yExtent/2)
sourceName.calcCoordinateList(xCentrePxl, yCentrePxl)

sourceName.calcCovarianceFit()

shapeletName.calcMoments( sourceName )

sourceName.setShapeletFromImage( self, shapeletName ):
shapeletName.saveShapelet( shapelet )

#------------------------------------------------------------------------#
shapeletName.getShapelet( shapelet )

modelName.getImageScale( filename, shapeletName )
modelName.setImageFromShapelet( shapeletName )

residualName.__dict__ = sourceName.__dict__.copy()
flatModelArray = shapeletName.calcShapeletModel( modelName )
flatSourceArray = sourceName.flattenExtendedSourceArray()
fudge = np.max(flatSourceArray)/np.max(flatModelArray)
flatModelArray = fudge*flatModelArray
flatResidualArray = flatSourceArray - flatModelArray

plotSource = sourceName.expandArray( flatSourceArray )
plotModel = modelName.expandArray( flatModelArray )
plotResidual = residualName.expandArray( flatResidualArray )

plt.figure(0)
plt_max = np.max(plotSource.flatten())
plt_min = np.min(plotSource.flatten())
plt.contour(plotSource,vmin=plt_min,vmax=plt_max)
plt.colorbar()
plt.figure(1)
plt.contour(plotModel,vmin=plt_min,vmax=plt_max)
plt.colorbar()
plt.figure(2)
plt.contour(plotResidual,vmin=plt_min,vmax=plt_max)
plt.colorbar()
plt.show()


