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
import matplotlib
import astropy.io.fits as pyfits
from scipy.special import eval_hermite
from tkinter import *
from tkinter import filedialog 
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.figure import Figure


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
                position_angle=90):

        self.yRA = yRA
        self.xDec = xDec
        self.moments = moments
        self.major = major
        self.minor = minor        
        self.position_angle = position_angle
        
#------------------------------------------------------------------------#
    def saveShapelet(self, source_name):
        filename = '../Models/' + source_name + '_moments.txt'
        np.savetxt(filename, self.moments)
        
        shapelet_parameters = ([self.yRA, self.xDec, self.major,
                    self.minor, self.position_angle])
                    
        filename = '../Models/' + source_name + '_parameters.txt'
        np.savetxt(filename, shapelet_parameters)

#------------------------------------------------------------------------#  
    def getShapelet(self, source_name):
        filename = '../Models/' + source_name + '_moments.txt'
        self.moments = np.loadtxt(filename)
        
        filename = '../Models/' + source_name + '_parameters.txt'
        shapelet_parameters = np.loadtxt(filename)

        self.yRA = shapelet_parameters[0]
        self.xDec = shapelet_parameters[1]
        self.major = shapelet_parameters[2]
        self.minor = shapelet_parameters[3]
        self.position_angle = shapelet_parameters[4]
        
#------------------------------------------------------------------------# 
    def calcMoments( self, source ):
        u_coords = source.pxl_coords_list[:,0]/self.major
        v_coords = source.pxl_coords_list[:,1]/self.minor
        
        n_max = source.max_basis_index
        u_basis = self.calcShapeletBasis(u_coords, self.major, n_max)
        v_basis = self.calcShapeletBasis(v_coords, self.minor, n_max)
        flattened_source = source.flattenExtendedSource()

        total_moments = int(0.5*(n_max+1)*(n_max+2))
        self.moments = np.zeros((total_moments, 3))
        
        k=0
        for n1 in range(n_max+1 ):
            for n2 in range( n_max+1 ):
                if (n1+n2) <= n_max:
                    basis_function=u_basis[n1,:]*v_basis[n2,:]
                    a_moment = (self.calcLeastSquares(basis_function,
                                          flattened_source))
                    self.moments[k,0] = n1
                    self.moments[k,1] = n2 
                    self.moments[k,2] = a_moment
                    k += 1
                    
#------------------------------------------------------------------------#
    def calcShapeletBasis( self, z_coords, beta, n_max ):
        total_coords = z_coords.shape[0]
        z_gauss = np.exp(0.5*np.power(z_coords,2))
        norm_const1 = m.sqrt(beta*m.sqrt(m.pi))
        shapelet_basis = np.zeros(( n_max+1, total_coords ))
                       
        for n in range( n_max+1 ):
            norm_const2 = m.sqrt(m.pow(2, n)*m.factorial(n))
            norm_constant = 1./(norm_const1*norm_const2)
            hermite_z = eval_hermite( n, z_coords )
            shapelet_basis[n, :] = norm_constant*np.multiply(
                                  hermite_z, z_gauss)
            
        return shapelet_basis

#------------------------------------------------------------------------#
    def calcLeastSquares(self, basis_function, flattened_source ):
        transpose_basis = np.transpose(basis_function)
        demoninator = np.dot(transpose_basis, basis_function)
        numerator = np.dot(transpose_basis, flattened_source)
        a_moment = numerator/demoninator
        return a_moment
        
#------------------------------------------------------------------------#
    def calcShapeletModel( self, model ):
        xExt = model.xextent
        yExt = model.yextent
        flattened_model = np.zeros((xExt*yExt, 1))

        total_moments = self.moments.shape[0]
        n_max = int(-3./2 + (m.sqrt(1.+8.*total_moments))/2)
        model.max_basis_index = n_max
               
        u_coords = model.pxl_coords_list[:,0]/self.major
        v_coords = model.pxl_coords_list[:,1]/self.minor
        u_basis = self.calcShapeletBasis(u_coords, self.major, n_max)
        v_basis = self.calcShapeletBasis(v_coords, self.minor, n_max)
        
        model.max_basis_index = n_max
            
        for k in range( total_moments ):
            n1 = int(self.moments[k, 0])
            n2 = int(self.moments[k, 1])
            a_moment = self.moments[k, 2]
            basis_function = u_basis[n1,:]*v_basis[n2,:]
            basis_function = np.reshape(xExt*yExt,1)
            flattened_model += a_moment*basis_function

        return flattened_model
               
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
    def __init__(self, RA=0.0, Dec=0.0, extended_source=None, yRApxl=0,
                 xDecpxl=0, pxl_coords_list=None, x_angscale=1.0,
                 y_angscale=1.0, xextent=100, yextent=100,
                 max_basis_index=5, normalised_mse=0.0, avg_ssim=1.0,
                 major_axis=30, minor_axis=30, posAng=0. ):
    
        self.RA = RA
        self.Dec = Dec
        self.extended_source = extended_source
        self.yRApxl = yRApxl
        self.xDecpxl = xDecpxl
        self.pxl_coords_list = pxl_coords_list
        self.x_angscale = x_angscale
        self.y_angscale = y_angscale
        self.xextent = xextent
        self.yextent = yextent
        self.max_basis_index = max_basis_index
        self.normalised_mse = normalised_mse
        self.avg_ssim = avg_ssim
        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.posAng = posAng
                
#------------------------------------------------------------------------#    
    def getFitsFile(self, filename):
        fitsFile = pyfits.open( filename )
        fitsHeader = fitsFile[0].header
        source = fitsFile[0].data
                
        # RA == y axis/column indices, Dec == x axis/row indices
        self.RA = np.radians(fitsHeader['CRVAL1'])
        self.Dec = np.radians(fitsHeader['CRVAL2'])
        self.x_angscale = np.radians(fitsHeader['CDELT2'])
        self.y_angscale = np.radians(fitsHeader['CDELT2'])
        self.xDecpxl = fitsHeader['CRPIX2']
        self.yRApxl = fitsHeader['CRPIX1']
        self.xextent = fitsHeader['NAXIS2']
        self.yextent = fitsHeader['NAXIS1']
        fitsFile.close()
        
        # Reference the extended source to the pxl (0,0)   
        self.changeRADecRef( 0, 0 ) 
        self.checkSourceData( source ) 
               
#------------------------------------------------------------------------#
    def checkSourceData(self, source):
        size_of_data = source.shape
        dim_of_data = len(size_of_data)
        
        if (dim_of_data == 2):
            self.extended_source = source
        elif (dim_of_data < 2):
            sys.exit()
        else:
             self.extended_source = source[0,0,:,:]

#------------------------------------------------------------------------#
    def resizeImageOfSource(self, new_xstart, new_ystart, new_xend, new_yend):
        self.xextent = new_xend-new_xstart
        self.yextent = new_yend-new_ystart
        new_source = np.zeros((xextent, yextent))

        for i in range(self.xextent):
            for j in range(self.yextent):
                x = int(i+new_xstart)
                y = int(j+new_ystart)
                new_source[i, j] = self.extended_source[x, y]
         
        self.changeRADecRef(new_xstart, new_ystart)
        self.yRApxl = 0
        self.xDecpxl = 0
        self.extended_source = new_source

#------------------------------------------------------------------------#       
    def changeRADecRef( self, new_xpxl, new_ypxl ):
        if (self.yRApxl != new_ypxl):
            self.RA = self.RA-(self.yRApxl-new_ypxl)*self.y_angscale
            self.yRApxl = new_ypxl
        if (self.xDecpxl != new_xpxl):
            self.Dec = self.Dec-(self.xDecpxl-new_xpxl)*self.x_angscale
            self.xDecpxl = new_xpxl

#------------------------------------------------------------------------#
    def calcCoordinateList(self, xcentre_pxl, ycentre_pxl):
        self.changeRADecRef( xcentre_pxl, ycentre_pxl )
        total_coordinates = self.xextent*self.yextent
        coordinates = np.zeros((total_coordinates, 2))
        xstart = -1*xcentre_pxl
        ystart = -1*ycentre_pxl
        
        k = 0
        for i in range( self.xextent ):
            for j in range( self.yextent ):
                 coordinates[k, 0] = xstart + i
                 coordinates[k, 1] = ystart + j
                 k += 1
        self.pxl_coords_list = coordinates

#------------------------------------------------------------------------#
    def flattenExtendedSource(self):

        ext_source = self.extended_source
        flattened_source = ext_source.flatten()

        return flattened_source

#------------------------------------------------------------------------#
    def expandArray( self, flat_array ):
        twoDarray = np.zeros((self.xextent, self.yextent))
        twoDarray = np.reshape( flat_array, (self.xextent, self.yextent))
        
        return twoDarray
        
#------------------------------------------------------------------------#
    def calcCovarianceFit( self ):
        sum_of_source = np.sum( np.sum( self.extended_source ))

        weighted_centre = self.calcWeightedCentre(sum_of_source)
        self.changeRADecRef( weighted_centre[0], weighted_centre[1] )
        
        covar_matrix = self.calcFluxMatrix( sum_of_source )
        (eigenval, eigenvect) = self.calcEigenValueDecomp( covar_matrix )
        
        self.calcMajorMinor( eigenval )
        self.posAng = m.atan2(eigenvect[1,1],eigenvect[0,1])
        self.calcTransformCoordinates( self.posAng )

#------------------------------------------------------------------------#
    def calcMajorMinor( self, eigenval ):
   
        if eigenval[1] < 0:
            self.major_axis = (0.9*self.xextent*
                                 m.pow(self.max_basis_index, -0.52))
        else:
            self.major_axis = np.sqrt(eigenval[1])
            
        if eigenval[0] < 0:
            self.minor_axis = self.major_axis
        else:
            self.minor_axis = np.sqrt(eigenval[0])
       
#------------------------------------------------------------------------#       
    def calcWeightedCentre(self, sum_of_source ):
        k = 0
        xcentre = 0.
        ycentre = 0.
        for i in range(self.xextent):
            for j in range(self.yextent):
                xcentre += self.extended_source[i,j]*self.pxl_coords_list[k,0]
                ycentre += self.extended_source[i,j]*self.pxl_coords_list[k,1]
                k += 1
       
        xcentre = int(xcentre/sum_of_source)
        ycentre = int(ycentre/sum_of_source)
        weighted_centre = ([xcentre, ycentre])
        return weighted_centre
        
#------------------------------------------------------------------------#
    def calcFluxMatrix(self, sum_of_source):
       # Matrix terms
       xx_value = 0.
       yy_value = 0.
       xy_value = 0.
       
       k = 0
       for i in range( self.xextent ):
           for j in range( self.yextent ):
               xx_value = xx_value + self.extended_source[i,j]* \
               self.pxl_coords_list[k, 0]*self.pxl_coords_list[k,0]
               yy_value = xx_value + self.extended_source[i,j]* \
               self.pxl_coords_list[k, 1]*self.pxl_coords_list[k,1]
               xy_value = xx_value + self.extended_source[i,j]* \
               self.pxl_coords_list[k, 0]*self.pxl_coords_list[k,1]
               k += 1
        
       xx_value = xx_value/sum_of_source
       yy_value = yy_value/sum_of_source
       xy_value = xy_value/sum_of_source
       covar_matrix = ([xx_value, xy_value, yy_value])
       return covar_matrix
       
#------------------------------------------------------------------------#
    def calcEigenValueDecomp(self, covar_matrix ):
        eigenmatrix = np.array([[covar_matrix[0],covar_matrix[1]],
        [covar_matrix[1],covar_matrix[2]]])
        
        (eigenval, eigenvect) = np.linalg.eig(eigenmatrix)
        return (eigenval, eigenvect)
       
#------------------------------------------------------------------------#
    def calcTransformCoordinates( self, PA ):
        transform_matrix=[[m.cos(PA), -m.sin(PA)],[m.sin(PA), m.cos(PA)]]
        coord_list=np.dot(transform_matrix,np.transpose(self. pxl_coords_list))
        self.pxl_coords_list = np.transpose(coord_list)
        
#------------------------------------------------------------------------#
    def getImageScale( self, filename, shapelet ):
        # Get user parameters to draw new image, otherwise draw from fits
        # file
       
#        if (useFits == 1):
#            self.getFitsFile(filename)        
#        else:
#            self.x_angscale = 1.0
#            self.y_angscale = 1.0
#            self.xextent = 100
#            self.yextent = 100
        self.getFitsFile(filename)
        self.RA = shapelet.yRA
        self.Dec = shapelet.xDec
           
        xcentre_pxl = m.floor(self.xextent/2)
        ycentre_pxl = m.floor(self.yextent/2)
        self.calcCoordinateList(xcentre_pxl, ycentre_pxl)
        self.calcTransformCoordinates( shapelet.position_angle )
       
#------------------------------------------------------------------------#
    def setShapeletFromImage( self, shapelet ):
        shapelet.yRA = self.RA
        shapelet.xDec = self.Dec
        shapelet.position_angle = self.posAng
        shapelet.major = self.major_axis*self.x_angscale
        shapelet.minor = self.minor_axis*self.y_angscale
        
#------------------------------------------------------------------------#
    def setImageFromShapelet( self, shapelet ):
        self.RA = shapelet.yRA
        self.Dec = shapelet.xDec
        self.posAng = shapelet.position_angle
        self.major_axis = shapelet.major/self.x_angscale
        self.minor_axis = shapelet.minor/self.y_angscale


#    def calcModelPerformance( self, flatResidualArray ):      
        
        
        
        
        
##########################################################################       
class ShapeletGUI(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.parent.title('Shapelets')
        self.fileLoad()
        self.ImageDisplay()
        self.SelectArea()
        self.CalcShapelet()
        self.SaveResult()
        self.Quitter()

#------------------------------------------------------------------------#
    def quitter():
        global root
        if askyesno('Verify', 'Do you really want to quit?'):
            root.destroy()

#------------------------------------------------------------------------#
    def fileLoad():
        filename = filedialog.askopenfilename(initialdir = "/",
                                              title = "Select file",
                                              filetypes =
                                              ("fits files","*.fits"))
        source.getFitsFile(filename)

#------------------------------------------------------------------------#
    def displayFits():
        fig = Figure(figsize(6,6))
        fits_image = fig.add_subplot(111)

        fits_data = source.extended_source
        plt_max = np.max(fits_data.flatten())
        plt_min = np.min(fits_data.flatten())
        contour(fits_data,vmin=plt_min,vmax=plt_max)
        colorbar()

        fits_image.set_title('Raw Fits Data')
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.get_tk_widget().pack()
        canvas.draw()
        
########################## MAIN PROGRAM ##################################

root = Tk()#------------------------------------------------------------------------#
shapelet = Shapelet()
source = SourceImage()
model = SourceImage()
residual = SourceImage()

app = ShapeletGUI(root)
root.mainloop()

source_name = 'FornaxA'

shapelet = Shapelet()
source = SourceImage()
model = SourceImage()
residual = SourceImage()


x_centre_pxl = m.floor(source.xextent/2)
y_centre_pxl = m.floor(source.yextent/2)
source.calcCoordinateList(x_centre_pxl, y_centre_pxl)

source.calcCovarianceFit()

shapelet.calcMoments(source)

source.setShapeletFromImage(shapelet)
shapelet.saveShapelet(source_name)

#------------------------------------------------------------------------#
shapelet.getShapelet(source_name)

model.getImageScale(filename, shapelet)
model.setImageFromShapelet(shapelet)

residual.__dict__ = source.__dict__.copy()
flattened_model = shapelet.calcShapeletModel( model )
flattened_source = source.flattenExtendedSource()
fudge = np.max(flattened_source)/np.max(flattened_model)
flattened_model = fudge*flattened_model
flatResidualArray = flattened_source - flattened_model

plotSource = source.expandArray( flattened_source )
plotModel = model.expandArray( flattened_model )
plotResidual = residual.expandArray( flatResidualArray )

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


