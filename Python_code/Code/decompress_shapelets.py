#! /usr/bin/env python
# code to decompress shapelet coefficients into a FITS image
# ang_res in degrees, ang_extent in arcmin, source name undecided

# def decompress_shapelets( )
import math as m
import numpy as np
import matplotlib.pyplot as plt
import astropy

################################################################### 
# Part 1: Loading the source info. This consists of;
# parameters = [ra, dec, major, minor, pa]
# parameters = [deg, deg, arcmins, arcmins, degrees] &
# coeffs = [351,3] = (n1, n2, moment).
###################################################################
# user inputs - proper I/O later

#==============================================================================
# while True:
# 	paras = raw_input('\nPlease enter parameter file (path): ')
# 	moments = raw_input('\nPlease enter coefficient file  (path): ')
# 	resolution = float(raw_input('\nPlease enter desired angular resolution (arcmins): '))
# 	image_size = float(raw_input('\nPlease enter the total size of the reconstruction (degrees): '))
# 	tot_coeffs = int(raw_input('\nPlease enter the number of coefficients you wish to include (1-351): '))
#==============================================================================
#	break

paras = 'VirA_paras.txt'
moments = 'VirA_351coeffs.txt'
tot_coeffs = 351
image_size = 0.1
resolution = 0.1

parameters = np.loadtxt(paras)
coeffs = np.loadtxt(moments)
pxl_width = m.radians(resolution/60)		# degrees				
extent = m.radians(image_size/resolution)	# degrees

######################################################################
# Part 2: The boring bit. All numerical values are changed to radians,
# all matrices/vectors are initialised.
######################################################################
ra = m.radians(parameters[0])
dec = m.radians(parameters[1])
b1 = m.radians(parameters[2]/60)
b2 = m.radians(parameters[3]/60)
pa = m.radians(parameters[4])+m.pi

# xpxls = image size in terms of number of pixels
# num_coeffs = number of moments
# n_max = max order of hermites
xpxls = int(m.floor(extent/pxl_width))	
num_coeffs = coeffs.shape[0]		 
n_max = int(max(max(coeffs[:,0]), max(coeffs[:,1])))
	
# create coordinate axis: assumes square image
xstart = (xpxls/2+1)*pxl_width
ystart = xstart
xaxis = np.zeros((xpxls,1))
yaxis = np.zeros((xpxls,1))
tot_length = xpxls*xpxls
ra_axis = np.zeros((xpxls,1))
dec_axis = np.zeros((xpxls,1))
coords = np.zeros((tot_length,2))
sky_coords = np.zeros((tot_length,2))
mat_coords = np.zeros((tot_length,2))
xrot = np.zeros((tot_length,1))
yrot = np.zeros((tot_length,1))
all_basis = np.zeros((tot_length, n_max+1, n_max+1))
model = np.zeros((tot_length, 1))
gauss = np.zeros((tot_length, 1))
Mpiece = np.zeros((tot_length, 1))
final_image = np.zeros((xpxls, xpxls))
basis_out = np.zeros((xpxls, xpxls))

#################################################################
# Part 3: Create and transform coordinate systems. mat_coords place
# (0,0) as the centre pixel and is used for basis calculations, 
# sky_coords are the true sky coordinates in degrees. A rotation is
# performed on mat_coords to account for the source position angle.
#################################################################
k=0;
for i in range(0,xpxls):
	ra_axis[i] = m.degrees(ra + xstart - i*pxl_width)
	dec_axis[i] = m.degrees(dec - ystart + i*pxl_width)
	xaxis[i] = -xpxls/2+i;
	for j in range(0,xpxls):
		yaxis[j] = -xpxls/2+j
		mat_coords[k,0] = xaxis[i]*pxl_width
		mat_coords[k,1] = yaxis[j]*pxl_width
		sky_coords[k,0] = ra_axis[i]
		sky_coords[k,1] = dec_axis[j]
		k=k+1

transform = [[m.cos(pa), -m.sin(pa)], [m.sin(pa), m.cos(pa)]]
temp_coords = np.dot(transform, mat_coords.transpose())
rot_coords = temp_coords.transpose()
xrot = rot_coords[:,0]/b1
yrot = rot_coords[:,1]/b2

####################################################################
# Part 4: Creates the 2D hermite basis functions and scales by the
# moment. Sums the result. Hermites are the static polynomial coeffs
# for the hermite functions (much faster than running the generating
# algorithm everytime). Output is in column form so final step is to
# create an image. 
####################################################################

hermites = np.loadtxt('/home/jthomps/Dropbox/Documents/Shapelets/Common/hermite_coeffs.txt')
gauss = np.exp(-0.5*(np.array(xrot)**2+np.array(yrot)**2))	
n1 = 0                                  
n2 = 0
for n1 in range(0,n_max):
    	n2 = 0
	while (n1+n2) <= n_max:
        	norm = m.sqrt(m.pow(2,n1+n2)*m.pi*b1*b2*m.factorial(n1)*m.factorial(n2))
               	k=0
            	h1=0.
            	h2=0.
            	while (k <= n1):		
            		h1 = h1+hermites[n1, k]*(np.array(xrot))**(n1-k)
            		k=k+1
            	k=0
		while (k <= n2):		
             		h2 = h2+hermites[n2, k]*(np.array(yrot))**(n2-k)
             		k=k+1
		all_basis[:, n1, n2] = gauss/norm*h1*h2;
		n2 = n2+1

for i in range(0,tot_coeffs):
 	n1 = int(coeffs[i,0])
    	n2 = int(coeffs[i,1])
    	f_hat = coeffs[i,2]
    	Mpiece[:,0] = all_basis[:,n1,n2]
    	model = model+f_hat*(Mpiece)

k=-1
for i in range(0, xpxls):
	for j in range(0, xpxls):
		k=k+1
		final_image[i,j] = model[k]
		basis_out[i,j] = all_basis[k,0,0]

plt.imshow(final_image)
plt.show()

####################################################################
# Part 5: Generates the output as a FITS file image 
#==============================================================================
# ####################################################################
# from astropy.io import fits
# mid = round(xpxls/2);
# ang_res = m.degrees(pxl_width)
# ra = ra_axis[mid]
# dec = dec_axis[mid]
# 
# hdu = fits.PrimaryHDU(final_image)
# hdulist = fits.HDUList([hdu])
# 
# hdu.header['BSCALE'] = float(1.0)
# hdu.header['BTYPE'] = 'Intensity'
# hdu.header['BUNIT'] = 'JY/PIXEL'
# hdu.header['CRPIX1'] = mid
# hdu.header['CDELT1'] = -1*ang_res
# hdu.header['CRVAL1'] = ra[0]
# hdu.header['CTYPE1'] = 'RA---SIN'
# hdu.header['CUNIT1'] = 'deg'
# hdu.header['CRPIX2'] = mid
# hdu.header['CDELT2'] = ang_res
# hdu.header['CRVAL2'] = dec[0]
# hdu.header['CTYPE2'] = 'DEC--SIN'
# hdu.header['CUNIT2'] = 'deg'
#==============================================================================
hdu.writeto('Shapelet_output.fits')





