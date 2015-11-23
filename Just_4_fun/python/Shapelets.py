#! /usr/bin/env python
# Shapelets.py
# This code takes an image and shapelifies it!

import math as m
import numpy as np
import pylab as lab
import matplotlib.pyplot as plt
from Shape_functions import *
from Data_manipulation import *
from scipy import misc

#def Shapelets()
# Inputs
generate_new_parameters = 0             # calc paras
generate_new_moments = 1                # calc coeffs
paras = 'Starship_paras.txt'
moms = 'Starship_351coeffs.txt'

# Outputs
create_plots = 1
save_moments = 0
save_paras = 0
source_name = 'Starship'
filetype = 'txt'

##############################################################################
## Code initailisations
n_max = 30                               # maximum shapelet order
coeffs = 0.5*(n_max+1)*(n_max+2)          # number of moments
coeffs_to_inc = 351                       # manually dictate how many
fact = 5

n_approx = 5                            # for beta searcher

shapes = np.zeros((5))
face = misc.imread('../common/Starship2.jpg')
npix_rows = face.shape[0]
npix_cols = face.shape[1]
nside = max(npix_rows, npix_cols)
npix = nside*nside
midpix = round(nside/2)
predata = np.zeros((nside, nside,1))
data = np.zeros((nside, nside))
coldata = np.zeros((npix))
ang_res = m.radians(1./60)
coords = np.zeros((npix,2))
precoords = np.zeros((npix,2))

if npix_rows < nside:
	mid = round((nside-npix_rows)/2)
	rowstart = int(mid)
if npix_rows == nside:
	rowstart = 0
if npix_cols < nside:
	mid = round((nside-npix_cols)/2)
	colstart = int(mid)
if npix_cols == nside:
	colstart = 0
	 
for i in range(0,npix_rows-1):
	for j in range(0,npix_cols-1):
		gray = rbg2gray(face[i,j,:])
		data[i+rowstart,j+colstart]=gray

k=-1
for i in range(0,nside-1):
	for j in range(0,nside-1):
		k+=1
		coords[k,0] = (i-midpix)*ang_res
		coords[k,1] = (j-midpix)*ang_res
		coldata[k] = data[i,j]

obj_size = npix_rows*ang_res/20

if generate_new_parameters == 0:
	beta2 = 0.9*pow(n_max,-0.52)*npix_rows/fact*m.degrees(ang_res)*60.
	beta1 = 0.9*pow(n_max,-0.52)*npix_cols/fact*m.degrees(ang_res)*60.
	shapes = [0,0,beta1, beta2,90]
if generate_new_moments == 0: 
    moments = np.loadtxt(moms)

#######################################################################
## Setting up the variables
    
#if generate_new_parameters == 1:  
#    (major, minor, PA, offsets) = cov_fit(coords, data)
#    PA = -PA
#    shapes[0] = m.degrees(offsets[1])
#    shapes[1] = m.degrees(offsets[0])
#    shapes[4] = np.degrees(PA)

PA = np.radians(shapes[4])

#coords[:,0]-=m.radians(shapes[1])
#coords[:,1]-=m.radians(shapes[0])

transform = ([[m.cos(PA), -m.sin(PA)], [m.sin(PA), m.cos(PA)]])
rot_coords = np.dot(transform,np.transpose(coords))
im_coords = np.transpose(rot_coords)

if generate_new_parameters == 1:
	beta1, beta2, a, temp = major_minor(im_coords,coldata,n_approx,obj_size,ang_res)  
	shapes[2] = np.degrees(beta1*60)
	shapes[3] = np.degrees(beta2*60)	

print shapes   
######################################################################
## Maths time...

beta1 = m.radians(shapes[2]/60)
beta2 = m.radians(shapes[3]/60)

if generate_new_moments == 1:
	moments = deconstruct(im_coords,coldata,beta1,beta2,n_max)

col_mod = reconstruct(im_coords,moments,beta1,beta2)

#######################################################################
## Visualisations

residuals = np.zeros((nside,nside))
final_image = np.zeros((nside, nside))	
colresid = np.zeros((npix))
colresid = coldata-col_mod

k=-1
for i in range(0, nside-1):
    for j in range(0, nside-1):
        k=k+1
        residuals[i,j] = colresid[k]
        final_image[i,j] = col_mod[k]

if create_plots ==1:
	plt.figure(0)
	plt_max = np.max(data.flatten())
	plt_min = np.min(data.flatten())
	cmap=plt.get_cmap(lab.gray())
	plt.imshow(data,vmin=plt_min,vmax=plt_max)
	plt.colorbar()
	plt.figure(1)
	plt.imshow(final_image,vmin=plt_min,vmax=plt_max)
	plt.colorbar()
	plt.figure(2)
	plt.imshow(residuals,vmin=plt_min,vmax=plt_max)
	plt.colorbar()        
	plt.show()
        
new_moms = ''
new_paras = ''
start = source_name
if save_moments == 1:
    middle = ''.join(['_',str(coeffs_to_inc),'coeffs'])
    new_moms = ''.join([start,middle,'.',filetype])
    np.savetxt(new,moms, moments)
if save_paras == 1:
    new_paras = ''.join([start,'_paras.',filetype])
    np.savetxt(new_paras, shapes)

