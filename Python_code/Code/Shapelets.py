#! /usr/bin/env python
# Shapelets.py
# This code will either take a fits file and compute/load parameters, moments, create a model and evalute the performance. Alternatively, it can be used to simply view a shapelet model.

import math as m
import numpy as np
import matplotlib.pyplot as plt
from Shape_functions import *
from Data_manipulation import *

#def Shapelets()
# Inputs
from_fits = 1                           # use a fits file
generate_new_parameters = 1             # calc paras
generate_new_moments = 1                # calc coeffs
filename = '../Fits_files/FnxA.fits'
paras = '../Models/ForA_paras.txt'
moms = '../Models/ForA_351coeffs.txt'

# Outputs
create_plots = 1                        # plot the results
save_moments = 0                        # save the results
save_paras = 0                          # save the parameters
min_coeffs = 1                          # minimise the number of coefficients used
cutoff_nmse = 10e-6
dir_name = '../Models/'
source_name = 'ForA'
filetype = 'txt'

##############################################################################
## Code initailisations
n_max = 25                                # maximum shapelet order
coeffs = 0.5*(n_max+1)*(n_max+2)          # number of moments
coeffs_to_inc = 351                       # manually dictate how many
obj_size = 8.0                            # in arcmin

n_approx = 5                            # for beta searcher

shapes = np.zeros((5))

if from_fits == 1:
    (source_info, data) = import_fits(filename)

if generate_new_parameters == 0:
    shapes = np.loadtxt(paras)
if generate_new_moments == 0: 
    moments = np.loadtxt(moms)


ext_source = np.radians(obj_size/60)
ang_res = max(abs(source_info[0,1]), abs(source_info[1,1]))  

#######################################################################
## Setting up the variables
    
coords0,sky_coords0,dataset0,coldata0=polish_data(source_info,data,ext_source, ang_res)

if generate_new_parameters == 0:
    dataset = dataset0
    coords = coords0
    coldata = coldata0
    sky_coords = sky_coords0

new_info = source_info.copy()
  
if generate_new_parameters == 1:  
    (major, minor, PA, offsets) = cov_fit(coords0, dataset0)
    PA = -PA
    shapes[0] = np.degrees(source_info[0,0])+m.degrees(offsets[1])
    shapes[1] = np.degrees(source_info[1,0])+m.degrees(offsets[0])
    shapes[4] = np.degrees(PA)
    new_info[0,0] = np.radians(shapes[0])
    new_info[1,0] = np.radians(shapes[1])
    new_info[0,2] += round(offsets[1]/ang_res)
    new_info[1,2] += round(offsets[0]/ang_res)
#np.savetxt('/home/jthomps/Dropbox/Documents/Shapelets/Image_shapes/Text/Vir_coords.txt',coords0)

(coords,sky_coords,dataset,coldata)=polish_data(new_info,data,ext_source,ang_res)

PA = np.radians(shapes[4])
transform = ([[m.cos(PA), -m.sin(PA)], [m.sin(PA), m.cos(PA)]])
rot_coords = np.dot(transform,np.transpose(coords))
im_coords = np.transpose(rot_coords)

if generate_new_parameters == 1:
	beta1, beta2, a, temp = major_minor(im_coords,coldata,n_approx,ext_source,ang_res)  
	shapes[2] = np.degrees(beta1*60)
	shapes[3] = np.degrees(beta2*60)	

######################################################################
## Maths time...

data_size = im_coords.shape
nside = int(m.sqrt(data_size[0]))
colresid = np.zeros((nside*nside))
final_image=np.zeros((nside,nside))
residuals = np.zeros((nside,nside))

beta1 = m.radians(shapes[2]/60)
beta2 = m.radians(shapes[3]/60)

if generate_new_moments == 1:
	moments = deconstruct(im_coords,coldata,beta1,beta2,n_max)

moments = moments[np.argsort(abs(moments[:,2]))]
moments = np.flipud(moments)
if min_coeffs == 1:
    tot_coeffs = minco(im_coords, coldata, moments, shapes[2], shapes[3], cutoff_nmse)
else:
    tot_coeffs = coeffs_to_inc

moms_inc = np.zeros((tot_coeffs, 3))
for i in range(0, tot_coeffs+1):
    moms_inc[i,:] = moments[i,:]

col_mod = reconstruct(im_coords,moms_inc,beta1,beta2)
colresid = coldata-col_mod

k=-1
for i in range(0, nside):
    for j in range(0, nside):
        k=k+1
        residuals[i,j] = colresid[k]
        final_image[i,j] = col_mod[k]

performance = simple_stats(dataset, final_image, residuals)
print performance

#######################################################################
## Visualisations
	
if create_plots ==1:
        plt.figure(0)
        plt_max = np.max(dataset.flatten())
        plt_min = np.min(dataset.flatten())
        plt.contour(dataset,vmin=plt_min,vmax=plt_max)
        plt.colorbar()
        plt.figure(1)
        plt.contour(final_image,vmin=plt_min,vmax=plt_max)
        plt.colorbar()
        plt.figure(2)
        plt.contour(residuals,vmin=plt_min,vmax=plt_max)
        plt.colorbar()        
	plt.show()
        
new_moms = ''
new_paras = ''
start = ''.join([dir_name,source_name])
if save_moments == 1:
    middle = ''.join(['_',str(coeffs_to_inc),'coeffs'])
    new_moms = ''.join([start,middle,'.',filetype])
    np.savetxt(new_moms, moms_inc)
if save_paras == 1:
    new_paras = ''.join([start,'_paras.',filetype])
    np.savetxt(new_paras, shapes)
if create_fits == 1:
    make_fitsfile(new_info, model) 

