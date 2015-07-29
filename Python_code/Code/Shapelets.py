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
filename = '../Fits_files/VirA.fits'
paras = '../Models/ForA_paras.txt'
moms = '../Models/VirA_351coeffs.txt'

# Outputs
create_plots = 1
save_moments = 0
save_paras = 0
create_fits = 0
dir_name = '../Models/'
source_name = 'VirA'
filetype = 'txt'

##############################################################################
## Code initailisations
n_max = 25                                # maximum shapelet order
coeffs = 0.5*(n_max+1)*(n_max+2)          # number of moments
coeffs_to_inc = 351                       # manually dictate how many
obj_size = 12.0

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
    
# default sizes for simple model viewing
if from_fits == 0:
    source_info[0,0] = shapes[0]
    source_info[1,0] = shapes[1]
    source_info[0,1] = ext_source/pxl_obj
    source_info[1,1] = source_info[0,1]
    source_info[0,2] = tot_pxl/2
    source_info[1,2] = source_info[0,2]
    source_info[0,3] = tot_pxl
    source_info[1,3] = source_info[0,3]

coords0,sky_coords0,dataset0,coldata0=polish_data(source_info,data,ext_source, ang_res)
    
if generate_new_parameters == 0:
    dataset = dataset0
    coords = coords0
    coldata = coldata0
    sky_coords = sky_coords0

new_info = source_info.copy()
  
if generate_new_parameters == 1:  
    (beta1, beta2, PA, offsets) = cov_fit(coords0, dataset0)
    PA = -PA
    shapes[0] = np.degrees(source_info[0,0])
    shapes[1] = np.degrees(source_info[1,0])
    shapes[4] = np.degrees(PA)
    new_info[0,0] = np.radians(shapes[0])
    new_info[1,0] = np.radians(shapes[1])
    
(coords,sky_coords,dataset,coldata)=polish_data(new_info,data,ext_source,ang_res)

PA = np.radians(shapes[4])
transform = ([[m.cos(PA), -m.sin(PA)], [m.sin(PA), m.cos(PA)]])
rot_coords = np.dot(transform,np.transpose(coords))
im_coords = np.transpose(rot_coords)

if generate_new_parameters == 1:
	beta1, beta2, a, temp = beta1_beta2(im_coords,coldata,n_approx,obj_size,ang_res)  
	shapes[2] = np.degrees(beta1*60)
	shapes[3] = np.degrees(beta2*60)	
   
######################################################################
## Maths time...

beta1 = m.radians(shapes[2]/60)
beta2 = m.radians(shapes[3]/60)

if generate_new_moments == 1:
	moments = deconstruct(im_coords,coldata,beta1,beta2,n_max)

col_mod = reconstruct(im_coords,moments,beta1,beta2)

performance = simple_stats(coldata,col_mod)
print performance

#######################################################################
## Visualisations
	
data_size = im_coords.shape
nside = int(m.sqrt(data_size[0]))
final_image=np.zeros((nside,nside))
residuals = np.zeros((nside,nside))
colresid = np.zeros((nside*nside))
colresid = coldata-col_mod

k=-1
for i in range(0, nside):
    for j in range(0, nside):
        k=k+1
        residuals[i,j] = colresid[k]
        final_image[i,j] = col_mod[k]

if create_plots ==1:
        figure(0)
        plt_max = np.max(dataset.flatten())
        plt_min = np.min(dataset.flatten())
        plt.contour(dataset,vmin=plt_min,vmax=plt_max)
        plt.colorbar()
        figure(1)
        plt.contour(final_image,vmin=plt_min,vmax=plt_max)
        plt.colorbar()
        figure(2)
        plt.contour(residuals,vmin=plt_min,vmax=plt_max)
        plt.colorbar()        
        
new_moms = ''
new_paras = ''
start = ''.join([dir_name,source_name])
if save_moments == 1:
    middle = ''.join(['_',str(coeffs_to_inc),'coeffs'])
    new_moms = ''.join([start,middle,'.',filetype])
    np.savetxt(new,moms, moments)
if save_paras == 1:
    new_paras = ''.join([start,'_paras.',filetype])
    np.savetxt(new_paras, shapes)
if generate_fits == 1:
    make_fitsfile(new_info, model) 

