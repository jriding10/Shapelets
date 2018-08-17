#! /usr/bin/env python
# Data_manipulation.py
# Contains a number of functions used to prepare the data for use by
# the main maths routinues. 
# def import_fits    :reads fits data and header
# def resize_data    :centers source, creates coordinates
# def cov_fit        :determines PA and flux center
# def make_fitsfile :creates a fits file containing the model image

import math as m
import numpy as np
from astropy.io import fits

##########################################################################
## import_fits: takes the specified fits filename, opens it and extracts
# the relevant information. Note that degrees are converted into radians.
#========================================================================# 

def import_fits(filename):
    # open file
    obj_info = np.zeros((2,4))
    obj_data = fits.open(filename)
    fits_info = obj_data[0].header
    alldata = obj_data[0].data

    # relevant coordinates from header
    ra = np.radians(fits_info['CRVAL1'])
    dec = np.radians(fits_info['CRVAL2'])
    ra_res = np.radians(fits_info['CDELT1'])
    dec_res = np.radians(fits_info['CDELT2'])
    ra_pxl = fits_info['CRPIX1']
    dec_pxl = fits_info['CRPIX2']

    ra_side = fits_info['NAXIS1']
    dec_side = fits_info['NAXIS2']

    data = alldata[0,0,:,:]

    # packaging row 0 = ra stuff, row 1 = dec stuff
    obj_info[0,:] = [ra, ra_res, ra_pxl, ra_side]
    obj_info[1,:] = [dec, dec_res, dec_pxl, dec_side]
   
    obj_data.close() 
    return (obj_info, data)

#######################################################################
## resize_data: responsible for resizing the data, transforming it into
# columns and creating meaningful axii.
#=====================================================================#

def resize_data(sinfo, in_data, new_size):

    fudge = 100.
    # ra corresponds to columns, dec to rows in matrix format
    num_rows = sinfo[1,3]-1
    num_cols = sinfo[0,3]-1
    ra_pxl = sinfo[0,2]
    dec_pxl = sinfo[1,2]

    new_num_rows = new_size[1,3]-1
    new_num_cols = new_size[0,3]-1
    new_ra_pxl = new_size[0,2]
    new_dec_pxl = new_size[1,2]

    # checking that the new size is not greater than the original
    if (sinfo[0,2]+nside/2) > side:
        nside = 2*(side-midpix_ra)
    if (sinfo[0,2]-nside/2) < 0:
        nside = 2*midpix_ra

    if (sinfo[1,2]+nside/2) > side:
        nside = 2*(side-midpix_dec)
    if (sinfo[1,2]-nside/2) < 0:
        nside = 2*midpix_dec
    
     # Initialise arrays
    nside = int(nside)
    midpix = int(nside/2)
    num_pxl = int(nside*nside)

    coords = np.zeros((num_pxl,2))
    rowaxis = np.zeros((nside,1))
    ra = np.zeros((nside,1))
    dec = np.zeros((nside,1))
    sky_coords = np.zeros((nside,2))
    out_data = np.zeros((nside, nside))
    col_data = np.zeros((num_pxl,1))

    for i in range(0,nside):
        rowaxis[i] = (i-midpix+1)*ang_res

    off1 = sinfo[1,2]-midpix
    off0 = sinfo[0,2]-midpix

    k=-1
    for i in range(0,nside):
        for j in range(0,nside):     
              k+=1
              coords[k,0]=rowaxis[i]
              coords[k,1]=rowaxis[j]
              out_data[i,j] = in_data[i+off1, j+off0]
              col_data[k,0] = out_data[i,j]  

    dec = np.degrees(rowaxis) + sinfo[1,0]*np.ones((nside,1))
    ra = -1*np.degrees(rowaxis) - sinfo[0,0]*np.ones((nside,1))
    sky_coords = np.concatenate((ra, dec), axis=1)

    return coords, sky_coords, out_data, col_data
    

#######################################################################
## Does a covariance fit of the data to determine flux center

def cov_fit(coords, data):
    S = sum(sum(data))
    nside = m.sqrt(coords.shape[0])
    nside = int(nside)
    npix = nside*nside
    offsets = (0.0,0.0)
    addmat = np.zeros((npix,2))

    x = 0
    y = 0
    Sxx = 0
    Sxy = 0
    Syy = 0

    k=-1
    for i in range(0,nside):
        for j in range(0,nside):
            k += 1
            x += data[i,j]*coords[k,0]
            y += data[i,j]*coords[k,1]

    x0 = (x/S)
    y0 = (y/S)

    offsets = (x0, y0)

    coords[:,0]-= x0
    coords[:,1]-= y0

    k=-1
    for i in range(0,nside):
        for j in range(0,nside):
            k=k+1
            Sxx = Sxx + data[i,j]*coords[k,0]*coords[k,0]
            Sxy = Sxy + data[i,j]*coords[k,0]*coords[k,1]
            Syy = Syy + data[i,j]*coords[k,1]*coords[k,1]

    a11 = Sxx/S
    a12 = Sxy/S
    a22 = Syy/S

    C = np.array([[a11, a12], [a12, a22]])

    (eigenval, eigenvect) = np.linalg.eig(C)

    minor = np.sqrt(eigenval[0])
    major = np.sqrt(eigenval[1])

    PA = m.atan2(eigenvect[1,1],eigenvect[0,1])

    return major, minor, PA, offsets

#######################################################################
## Makes a fits file of the model

def make_fitsfile(sky_coords, model):

    data_size = model.shape
    nside = data_size[0]
    mid = round(nside/2)
    ra = ra_axis[mid]
    dec = dec_axis[mid]
    ang_res = abs(sky_coords[0,1]-sky_coords[0,0])
    ang_res = m.degrees(ang_res)    

    hdu = fits.PrimaryHDU(final_image)
    hdulist = fits.HDUList([hdu])

    hdu.header['BSCALE'] = float(1.0)
    hdu.header['BTYPE'] = 'Intensity'
    hdu.header['BUNIT'] = 'JY/PIXEL'
    hdu.header['CRPIX1'] = mid
    hdu.header['CDELT1'] = -1*ang_res
    hdu.header['CRVAL1'] = ra[0]
    hdu.header['CTYPE1'] = 'RA---SIN'
    hdu.header['CUNIT1'] = 'deg'
    hdu.header['CRPIX2'] = mid
    hdu.header['CDELT2'] = ang_res
    hdu.header['CRVAL2'] = dec[0]
    hdu.header['CTYPE2'] = 'DEC--SIN'
    hdu.header['CUNIT2'] = 'deg'
    hdu.writeto('Shapelet_output.fits')

    return 
