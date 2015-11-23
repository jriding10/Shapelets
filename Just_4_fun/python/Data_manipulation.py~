#! /usr/bin/env python
# Data_manipulation.py
# Contains a number of functions used to prepare the data for use by
# the main maths routinues. 
# def userIO        :inputs user choices
# def import_fits    :reads fits data and header
# def polish_data    :centers source, creates coordinates
# def cov_fit        :determines PA and flux center
# def make_fitsfile :creates a fits file containing the model image

import math as m
import numpy as np
import pyfits

#######################################################################
## Does a covariance fit of the data to determine flux center

def cov_fit(coords, data):
    S = sum(sum(data))
    npix_row = data.shape[0]
    npix_col = data.shape[1]
    npix = npix_row*npix_col
    offsets = (0.0,0.0)
    addmat = np.zeros((npix,2))

    x = 0
    y = 0
    Sxx = 0
    Sxy = 0
    Syy = 0

    k=-1
    for i in range(0,npix_row):
        for j in range(0,npix_col):
            k += 1
            x += data[i,j]*coords[k,0]
            y += data[i,j]*coords[k,1]

    x0 = (x/S)
    y0 = (y/S)

    offsets = (x0, y0)

    coords[:,0]-= x0
    coords[:,1]-= y0

    k=-1
    for i in range(0,npix_row):
        for j in range(0,npix_col):
            k+=1
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

######################################################################
def rbg2gray(rgb):

    r, g, b = rgb[0], rgb[1], rgb[2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray  

