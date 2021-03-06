#! /usr/bin/env python
# Useful.py

import io
import math as m
import numpy as np
import astropy
import matplotlib.pyplot as plt

#######################################################################
##  To make squares out of columns

def make_image(coldata):
	
	data_size = coldata.shape
	nside = int(m.sqrt(data_size[0]))
	new_image = np.zeros((nside,nside))

	k=-1
	for i in range(0,nside):
		for j in range(0,nside):
			k+=1
			new_image[i,j]=coldata[k]

	return new_image

#######################################################################
## Display coefficients

def view_coeffs(moments):

    data_size = moments.shape
    npix = int(data_size[0])
    none = max(moments[:,0])
    ntwo = max(moments[:,1])
    nmax = max(none, ntwo)
    
    coeffs = np.zeros((nmax+1,nmax+1))
    
    for i in range(0,npix):
        n1 = moments[i,0]
        n2 = moments[i,1]
        coeffs[n1, n2]=moments[i,2]
        
    return coeffs
    

 
    


