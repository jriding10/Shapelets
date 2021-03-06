#! /usr/bin/env python
# Shape_functions.py
# Contains a number of functions used in the shapelet main code

import math as m
import numpy as np

#######################################################################
## Determines the shapelet moments

def deconstruct(coords, coldata, beta1, beta2, nmax):

    min_percent = 0.0    # zero small coefficients
    tot_moms = 0.5*(nmax+1)*(nmax+2)
    f_coeffs = np.zeros((tot_moms,3))
    npix = coords.shape[0]
    Mpiece = np.zeros((npix,1))
    M = basis(coords, beta1, beta2, nmax)

    i=-1
    for n1 in range(0,nmax):
        for n2 in range(0,nmax):
            if (n1+n2)<=nmax:            
                Mpiece[:,0] = M[:,n1,n2]
                demon = np.dot(np.transpose(Mpiece),Mpiece)
                numer = np.dot(np.transpose(Mpiece),coldata)
                fhat = numer/demon
                if (n1+n2) == 0:
                    fmax = fhat
                if abs(fhat) > min_percent*fmax:
                    i+=1
                    f_coeffs[i,0]=n1
                    f_coeffs[i,1]=n2
                    f_coeffs[i,2]=fhat
                else:
                    i+=1
                    f_coeffs[i,0]=n1
                    f_coeffs[i,1]=n2
                    f_coeffs[i,2]=0.0    
    
    return f_coeffs

#######################################################################
## Creates the orthogonal basis functions

def basis(coords,b1,b2,nmax):

    hermites = np.loadtxt('../common/hermite_coeffs.txt')
    xrot = coords[:,0]/b1
    yrot = coords[:,1]/b2
    npix = coords.shape[0]
    all_basis = np.zeros((npix,nmax+1,nmax+1))    

    gauss = np.exp(-0.5*(np.array(xrot)**2+np.array(yrot)**2))    
    n1 = 0                                  
    n2 = 0
    for n1 in range(0,nmax):
        n2 = 0
        while (n1+n2) <= nmax:
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

    return all_basis

#######################################################################
## Calculates some simple statistics about the data - NMSE & PSNR

def simple_stats(coldata, colmod):

    data_size = coldata.shape
    npix = data_size[0]
    performance = np.zeros((2))
    norm1 = 0
    norm2 = 0
    mse = 0
    PSNR = 0

    for i in range(0,npix-1):
           z = coldata[i] - colmod[i]
           norm1+=coldata[i]
           norm2+=colmod[i]
           mse+=z*z

    NMSE = mse/(norm1*norm2)

    G = np.max(coldata)

    PSNR = 20*np.log10(G)+20*np.log10(npix)-10*np.log10(mse)
    performance[0] = NMSE
    performance[1] = PSNR

    return performance

#######################################################################
## Handles the function calls of beta1_beta2

def quick_model(coords,coldata,beta1,beta2,n):

    npix = coords.shape[0]
    tot_moms = 0.5*(n-1)*(n-2)        
    temp_moments = np.zeros((tot_moms,3))
    temp_mod = np.zeros((npix,1))

    temp_moments = deconstruct(coords,coldata,beta1,beta2,n)
    temp_mod = reconstruct(coords,temp_moments,beta1,beta2)
    performance = simple_stats(coldata,temp_mod)

    return performance[0]

#######################################################################
## reconstructs the model

def reconstruct(coords,moments,beta1,beta2):

    npix = coords.shape[0]
    tot_coeffs = moments.shape[0]
    n = int(max(moments[:,0]))
    Mpiece = np.zeros((npix))
    colmod = np.zeros((npix))

    M = basis(coords,beta1,beta2,n)
    
    for i in range(0,tot_coeffs):
        n1 = int(moments[i,0])
        n2 = int(moments[i,1])
        f_hat = moments[i,2]
        Mpiece[:] = M[:,n1,n2]
        colmod += f_hat*(Mpiece)
    
    return colmod

#######################################################################
## Determines how many coefficients are needed to match the mse of the noise around the source.

def minco(coords,coldata,moments,beta1,beta2,mse):

    npix = coldata.shape[0]
    Mpiece = np.zeros((npix,1))
    factor = 10
    colmod = np.zeros((npix,1))
    
    none = max(moments[:,0])
    ntwo = max(moments[:,1])
    n_max = int(max(none, ntwo))

    M = basis(coords,beta1,beta2,n_max)

    resid_nmse = 100
    target=mse/factor
    i=0

    while (resid_nmse > target) and (i<350):
        i+=1
        n1 = moments[i,0]
        n2 = moments[i,1]
        fhat = moments[i,2]
        Mpiece[:,0] =M[:,n1,n2]
        colmod += fhat*(Mpiece)

        z=0
        for k in range(0,npix-1):
           z = coldata[k] - colmod[k]
           mse+=z*z

        resid_nmse = mse/npix
        print i, resid_nmse

    tcoeffs = i+1
    print tcoeffs
    return tcoeffs

#######################################################################

def major_minor(coords,coldata,n,ext_source,ang_res):
    
    nside = m.sqrt(coldata.shape[0])
    max_iters = 2000
    fine = 10e6        # precision of the two stages
    coarse = 10e4
    temp = np.ones((max_iters,3))
    done = 0 
    second = 0
    imsig = 100

    max_beta = 0.9*pow(n,-0.52)*nside*ang_res
    min_beta = ang_res/2
    guess = 0.9*pow(n,-0.52)*ext_source
                          
    # Initialise the loop
    prev_beta1 = guess  
    prev_beta2 = guess  
    MSE = imsig                                
    prev_MSE = 2*imsig  
    beta2 = prev_beta2 
    change = 0                                 
    c1=0                                       
    c2=0                         
    step = 5*ang_res 
    precise = coarse
    print max_beta, min_beta, guess

    l=-1
    # Coarse loop to isolate beta 1 and 2: Minimise beta1, then beta2. If beta1 and beta2 are the same as last minimisation loop, convergence is assumed.
    while (done == 0 & l<max_iters-1):
        if c1 != 1:
            beta1 = prev_beta1+5*step
            if beta1 > max_beta:
                beta1 = max_beta
            while (MSE < prev_MSE):
                l+=1 
                prev_MSE = MSE
                beta1 = beta1-step
                if beta1 < 0:
                    beta1 = max_beta
                MSE = quick_model(coords,coldata,beta1,beta2,n) 
                temp[l,0] = beta1
                temp[l,1] = beta2
                temp[l,2] = MSE

    # is it a minimum? could it be a global minimum?
        if round((beta1+step)*precise) == round(prev_beta1*precise):
            beta1=prev_beta1
            c1=1
        else:
            prev_beta1 = beta1+step
            beta1=prev_beta1
            c2=0
            c1=0

            prev_MSE = 2*imsig;
            MSE = imsig;
    
        if c2 != 1:
            beta2= prev_beta2+5*step
            if beta2 > max_beta:
                beta2 = max_beta
            while (MSE < prev_MSE):
                l+=1            
                prev_MSE = MSE
                beta2 = beta2-step 
                if beta2 < 0:
                    beta2 = max_beta
                MSE = quick_model(coords,coldata,beta1,beta2,n)
                temp[l,0] = beta1
                temp[l,1] = beta2
                temp[l,2] = MSE
    
    # is it a minimum? is it a global minimum?
        if round((beta2+step)*precise) == round(prev_beta2*precise):
            beta2=prev_beta2
            c2 = 1
        else:
            prev_beta2 = beta2+step
            beta2=prev_beta2
            c2=0
            c1=0

        done = c1+c2
        prev_MSE = 2*imsig
        MSE = imsig 
        if (second == 0 & done == 2):
            done = 0
            second = 1
            step = ang_res
            precise = fine
        if (second == 1 & done == 2):
            done == 2 
    
    temp=temp[np.argsort(temp[:,2])]
    beta1 = temp[0,0]
    beta2 = temp[0,1]

    return (beta1, beta2, l, temp)

