#! /usr/bin/env python
# Shape_functions.py
# Contains a number of functions used in the shapelet main code

import io
import math as m
import numpy as np
import astropy
import pyfits

#######################################################################
## Attempts to optimise beta1, beta2: non-optimal (and dumb - see temp.py for a smarter version) routinue for 2 dependent varibles

def beta1_beta2(coords,coldata,n,ext_source,ang_res):
    
    data_size = coldata.shape
    nside = m.sqrt(data_size[0])
    fine = 10e6        # precision of the two stages
    coarse = 10e4
    iters = 200
    temp = np.ones((iters,3))
    done = 0 
    first = 1
    second = 0

    # Coarse loop to isolate beta 1 and 2: Minimise beta1, then beta2. If beta1 and beta2 are the same as last minimisation loop, convergence is assumed. Fine loop just adds extra decimal places.
    max_beta = 0.9*pow(n,-0.52)*ext_source*ang_res*4                          
    prev_beta1 = max_beta
    prev_beta2 = max_beta
    coarse_step = 4
    step_size = (max_beta-ang_res)/coarse_step

    l=-1
    while done == 0 & l<iters-1:
        if first == 1:
            precise = coarse
            step_size = (max_beta-ang_res)/coarse_step
            beta1_start = ang_res
            beta2_start = ang_res
            beta1 = beta1_start
            beta2 = beta2_start           
        for step_number in range(0,9):
            l+=1
            if l>iters-1:
                break
            beta1 += step_size
            MSE = quick_model(coords,coldata,beta1,beta2,n) 
            temp[l,0] = beta1
            temp[l,1] = beta2
            temp[l,2] = MSE
            
        temp = temp[np.argsort(temp[:,2])]
        beta1 = temp[0,0]
        beta2 = temp[0,1]
        
        if round(prev_beta1*precise) == round(beta1*precise):
            c1 = 1
        else:
            prev_beta1 = beta1
            beta2 = beta2_start
            c1 = 0
        
        for step_number in range(0,9):
            l+=1
            if l>iters-1:
                break
            beta2 += step_size
            MSE = quick_model(coords,coldata,beta1,beta2,n)
            temp[l,0] = beta1
            temp[l,1] = beta2
            temp[l,2] = MSE
            
        temp = temp[np.argsort(temp[:,2])]
        beta1 = temp[0,0]
        beta2 = temp[0,1]
        
        if round(prev_beta2*precise) == round(beta2*precise):
            c2 = 1
        else:
            prev_beta2 = beta2
            beta1 = beta1_start
            c2 = 0
        if first == 1:
            finish = c1+c2
            if finish == 2:
                first = 0
                precise = fine
                step_size = ang_res
                if beta1-5*ang_res < 0:
                    beta1 = 6*ang_res
                if beta2-5*ang_res < 0:
                    beta2 = 6*ang_res
                beta1_start = beta1-5*ang_res
                beta2_start = beta2-5*ang_res
                beta1 = beta1_start
                beta2 = beta2_start
                second = 1
                c1 = 0
                c2 = 0
        if second == 1:   
            done = c1+c2
   
    temp=temp[np.argsort(temp[:,2])]
    beta1 = temp[0,0]
    beta2 = temp[0,1]
   
    return (beta1, beta2, l, temp)

#######################################################################
## Determines the shapelet moments

def deconstruct(coords, coldata, beta1, beta2, nmax):

    min_percent = 0.0    # zero small coefficients
    tot_moms = 0.5*(nmax+1)*(nmax+2)
    f_coeffs = np.zeros((tot_moms,3))
    data_size = coldata.shape
    npix = data_size[0]
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
                    i=i+1
                    f_coeffs[i,0]=n1
                    f_coeffs[i,1]=n2
                    f_coeffs[i,2]=fhat
                else:
                    i=i+1
                    f_coeffs[i,0]=n1
                    f_coeffs[i,1]=n2
                    f_coeffs[i,2]=0.0    
    
    return f_coeffs

#######################################################################
## Creates the orthogonal basis functions

def basis(coords,b1,b2,nmax):

    hermites = np.loadtxt('hermite_coeffs.txt')
    xrot = coords[:,0]/b1
    yrot = coords[:,1]/b2
    data_size = coords.shape
    npix = data_size[0]
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

    PSNR = 20*np.log10(G)-10*np.log10(NMSE)
    performance[0] = NMSE
    performance[1] = PSNR

    return performance

#######################################################################
## Handles the function calls of beta1_beta2

def quick_model(coords,coldata,beta1,beta2,n):

    data_size = coords.shape
    npix = data_size[0]
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

    data_size = coords.shape
    npix = data_size[0]
    mom_size = moments.shape
    tot_coeffs = mom_size[0]
    n = int(max(moments[:,0]))
    Mpiece = np.zeros((npix,1))
    colmod = np.zeros((npix,1))

    M = basis(coords,beta1,beta2,n)
    
    for i in range(0,tot_coeffs):
        n1 = int(moments[i,0])
        n2 = int(moments[i,1])
        f_hat = moments[i,2]
        Mpiece[:,0] = M[:,n1,n2]
        colmod += f_hat*(Mpiece)
    
    return colmod

#######################################################################

def major_minor(coords,coldata,n,ext_source,ang_res):
    
    data_size = coldata.shape
    nside = m.sqrt(data_size[0])
    max_iters = 100
    fine = 10e6        # precision of the two stages
    coarse = 10e4
    temp = np.ones((10,3))
    done = 0 
    first = 1
    second = 0

    # Coarse loop to isolate beta 1 and 2: Minimise beta1, then beta2. If beta1 and beta2 are the same as last minimisation loop, convergence is assumed.
    max_beta = 0.9*pow(n,-0.52)*nside/2.*ang_res
    min_beta = ang_res                          
    prev_beta1 = min_beta  
    prev_beta2 = min_beta  
    l=-1

    while done == 0:
        imsig = 100                                
        MSE = imsig                                
        prev_MSE = 2*imsig  
        beta2 = prev_beta2 
        change = 0                                 
        c1=0                                       
        c2=0                         
        if first == 1:
            step = 10*ang_res 
            precise = coarse
        if second == 1:
            step = ang_res  
            precise = fine
    
        while (change < 2 & l<max_iters):
            if c1 != 1:
                beta1 = prev_beta1-5*step
                if beta1 < 0:
                    beta1 = min_beta
            while (MSE < prev_MSE):
                l+=1 
                prev_MSE = MSE
                beta1 = beta1+step
            if beta1 > max_beta:
                beta1 = min_beta
            MSE = quick_model(coords,coldata,beta1,beta2,n) 
            temp[l,0] = beta1
            temp[l,1] = beta2
            temp[l,2] = MSE

    # is it a minimum? could it be a global minimum?
            if round((beta1-step)*precise) == round(prev_beta1*precise):
                beta1=prev_beta1
                c1=1
            else:
                prev_beta1 = beta1-step
                beta1=prev_beta1
                c2=0
                c1=0

            prev_MSE = 2*imsig;
            MSE = imsig;
    
            if c2 != 1:
                beta2= prev_beta2-5*step
                if beta2 < 0:
                    beta2 = min_beta
            while (MSE < prev_MSE):
                l+=1            
                prev_MSE = MSE
                beta2 = beta2+step 
                if beta2 > max_beta:
                    beta2 = min_beta
            MSE = quick_model(coords,coldata,beta1,beta2,n)
            temp[l,0] = beta1
            temp[l,1] = beta2
            temp[l,2] = MSE
    
    # is it a minimum? is it a global minimum?
        if round((beta2-step)*precise) == round(prev_beta2*precise):
            beta2=prev_beta2
            c2 = 1
        else:
            prev_beta2 = beta2-step
            beta2=prev_beta2
            c2=0
            c1=0

        change = c1+c2
        prev_MSE = 2*imsig
        MSE = imsig 

        # because I don't entirely trust my algorythm to find the optimal..    
        temp=temp[np.argsort(temp[:,2])]
        prev_beta1 = temp[0,0]
        prev_beta2 = temp[0,1]
        beta2 = prev_beta2
        if second == 1:
            done = 1
        else:
            second = 1
            first = 0
   
    l+=1
    temp[l,0]=beta1
    temp[l,1]=beta2
    temp[l,2]=MSE   
    
    temp=temp[np.argsort(temp[:,2])]
    beta1 = temp[0,0]
    beta2 = temp[0,1]

    return (beta1, beta2, l, temp)