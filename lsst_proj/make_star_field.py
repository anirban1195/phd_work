#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 12:34:49 2022

@author: dutta26
"""


from astropy.io import fits
import numpy as np
nStar = 3000
img = np.random.normal(20, np.sqrt(20), (10000, 10000))
sigxArr=[]
for j in range(nStar):
    x_cent= np.random.randint(50, 9950)
    y_cent= np.random.randint(50, 9950)
    
    flux = np.random.randint(100,10000)
    sizex=sizey=100
    sigx = np.random.normal(6, np.sqrt(6/flux)) + (x_cent/10000) + (y_cent/10000)  + 2*np.exp( - 0.5*(( x_cent-5000)/2000)**2 - 0.5*(( y_cent-5000)/2000)**2)
   
    
    sigy = np.random.normal(6, np.sqrt(6/flux)) + (x_cent/10000) + (y_cent/10000)  + 2*np.exp( - 0.5*(( x_cent-5000)/2000)**2 - 0.5*(( y_cent-5000)/2000)**2)
    
    
    sigxy = np.random.normal(0.5, np.sqrt(0.5/flux))
    sigxArr.append(sigx)
    
    muArr= [sizex/2.0-0.5, sizey/2.0]
    cov = [[sigx, sigxy], [sigxy, sigy]]
    const = flux
    x, y = np.random.multivariate_normal(muArr, cov, const).T
    x = np.int32(np.round(x))
    y = np.int32(np.round(y))
    obj = np.zeros((sizex,sizey))
    np.add.at(obj, (y,x), 1)
    img[int(y_cent-sizey/2): int(y_cent+sizey/2), int(x_cent-sizex/2):int( x_cent+sizex/2)] += obj
    
hdu = fits.PrimaryHDU(img)  
hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/star_field.fits', overwrite=True)
                
print (max(sigxArr), min(sigxArr), np.std(sigxArr))
    