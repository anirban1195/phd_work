#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 11:38:21 2023

@author: dutta26
"""

import numpy as np
import sys
import getgamma 
from astropy.io import fits
f=open('/scratch/bell/dutta26/wiyn_sim/master.cat')
content = f.readlines()
f.close()


z_mass = 0.3
cent_ra =  328.3941
cent_dec = 17.6697

nSources = 16000
bkg = 100
img = np.random.normal(bkg, np.sqrt(bkg), (10000, 10000))
psf_xx = 2
psf_yy = 2
psf_xy = 0.0
sizex = sizey = 100
for j in range(nSources):
    print (j)
    #Simulate star
    if(np.random.random() <= 0.3):
        flux = int(np.random.normal(50000, 5000))
        rand1 = np.random.normal(0,0.1)
        rand2 = np.random.normal(0,0.1)
        rand3 = np.random.normal(0,0.01)
        cov = [[2*psf_xx, 2*psf_xy+rand3],[2*psf_xy+rand3 , 2*psf_yy+rand2]]
        
        xLoc = int(np.random.uniform(100, 9950))
        yLoc = int(np.random.uniform(100, 9950))
        muArr= [sizex/2.0-0.5, sizey/2.0-0.5]
        
        x, y = np.random.multivariate_normal(muArr, cov, flux).T
        x = np.int32(np.round(x))
        y = np.int32(np.round(y))
        obj = np.zeros((sizex,sizey))
        np.add.at(obj, (y,x), 1)
        img[yLoc-50:yLoc+50, xLoc-50:xLoc+50] += obj
    else: #If galaxy
        flux = int(np.random.normal(25000, 500))
        rand1 = np.random.normal(0,0.1)
        rand2 = np.random.normal(0,0.1)
        rand3 = np.random.normal(0,0.01)
        xLoc = int(np.random.uniform(100, 9950))
        yLoc = int(np.random.uniform(100, 9950))
        
        ra = cent_ra - (yLoc - 5000)* (0.11/3600)
        dec = cent_dec + (xLoc - 5000)* (0.11/3600)
        g1, g2, kappa = getgamma.getGamma(cent_ra, cent_dec, 0.3, ra, dec, 2, 0.25)
        #g1=g1*2
        #g2=g2*2
        
        galxx = np.random.uniform(10,20)
        galyy = np.random.uniform(10,20)
        avg = abs(galxx+galyy)/2
        galxy = np.random.uniform(-0.5*avg, +0.5*avg)
        cov1 = [[2*galxx, 2*galxy],[2*galxy , 2*galyy]]
        shear_mat = [[1-g1, -g2],[-g2, 1+g1]]
        psf = [[2*psf_xx, 2*psf_xy+rand3],[2*psf_xy+rand3 , 2*psf_yy+rand2]]
        final_cov = np.add(cov1,psf)
        muArr= [0, 0]
        x, y = np.matmul(np.random.multivariate_normal(muArr, final_cov, flux), shear_mat).T
        x = np.int32(np.round(x)+sizex/2.0-0.5)
        y = np.int32(np.round(y)+sizey/2.0-0.5)
        obj = np.zeros((sizex,sizey))
        np.add.at(obj, (y,x), 1)
        img[yLoc-50:yLoc+50, xLoc-50:xLoc+50] += obj
        
# =============================================================================
# hdu = fits.PrimaryHDU(img)
# hdu.writeto('/scratch/bell/dutta26/wiyn_sim/temp.fits', overwrite=True)
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'RADESYS', value='ICRS')
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CUNIT1', value='deg')
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CUNIT2', value='deg')  
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CTYPE1', value='DEC--TAN')
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CTYPE2', value='RA---TAN')   
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CRVAL2', value=328.3941)
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CRVAL1', value=17.6697)   
#         
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CRPIX1', value=5000)
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CRPIX2', value=5000)   
# 
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CD2_2', value=-3.055543857045E-05)
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CD1_2', value=0)   
# 
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CD2_2', value=0)
# fits.setval('/scratch/bell/dutta26/wiyn_sim/temp.fits' , 'CD1_1', value=3.055543857045E-05)   
# =============================================================================
f=fits.open('/scratch/bell/dutta26/wiyn_sim/temp.fits', mode='update')
f[0].data= img
f.flush()