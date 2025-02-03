#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 09:25:58 2023

@author: dutta26
"""

import numpy as np
from astropy.io import fits
import helper_phosim,correct_wcs, helper
import subprocess,sys, helper
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt

lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
# =============================================================================
# catalog = "/scratch/bell/dutta26/wiyn_sim/coadd_sim.cat"
# coadd_file = '/scratch/bell/dutta26/wiyn_sim/wted_coadds/ir_coadd_wt.fits'
# outFile = '/home/dutta26/codes/wiyn_wl_sim/source_list.pk1'
# =============================================================================

# =============================================================================
# catalog = "/scratch/bell/dutta26/wiyn_sim/coadd_sim_perfectGauss.cat"
# coadd_file = '/scratch/bell/dutta26/wiyn_sim/temp.fits'
# outFile = '/home/dutta26/codes/wiyn_wl_sim/source_list_perfectGauss.pk1'
# =============================================================================
def seperate(catalog, coadd_file, outFile, plotLoc, psfMin, psfMax, fluxMin, fluxMax, psfMin_2, psfMax_2, fluxMin_2, fluxMax_2):
    f=fits.open(coadd_file)
    data = f[0].data
    f.close()
    sizex, sizey= np.shape(data)
    #Read the catalog 
    with open(catalog) as f:
        content = f.readlines()
    
    Raindex=5
    Decindex=6
    xxIndex = 7
    yyIndex = 8
    xyIndex = 9
    
    #Make an array containing x and y indices 
    raList =[]
    decList =[]
    sex_xx_list = []
    sex_yy_list = []
    sex_xy_list = []
    for j in range(len(content)):
        #Skip the comment lines
        if((content[j].split()[0]) == '#'):
            continue
        
        raList.append(float(content[j].split()[Raindex])) 
        decList.append(float(content[j].split()[Decindex])) 
        sex_xx_list.append(float(content[j].split()[xxIndex]))
        sex_yy_list.append(float(content[j].split()[yyIndex]))
        sex_xy_list.append(float(content[j].split()[xyIndex]))
    
    xList, yList = helper.convertToXY(decList, raList, coadd_file)
    fluxArr=[]
    sizeArr=[]
    starArr=[]
    
    #Now run all through measure 
    for j in range(len(raList)):
        print (j)
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        cut = data[y-25: y+25, x-25: x+25]
        if(x<30 or y<30 or x>(sizex-30) or y>(sizey-30)):
            fluxArr.append(0)
            sizeArr.append(0)
            starArr.append(0)
            continue
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2)
        if(fluxArr == None or sigxx==None or sigyy==None or np.isnan(psf)):
            fluxArr.append(0)
            sizeArr.append(0)
            starArr.append(0)
            continue
        fluxArr.append(flux)
        sizeArr.append(np.sqrt(sigxx+sigyy))
        if(np.log10(flux)>np.log10(fluxMin) and np.log10(flux)<np.log10(fluxMax) and psf>psfMin and psf<psfMax):
            starArr.append(1)
        elif(np.log10(flux)>np.log10(fluxMin_2) and np.log10(flux)<np.log10(fluxMax_2) and psf>psfMin_2 and psf<psfMax_2):
            starArr.append(2)
        else:
            starArr.append(0)
            
    #Make a pandas dataframe for the sources . 1 for star and 0 for galaxy
    a=np.array([raList, decList, sex_xx_list, sex_yy_list, sex_xy_list, starArr])
    b=np.swapaxes(a, 0,1)
    df_source = pd.DataFrame(b,  columns = ['ra', 'dec', 'sex_xx', 'sex_yy', 'sex_xy', 'star_bool'])
    
    #df_source.to_pickle(outFile)     
    plt.subplot()
    plt.yscale('log')
    plt.plot(sizeArr, fluxArr,'b.', markersize= 2)  
    plt.xlabel('Size (in pixels)')
    plt.ylabel('Counts')
    plt.plot([psfMax, psfMin], [fluxMax, fluxMax], 'r-',markersize= 2)
    plt.plot([psfMax, psfMin], [fluxMin, fluxMin], 'r-',markersize= 2)
    plt.plot([psfMax, psfMax], [fluxMax, fluxMin], 'r-',markersize= 2)
    plt.plot([psfMin, psfMin], [fluxMax, fluxMin], 'r-',markersize= 2)
    
    plt.plot([psfMax_2, psfMin_2], [fluxMax_2, fluxMax_2], 'k-',markersize= 2)
    plt.plot([psfMax_2, psfMin_2], [fluxMin_2, fluxMin_2], 'k-',markersize= 2)
    plt.plot([psfMax_2, psfMax_2], [fluxMax_2, fluxMin_2], 'k-',markersize= 2)
    plt.plot([psfMin_2, psfMin_2], [fluxMax_2, fluxMin_2], 'k-',markersize= 2)
    
    plt.savefig(plotLoc+'flux_vs_size.png')
    plt.close()
