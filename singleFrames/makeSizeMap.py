#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 12:28:03 2021

@author: dutta26
"""

import subprocess
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord   
from astropy.coordinates import FK5 
import measure_pythonV
from astropy.stats import sigma_clipped_stats
from scipy.interpolate import griddata
import sys,helper,os

def makeSizeMap(size):
    #Read the starList

    star_catalog = '/home/dutta26/codes/makeCuts/starList.txt'
    swarpLoc = '/home/dutta26/apps/bin/bin/'
    star_raArr=[]
    star_decArr =[]
    with open(star_catalog) as f:
        content = f.readlines()
    
    for j in range(len(content)):
        star_raArr.append(float(  (content[j].split(','))[0]) )
        star_decArr.append(float(  (content[j].split(','))[1]) )
        
    
        
    
    filename = '/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits'
    f=fits.open(filename)
    img = np.array(f[0].data)
    w=wcs.WCS(f[0].header)
    sizey, sizex = np.shape(img)
    f.close()
    
    temp = np.zeros((len(star_raArr), 2), dtype = np.float32)
    temp[:,0] = star_raArr
    temp[:,1] = star_decArr
    
    
    temp1 = w.wcs_world2pix(temp, 0)
    star_xList = temp1[:,0]
    star_yList = temp1[:,1]
    if(size < 6):
        size = 6
    sizeList = np.ones((len(star_yList)) , dtype =np.float32)*size*1.5
    indexList = helper.detectBadListV(filename, star_xList, star_yList, sizeList)
    
    star_xx=[]
    star_yy=[]
    star_xy=[]
    star_flux=[]
    star_xList1=[]
    star_yList1 =[]
    star_err =[]
    for j in range(len(star_xList)):
        
        if(j not in indexList):
            continue
        y=int(star_yList[j])
        x=int(star_xList[j])
        cut = img[y - 20: y+20, x-20:x+20]
        #Check if the cuts are good ti be used
        flag = helper.detectBad(cut)
        if(flag == 1):
            continue 
            
        flux, mux, muy, e1, e2, back, size, sigxx, sigyy, sigxy = measure_pythonV.measureReturnGuess(cut)
        
        if(sigxx == None or flux == None or sigxx<0 or sigyy<0):
            continue
        ellip = np.sqrt(e1**2 + e2**2)
        star_xx.append(sigxx)
        star_yy.append(sigyy)
        star_xy.append(sigxy)
        star_flux.append(flux)
        star_xList1.append(x)
        star_yList1.append(y)
        star_err.append(0)
        del cut
    
    
    
    fact = 5
    star_xx = np.array(star_xx)
    star_yy = np.array(star_yy)  
    star_xy = np.array(star_xy)  
    star_e1 = (star_xx - star_yy)/(star_xx+star_yy)
    star_e2 = (2*star_xy)/(star_xx+star_yy)
    #sys.exit()     
    star_xList1 = np.array(star_xList1) /fact
    star_yList1 = np.array(star_yList1) /fact 
    
    star_sizeList = np.sqrt(star_xx  + star_yy)
    
    #Get the rage of values of x and y
    xShape = int(sizex/fact) + 10
    yShape = int(sizey/fact) + 10
    
    points=[]
    values =[] 
    grid_y, grid_x = np.meshgrid(np.linspace(0, yShape-1, yShape), np.linspace(0, xShape-1, xShape))
    
    for j in range(len(star_sizeList)- 1):
        x = int(round(star_xList1[j]))
        y = int(round(star_yList1[j]))
        
        points.append([y,x])
        values.append( star_xx[j])
    
    
    grid_z1 = griddata(points, values, (grid_y, grid_x), method='nearest')    
    grid_z1 = np.array(grid_z1, dtype = np.float32).T
    hdu = fits.PrimaryHDU(grid_z1) 
    hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/sigxx.fits', clobber=True)  
    
    
    #Get the rage of values of x and y
    xShape = int(sizex/fact) + 10
    yShape = int(sizey/fact) + 10
    
    points=[]
    values =[] 
    grid_y, grid_x = np.meshgrid(np.linspace(0, yShape-1, yShape), np.linspace(0, xShape-1, xShape))
    
    for j in range(len(star_sizeList)- 1):
        x = int(round(star_xList1[j]))
        y = int(round(star_yList1[j]))
        
        points.append([y,x])
        values.append( star_yy[j])
    
    
    grid_z1 = griddata(points, values, (grid_y, grid_x), method='nearest')    
    grid_z1 = np.array(grid_z1, dtype = np.float32).T
    hdu = fits.PrimaryHDU(grid_z1) 
    hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/sigyy.fits', clobber=True)  
    
    
    #Get the rage of values of x and y
    xShape = int(sizex/fact) + 10
    yShape = int(sizey/fact) + 10
    
    points=[]
    values =[] 
    grid_y, grid_x = np.meshgrid(np.linspace(0, yShape-1, yShape), np.linspace(0, xShape-1, xShape))
    
    for j in range(len(star_sizeList)- 1):
        x = int(round(star_xList1[j]))
        y = int(round(star_yList1[j]))
        
        points.append([y,x])
        values.append( star_xy[j])
    
    
    grid_z1 = griddata(points, values, (grid_y, grid_x), method='nearest')    
    grid_z1 = np.array(grid_z1, dtype = np.float32).T
    hdu = fits.PrimaryHDU(grid_z1) 
    hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/sigxy.fits', clobber=True)  
    










