#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 09:03:56 2021

@author: anirban
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

def makeSizeMap(filename)
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
        
    
        
    #Convert star ra dec to x and y 
    #filename = '/scratch/halstead/d/dutta26/abell_2390/r/20191005T225729.3_Abell2390_odi_r.8487.fits'
    f= open('/home/dutta26/psf_imgList.ascii', 'w+')
    f.write(filename)
    f.close()
    swarpCommand = './swarp @/home/dutta26/psf_imgList.ascii -c /home/dutta26/default_psfMap.swarp -IMAGEOUT_NAME /scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits'
    process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
    output, error = process.communicate()
    
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
    indexList = helper.detectBad1(filename, star_xList, star_yList)
    
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
        values.append( star_sizeList[j])
    
    
    grid_z1 = griddata(points, values, (grid_y, grid_x), method='nearest')    
    grid_z1 = np.array(grid_z1, dtype = np.float32).T
    hdu = fits.PrimaryHDU(grid_z1) 
    hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/residue.fits', clobber=True)     
    
    










