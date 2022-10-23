#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 10:15:37 2021

@author: dutta26
"""


from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import sys
from scipy.ndimage import rotate
import pandas as pd
import os
import helper
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import math,sys
import subprocess

def make(band):
    #filename = '/scratch/halstead/d/dutta26/m_38/coadd.fits'
    #filename = str(sys.argv[1])
    #Read the catalog 
    catalog = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/test.cat'
    #catalog = str(sys.argv[2])
    with open(catalog) as f:
        content = f.readlines()
        
    Raindex=5
    Decindex=6
    
    #Make an array containing x and y indices 
    raList =[]
    decList =[]
    for j in range(len(content)):
        #Skip the comment lines
        if((content[j].split()[0]) == '#'):
            continue
        raList.append(float(content[j].split()[Raindex])) 
        decList.append( float(content[j].split()[Decindex])) 
        
    # =============================================================================
    # tgtRa = 328.3208
    # tgtDec = 17.5833
    # for j in range(len(raList)):
    #     if(abs(tgtRa- raList[j])< 0.002 and abs(tgtDec- decList[j]) < 0.002):
    #         print (j)
    # =============================================================================
    #indices = [18430, 18131, 12933, 23442]
    indices = [18430, 23442]
    raList1  = []
    decList1 =[]
    for j in range(len(indices)):
        raList1.append(raList[indices[j]])
        decList1.append(decList[indices[j]])
    
    imgLocList = ['/scratch/halstead/d/dutta26/abell_2390/' +str(band)+'/' ]
    
    swarpLoc = '/home/dutta26/apps/bin/bin/'
    xList, yList = helper.convertToXY(raList1, decList1, '/scratch/halstead/d/dutta26/abell_2390/abell_' +str(band)+'_coadd_wted.fits')
    f=fits.open('/scratch/halstead/d/dutta26/abell_2390/abell_' +str(band)+'_coadd_wted.fits')
    data = np.array(f[0].data)
    f.close()
    sizeList =[]
    sigxxList =[]
    sigyyList =[]
    sigxyList =[]
    #Create cutouts and run the measure algorithm 
    for j in range(len(xList)):
        #print (j)
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        cut = data[y-25: y+26, x-25: x+26]
        
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measureReturnGuess(cut) 
        sizeList.append (3*psf)
        sigxxList.append(sigxx)
        sigyyList.append(sigyy)
        sigxyList.append(sigxy)
        print (flux, mux, muy, np.sqrt(e1**2 + e2**2), bkg, psf, e1, e2, sigxx, sigyy, sigxy)
    
    store = np.zeros((120, 4, 15), dtype = np.float64)
    counter = 0
    for loc in imgLocList:
        fileList = os.listdir(loc)
        for files in fileList:
            if '.weight.fits' in files:
                continue
            print (loc+files)
            
            
            
            f=fits.open(loc+files)
            back_h = float((f[0].header)['SKYBG'])
            seeing = float((f[0].header)['SEEING'])
            zp = float((f[0].header)['MAGZERO'])
            fwhm = float((f[0].header)['FWHM_FLT'])
            mjd = float((f[0].header)['MJD-MID'])
            airmass = float((f[0].header)['AIRMASS'])
            mphase = float((f[0].header)['MOONPHSE'])
            mAngle = float((f[0].header)['MOON_D'])
            
            f.close()
            wt = 100*np.power(10,(zp-25)/2.5)/((seeing)**2 * back_h)
            
    # =============================================================================
    #         if(counter != 57):
    #             counter += 1
    #             continue
    # =============================================================================
            f= open('/home/dutta26/psf_imgList.ascii', 'w+')
            f.write(loc+files)
            f.close()
            swarpCommand = './swarp @/home/dutta26/psf_imgList.ascii -c /home/dutta26/default_psfMap.swarp -IMAGEOUT_NAME /scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits'
            process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
            output, error = process.communicate()
            #First measure in the coadd 
    
            xList, yList = helper.convertToXY(raList1, decList1, '/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
            indexList= helper.detectBad('/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits', xList, yList,sizeList )
            #xList1,yList1, sizeList1= helper.createNewList(xList, yList, sizeList, indexList)
            
            f=fits.open('/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
            data = np.array(f[0].data)
            data[np.isnan(data)] = 0
            sizex, sizey = np.shape(data)
            f.close()
            
            print (indexList)
            
            for j in range(len(xList)):
    # =============================================================================
    #             if(j== 1):
    #                 break
    # =============================================================================
                #print (j)
                if(j not in indexList):
                    continue
                x = int(round( xList[j]))
                y = int(round(yList[j]))
                if(x<0 or x> sizex):
                    continue
                if(y<0 or y> sizey):
                    continue
                sigma = int(round(1.5*sizeList[j]))
                if(sigma< 20):
                    sigma = 20
                if(sigma>50):
                    sigma=50
                cut = data[y-sigma : y+sigma , x-sigma : x+sigma]    
                
                guessmux = xList[j]- x
                guessmuy = yList[j] - y
                
                hdu = fits.PrimaryHDU(cut)
                
                
                #flux, mux, muy, e1, e2, bkg, psf = measure_pythonV.measureSingleIter(cut, np.sqrt(sigxxList[j]*2), np.sqrt(sigyyList[j]*2), sigxyList[j]*2, guessmux, guessmuy) 
                flux, mux, muy, e1, e2, bkg, psf = measure_pythonV.measure(cut)
                #flux, mux, muy, e1, e2, bkg, psf = measure_pythonV.measureSingleIter(cut, np.sqrt(sigxxList[j]+sigyyList[j]), np.sqrt(sigxxList[j]+sigyyList[j]), 0, guessmux, guessmuy) 
                print ('###############################')
                print (flux, mux, muy, e1, e2, bkg, psf)
                print ('###############################')
                
                if(flux == None ):
                    continue
                store[counter, j, :] = flux, mux, muy, e1, e2, bkg, psf , mjd, back_h, seeing, zp, fwhm, airmass, mphase, mAngle
                hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/guess0/'+str(files[0:21])+'_' +str(round(flux,2))+'.fits', overwrite = True)
            counter = counter + 1    
            #if(counter >= 9):
            #    break
                
        #if(counter >= 9):   
        #    break
        #break
            #Now open the image and make the cutouts 
            
    np.save('/home/dutta26/test11_' +str(band)+'.npy', store)
    
#make('g')
make('r')
#make('i')