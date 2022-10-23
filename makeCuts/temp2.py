#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 22:19:31 2021

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
from scipy.signal import convolve2d as conv2
from skimage import color, data, restoration

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
    fluxIndex = 1
    
    #Make an array containing x and y indices 
    raList =[]
    decList =[]
    fluxList =[]
    for j in range(len(content)):
        #Skip the comment lines
        if((content[j].split()[0]) == '#'):
            continue
        raList.append(float(content[j].split()[Raindex])) 
        decList.append(float(content[j].split()[Decindex])) 
        fluxList.append(float(content[j].split()[fluxIndex])) 
        
    fluxList = np.array(fluxList)
    #indices = [18430, 18131, 12933, 23442]
    #indices = [18430, 23442, 18131, 12933]
    indices = np.where(fluxList>100)[0]
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
    failedArr=[]
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
#             if(counter != 12):
#                 counter += 1
#                 continue
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
            success = 10 
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
                
                
                if(j== 0 or j == 1):
                    flux, mux, muy, e1, e2, bkg, psf = measure_pythonV.measure(cut)
                    if(flux == None ):
                        failedArr.append(str(j)+' _ '+str(files[0:21]))
                        hdu = fits.PrimaryHDU(cut)
                        hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/guess0/'+str(j) + '_' + str(files[0:21])+'_' +str(mjd)+'.fits', overwrite = True)
                        continue
                    flux_current, mux_current, muy_current, e1_current, e2_current, bkg_current, psf_current, sigxx_current, sigyy_current, sigxy_current = measure_pythonV.measureReturnGuess(cut)
                    success = j
                    print ('###############################')
                    print (flux, mux, muy, e1, e2, bkg, psf)
                    print ('###############################')
                if( (success == 1 or success == 0) and (j== 2 or j == 3)):  
                    

                    
                    #First deconvolve original gal psf from star psf 
                    x = np.linspace(0, 100-1, 100)
                    y = np.linspace(0, 100-1, 100)
                    x= x -100/2.0 + 0.5 
                    y= y -100/2.0 + 0.5 
                    #mux=muy= np.random.normal(29.5, 0.15)
                    x, y = np.meshgrid(x, y)
                    
                    
                    sigmaxx_galDecon = (sigxxList[j] - sigxxList[success])
                    sigmayy_galDecon = (sigyyList[j] - sigyyList[success])
                    sigmaxy_galDecon = (sigxyList[j] - sigxyList[success])
                    alphax = np.sqrt(2*sigmaxx_galDecon)
                    alphay = np.sqrt(2*sigmayy_galDecon)
                    alphaxy = 2*sigmaxy_galDecon
                    A =10e4/(2*np.pi*alphax*alphay*np.sqrt(1- (alphaxy/(alphax*alphay))**2 ))
                    arb_const = 2*(1 - (alphaxy/ (alphax*alphay))**2 )
                    
                    deconvolved_galPsf = (A * np.exp(-((x)**2/(arb_const*alphax**2)+ (y)**2/(arb_const*alphay**2)  - 2*alphaxy*(y)*(x)/(arb_const*alphax**2 * alphay**2 ) )))
                    hdu = fits.PrimaryHDU(deconvolved_galPsf)
                    hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/deconGal.fits', overwrite = True)
                    
                    #Now convolve the deconvolved gal psf with this frames psf
                    sigmaxx_galCurrent = (sigxxList[j] - sigxxList[success]+sigxx_current)
                    sigmayy_galCurrent = (sigyyList[j] - sigyyList[success]+ sigyy_current)
                    sigmaxy_galCurrent = (sigxyList[j] - sigxyList[success] + sigxy_current)
                    
                    
                    alphax = np.sqrt(2*sigmaxx_galCurrent)
                    alphay = np.sqrt(2*sigmayy_galCurrent)
                    alphaxy = sigmaxy_galCurrent*2
                    A =10e4/(2*np.pi*alphax*alphay*np.sqrt(1- (alphaxy/(alphax*alphay))**2 ))
                    arb_const = 2*(1 - (alphaxy/ (alphax*alphay))**2 )
                    
                    galShape = (A * np.exp(-((x)**2/(arb_const*alphax**2)+ (y)**2/(arb_const*alphay**2)  - 2*alphaxy*(y)*(x)/(arb_const*alphax**2 * alphay**2 ) )))
                    
                    hdu = fits.PrimaryHDU(galShape)
                    hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/galConPsf.fits', overwrite = True)
                    
                    #flux_gal, mux_gal, muy_gal, e1_gal, e2_gal, bkg_gal, psf_gal, sigxx_gal, sigyy_gal, sigxy_gal = measure_pythonV.measureReturnGuess(galShape*1000)
                    
                    attemptNo = 0
                    while(attemptNo <= 0):
                        #Now measure the galaxy shape with the new params
                        flux, mux, muy, e1, e2, bkg, psf = measure_pythonV.measureSingleIter(cut, alphax, alphay, alphaxy, guessmux, guessmuy, attemptNo) 
                        attemptNo += 1
                        if (flux != None):
                            break
                    
                    if(flux == None ):
                        failedArr.append(str(j)+' _ '+str(files[0:21]))
                        hdu = fits.PrimaryHDU(cut)
                        hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/guess0/'+str(j) + '_' + str(files[0:21])+'_' +str(mjd)+'.fits', overwrite = True)
                        continue
                    print ('###############################')
                    print (flux, mux, muy, e1, e2, bkg, psf)
                    print ('###############################')
                
                
                
                
                store[counter, j, :] = flux, mux, muy, e1, e2, bkg, psf , mjd, back_h, seeing, zp, fwhm, airmass, mphase, mAngle
                
                
            counter = counter + 1    
            #if(counter >= 1):
            #    break
                
        #if(counter >= 9):   
        #    break
        #break
            #Now open the image and make the cutouts 
            
    np.save('/home/dutta26/test14_' +str(band)+'.npy', store)
    return failedArr
    
#print(make('g'))
#print (make('r'))
print(make('i'))