#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 09:23:49 2021

@author: dutta26
"""


from astropy.io import fits
import matplotlib.pyplot as plt
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
    indices = np.where((fluxList>10) & (fluxList<10000))[0]
    print (len(indices))
    raList1  = []
    decList1 =[]
    for j in range(len(indices)):
        raList1.append(raList[indices[j]])
        decList1.append(decList[indices[j]])
    
    imgLocList = ['/scratch/halstead/d/dutta26/abell_2390/' +str(band)+'/' ]
    print (len(raList1))
    swarpLoc = '/home/dutta26/apps/bin/bin/'
    xList, yList = helper.convertToXY(raList1, decList1, '/scratch/halstead/d/dutta26/abell_2390/abell_' +str(band)+'_coadd_wted.fits')
    f=fits.open('/scratch/halstead/d/dutta26/abell_2390/abell_' +str(band)+'_coadd_wted.fits')
    data = np.array(f[0].data)
    f.close()
    sizeList_star =[]
    sigxxList_star =[]
    sigyyList_star =[]
    sigxyList_star =[]
    raList_star =[]
    decList_star=[]
    fluxList_star =[]
    muxList_star=[]
    muyList_star =[]
    e1List_star =[]
    e2List_star=[]
    bkgList_star=[]
    psfList_star=[]
    
    sizeList_gal =[]
    sigxxList_gal =[]
    sigyyList_gal =[]
    sigxyList_gal =[]
    raList_gal =[]
    decList_gal=[]
    fluxList_gal =[]
    muxList_gal=[]
    muyList_gal =[]
    e1List_gal =[]
    e2List_gal=[]
    bkgList_gal=[]
    psfList_gal=[]
    
    test_flux=[]
    test_psf=[]
    #Create cutouts and run the measure algorithm 
    ellipList=[]
    e1_temp=[]
    e2_temp=[]
    counter = 0
    for j in range(len(xList)):
        #print (j)
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        cut = data[y-25: y+26, x-25: x+26]
        
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measureReturnGuess(cut) 
        
        if(flux == None):
             continue
         
# =============================================================================
#         if(np.sqrt(e1**2 + e2**2) > 0.8):
#             hdu = fits.PrimaryHDU(cut)
#             hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/guess0/coadd'+str(j) +'.fits', overwrite = True)
#         ellipList.append(np.sqrt(e1**2 + e2**2))
#         e1_temp.append(e1)
#         e2_temp.append(e2)
#         test_flux.append(flux)
#         test_psf.append(psf)
# =============================================================================
        if(flux>1500 and flux<5000 and np.sqrt(sigxx+sigyy)>3.2 and np.sqrt(sigxx+sigyy)<3.7):
            sizeList_star.append (3*psf)
            sigxxList_star.append(sigxx)
            sigyyList_star.append(sigyy)
            sigxyList_star.append(sigxy)
            raList_star.append(raList1[j])
            decList_star.append(decList1[j])
            fluxList_star.append(flux)
            muxList_star.append(mux)
            muyList_star.append(muy)
            e1List_star.append(e1)
            e2List_star.append(e2)
            bkgList_star.append(bkg)
            psfList_star.append(psf)
        else:
            sizeList_gal.append (3*psf)
            sigxxList_gal.append(sigxx)
            sigyyList_gal.append(sigyy)
            sigxyList_gal.append(sigxy)
            raList_gal.append(raList1[j])
            decList_gal.append(decList1[j])
            fluxList_gal.append(flux)
            muxList_gal.append(mux)
            muyList_gal.append(muy)
            e1List_gal.append(e1)
            e2List_gal.append(e2)
            bkgList_gal.append(bkg)
            psfList_gal.append(psf)
            
            
        #print (flux, mux, muy, np.sqrt(e1**2 + e2**2), bkg, psf, e1, e2, sigxx, sigyy, sigxy)
    print (len(raList_gal), len(raList_star))
    
    
    
    
    store = np.zeros((150, len(raList_star)+len(raList_gal), 15), dtype = np.float64)
# =============================================================================
#     corr_sigxx= sigxxList_gal - np.median(sigxxList_star)
#     corr_sigyy = sigyyList_gal - np.median(sigyyList_star)
#     corr_sigxy = sigxyList_gal - np.median(sigxyList_star)
#     e1List_gal = (corr_sigxx - corr_sigyy)/(corr_sigxx + corr_sigyy)
#     e2List_gal = 2*corr_sigxy / (corr_sigxx + corr_sigyy)
#     psfList_gal = np.sqrt(corr_sigxx + corr_sigyy)
# =============================================================================
    
    #return e1_temp, e2_temp
    #return e1List_gal, e2List_gal
    #Make the forst values the values from coadd
    for j in range(len(raList_star)):
        store[0, j, :] = fluxList_star[j], muxList_star[j], muyList_star[j], sigxxList_star[j], sigyyList_star[j], bkgList_star[j], sigxyList_star[j], 11, 11, 11, 11, 11, 11, 11, 11
    
    
    for j in range(len(raList_gal)):
        store[0, j+len(raList_star), :] = fluxList_gal[j], muxList_gal[j], muyList_gal[j], sigxxList_gal[j], sigyyList_gal[j], bkgList_gal[j], sigxyList_gal[j], 11, 11, 11, 11, 11, 11, 11, 11
    
    counter = 1
    failedArr=[]
    for loc in imgLocList:
        fileList = os.listdir(loc)
        for files in fileList:
            if '.weight.fits' in files:
                continue
            #print (loc+files)
# =============================================================================
#             if('20171109T203729.8_abell_2390_odi_r' not in files):
#                 continue
# =============================================================================
            
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
    
            xList_star, yList_star = helper.convertToXY(raList_star, decList_star, '/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
            indexList_star= helper.detectBad('/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits', xList_star, yList_star,sizeList_star )
            #xList1,yList1, sizeList1= helper.createNewList(xList, yList, sizeList, indexList)
            
            xList_gal, yList_gal = helper.convertToXY(raList_gal, decList_gal, '/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
            indexList_gal= helper.detectBad('/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits', xList_gal, yList_gal,sizeList_gal )
            
            
            
            f=fits.open('/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
            data = np.array(f[0].data)
            data[np.isnan(data)] = 0
            sizex, sizey = np.shape(data)
            f.close()
            
            current_xx =[]
            current_yy = []
            current_xy = []
            #First process the stars 
            for j in range(len(xList_star)):
                #print ('Star'+str(j))
                if((j not in indexList_star) or (sizeList_star[j] == None)):
                    continue
                x = int(round( xList_star[j]))
                y = int(round(yList_star[j]))
                if(x<0 or x> sizex):
                    continue
                if(y<0 or y> sizey):
                    continue
                sigma = int(round(1.5*sizeList_star[j]))
                if(sigma< 20):
                    sigma = 20
                if(sigma>50):
                    sigma=50
                cut = data[y-sigma : y+sigma , x-sigma : x+sigma]    
                
                guessmux = xList_star[j]- x
                guessmuy = yList_star[j] - y
                
                flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measureReturnGuess(cut)
                                
                
                if(flux == None ):
                    failedArr.append(str(j)+' _ '+str(files[0:21]))
                    continue
                
                
# =============================================================================
#                 if(e2<-1 or e2>1):
#                     hdu = fits.PrimaryHDU(cut)
#                     hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/guess0/star'+str(j) +'.fits', overwrite = True)
# 
# =============================================================================
                
                current_xx.append(sigxx)
                current_yy.append(sigyy)
                current_xy.append(sigxy)
                store[counter, j, :] = flux, mux, muy, sigxx, sigyy, bkg, sigxy , mjd, back_h, seeing, zp, fwhm, airmass, mphase, mAngle
                
            #Now do the galaxies
            avg_sigxx= np.nanmedian(current_xx)
            avg_sigyy= np.nanmedian(current_yy)
            avg_sigxy= np.nanmedian(current_xy)
            
# =============================================================================
#             print (avg_sigxx, np.median(sigxxList_star))
#             print (avg_sigyy, np.median(sigyyList_star))
#             print (avg_sigxy, np.median(sigxyList_star))
#             print (len(current_xx), len(current_yy), len(current_yy), len(xList_star), len(indexList_star))
#             
# =============================================================================
            #return
        
            for j in range(len(xList_gal)):
                #print ('Gal'+str(j))
                if((j not in indexList_gal) or (sizeList_gal[j] == None)):
                    continue
                x = int(round( xList_gal[j]))
                y = int(round(yList_gal[j]))
                if(x<0 or x> sizex):
                    continue
                if(y<0 or y> sizey):
                    continue
                sigma = int(round(1.5*sizeList_gal[j]))
                if(sigma< 20):
                    sigma = 20
                if(sigma>50):
                    sigma=50
                cut = data[y-sigma : y+sigma , x-sigma : x+sigma]    
                
                guessmux = xList_gal[j]- x
                guessmuy = yList_gal[j] - y
                
                
                
                #Fist see if normal convergence works.
                flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measureReturnGuess(cut)
                
                
                if(flux == None ):
                    if((2* ((sigxxList_gal[j]) - np.median(sigxxList_star) + avg_sigxx))< 0 ):
                       failedArr.append(str(j)+' _ '+str(files[0:21]))
                       continue
                   
                    if((2* ((sigyyList_gal[j]) - np.median(sigyyList_star) + avg_sigyy))< 0 ):
                        failedArr.append(str(j)+' _ '+str(files[0:21]))
                        continue
                    
                    alphax = np.sqrt(2* ((sigxxList_gal[j]) - np.median(sigxxList_star) + avg_sigxx))
                    alphay = np.sqrt(2* ((sigyyList_gal[j]) - np.median(sigyyList_star) + avg_sigyy))
                    alphaxy = (2* ((sigxyList_gal[j]) - np.mean(sigxyList_star) + avg_sigxy))
# =============================================================================
#                     if(j==7423):
#                         print (alphax, alphay, alphaxy, sigxxList_gal[j], np.median(sigxxList_star), avg_sigxx)
#                         hdu = fits.PrimaryHDU(cut)
#                         hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/guess0/'+str(j) + '.fits', overwrite = True)
#                         return
# =============================================================================
                    flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measureSingleIter(cut, alphax, alphay, alphaxy, guessmux, guessmuy, 0) 
                    
                    if(flux == None):
                        failedArr.append(str(j)+' _ '+str(files[0:21]))
                        continue
# =============================================================================
#                     elif(e2<-1 or e2>1):
#                         print ('*********************')
#                         print (flux, e1, e2, psf, alphax, alphay, alphaxy, guessmux, guessmuy, j)
#                         print ('*********************')
#                         hdu = fits.PrimaryHDU(cut)
#                         hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/guess0/gal_1_'+str(j) +'.fits', overwrite = True)
#                 
#                 elif(e2<-1 or e2>1):
#                     hdu = fits.PrimaryHDU(cut)
#                     hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/guess0/gal'+str(j) +'.fits', overwrite = True)
# 
# =============================================================================
                #e1 = ((sigxx - avg_sigxx) - (sigyy - avg_sigyy))/ ((sigxx - avg_sigxx) + (sigyy - avg_sigyy))
                #e2 = 2* (sigxy- avg_sigxy)/ ((sigxx - avg_sigxx) + (sigyy - avg_sigyy))
                #psf = np.sqrt((sigxx - avg_sigxx) + (sigyy - avg_sigyy))
                store[counter, j+len(raList_star), :] = flux, mux, muy, sigxx, sigyy, bkg, sigxy , mjd, back_h, seeing, zp, fwhm, airmass, mphase, mAngle
                
            
                
            counter = counter + 1    
            #if(counter >= 1):
            #    break
                
        #if(counter >= 9):   
        #    break
        #break
            #Now open the image and make the cutouts 
            
    np.save('/home/dutta26/test18_' +str(band)+'.npy', store)
    return failedArr
    
    
    
a,b=make('r')  