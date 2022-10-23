#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 10:50:10 2021

@author: anirban
"""

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
#from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord   
from astropy.coordinates import FK5 
import measure_pythonV
import math,sys
import subprocess
import helper1,makeSizeMap, helper, os
from astropy.stats import sigma_clipped_stats

#Read the catalog 
catalog = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/test.cat'
star_catalog = '/home/dutta26/codes/makeCuts/starList.txt'




raArr =[] 
decArr = []
fluxArr = []
xPosArr =[]
yPosArr =[]
sizeArr = []

#Read catalog
with open(catalog) as f:
    content = f.readlines()
    
#Remove stars and large objects from catalog    
indexList = helper1.removeLarge(np.arange(len(content)), catalog)
indexList = helper1.removeStars(indexList, catalog)

for j in range(len(content)):
    if((content[j].split())[0] == '#' or (j not in indexList)):
       continue
    raArr.append(float(  (content[j].split())[5]) )
    decArr.append(float(  (content[j].split())[6]) )
    fluxArr.append(float(  (content[j].split())[1]) )
    xPosArr.append(float(  (content[j].split())[3]) )
    yPosArr.append(float(  (content[j].split())[4]) )
    size1 = float((content[j].split())[7])
    size2 = float((content[j].split())[8])
    sizeArr.append(max(size1, size2))


star_raArr=[]
star_decArr =[]
with open(star_catalog) as f:
    content = f.readlines()

for j in range(len(content)):
    star_raArr.append(float(  (content[j].split(','))[0]) )
    star_decArr.append(float(  (content[j].split(','))[1]) )


#Convert star ra dec to x and y 
f=fits.open('/scratch/halstead/d/dutta26/abell_2390/abell_ir_coadd_wted1.fits')
img = np.array(f[0].data)
w=wcs.WCS(f[0].header)
sizex, sizey = np.shape(img)
f.close()

temp = np.zeros((len(star_raArr), 2), dtype = np.float32)
temp[:,0] = star_raArr
temp[:,1] = star_decArr


temp1 = w.wcs_world2pix(temp, 0)
star_xList = temp1[:,0]
star_yList = temp1[:,1]

#Now measure star properties and store them 
star_xx=[]
star_yy=[]
star_xy=[]
star_flux=[]
star_xList1=[]
star_yList1 =[]
for j in range(len(star_xList)):
    y=int(star_yList[j])
    x=int(star_xList[j])
    cut = img[y - 15: y+15, x-15:x+15]
    if(helper.vert_stripe(cut) == 1):
        continue
    flux, mux, muy, e1, e2, back, size, sigxx, sigyy, sigxy = measure_pythonV.measure(cut)
    if(sigxx == None):
        continue
    star_xx.append(sigxx)
    star_yy.append(sigyy)
    star_xy.append(sigxy)
    star_flux.append(flux)
    star_xList1.append(x)
    star_yList1.append(y)
    del cut

star_xx = np.array(star_xx)
star_yy = np.array(star_yy)  
star_xy = np.array(star_xy)   
star_xList1 = np.array(star_xList1) 
star_yList1 = np.array(star_yList1) 
star_store = np.zeros((len(star_xx), 6))
star_store[:,0] = star_xList1
star_store[:,1] = star_yList1
star_store[:,2] = star_xx
star_store[:,3] = star_yy
star_store[:,4] = star_xy

#Now measure objects and correct for PSF
store = np.zeros((len(xPosArr), 12), dtype =np.float32)    
cnt=cnt1 = 0
for j in range(len(xPosArr)):
    x = int(xPosArr[j])
    y = int(yPosArr[j])
    size = int(sizeArr[j])
    if(size < 6.5):
        size = 6.5
    cut = img[y - int(2.5*size): y+int(2.5*size), x-int(2.5*size):x+int(2.5*size)]
    flux, mux, muy, e1, e2, back, size, sigxx, sigyy, sigxy = measure_pythonV.measure(cut)
    if(flux == None or sigxx == None or sigxy == None or mux == None):
        continue
    temp = np.copy(star_store)
    temp[:,0] = temp[:, 0]-x
    temp[:,1] = temp[:, 1]-y
    temp[:,5 ]= np.sqrt( temp[:, 0] **2 +  temp[:, 1] **2 )
    temp = temp[temp[:,5].argsort()]
    #If same star then continue 
    if(temp[0, 5]< 5):
        cnt += 1
        continue
    avgSigxx = np.sum( temp[0:10, 2] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
    avgSigyy = np.sum( temp[0:10, 3] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
    avgSigxy = np.sum( temp[0:10, 4] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
    
    store[j,:] = flux, mux, muy, e1, e2, back, size, sigxx-avgSigxx, sigyy - avgSigyy, sigxy-avgSigxy, raArr[j], decArr[j] 

    del cut
    
np.save( '/home/dutta26/allObj_params.npy', store)
del img    

store =np.load('/home/dutta26/allObj_params.npy')
swarpLoc = '/home/dutta26/apps/bin/bin/'
bandList = ['z']
for band in bandList:
    storei = np.zeros((len(xPosArr), 200, 25) , dtype = np.float32)
    bandLoc = '/scratch/halstead/d/dutta26/abell_2390/'+band
    frameCount = 0
    for file in os.listdir(bandLoc):
# =============================================================================
#         if('00832.7' not in file):
#             continue
# =============================================================================
        #Skip weight files
        if('.weight' in file):
            continue
        
        #Swarp the file 
        filename = bandLoc +'/'+ file
        f=fits.open(filename)
        back_h = float((f[0].header)['SKYBG'])
        seeing = float((f[0].header)['SEEING'])
        zp = float((f[0].header)['MAGZERO'])
        fwhm = float((f[0].header)['FWHM_FLT'])
        mjd = float((f[0].header)['MJD-MID'])
        airmass = float((f[0].header)['AIRMASS'])
        mphase = float((f[0].header)['MOONPHSE'])
        mAngle = float((f[0].header)['MOON_D'])
        f.close()
        
        f= open('/home/dutta26/psf_imgList.ascii', 'w+')
        f.write(filename)
        f.close()
        swarpCommand = './swarp @/home/dutta26/psf_imgList.ascii -c /home/dutta26/default_psfMap.swarp -IMAGEOUT_NAME /scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits'
        process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
        output, error = process.communicate()
        
        makeSizeMap.makeSizeMap(seeing/0.11)
        
        f=fits.open('/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
        img = np.array(f[0].data)
        w=wcs.WCS(f[0].header)
        imgsizey, imgsizex = np.shape(img)
        f.close()
        
        f=fits.open('/scratch/halstead/d/dutta26/abell_2390/sigxx.fits')
        sigxx_Map = np.array(f[0].data)
        f.close()
        f=fits.open('/scratch/halstead/d/dutta26/abell_2390/sigyy.fits')
        sigyy_Map = np.array(f[0].data)
        f.close()
        f=fits.open('/scratch/halstead/d/dutta26/abell_2390/sigxy.fits')
        sigxy_Map = np.array(f[0].data)
        f.close()
        
        temp = np.zeros((len(xPosArr), 2), dtype = np.float32)
        temp[:,0] = store[:,10]
        temp[:,1] = store[:,11]
        
        
        temp1 = w.wcs_world2pix(temp, 0)
        xFrameLoc = temp1[:,0]
        yFrameLoc = temp1[:,1]
        #del temp, temp1 
        #Now do the non stars 
        for j in range(len(xPosArr)):
        #for j in [6000, 10415, 11000]:
            x = int(xFrameLoc[j])
            y = int(yFrameLoc[j])
            if(x<0 or y<0 or y>imgsizey or x>imgsizex):
                continue
            
            sigxx_temp = sigxx_Map[int(y/5), int(x/5)]
            sigxy_temp = sigxy_Map[int(y/5), int(x/5)]
            sigyy_temp = sigyy_Map[int(y/5), int(x/5)]
            
            gsizex = store[j,7] + sigxx_temp
            gsizey =  store[j,8] + sigyy_temp
            gsizexy = store[j,9] + sigxy_temp
            if((gsizex + gsizey) < 0 ):
                continue
            gsize = np.sqrt(gsizex + gsizey)
            if(gsize < 6.5):
                gsize = 6.5
            cut = img[y - int(2.5*gsize): y+int(2.5*gsize), x-int(2.5*gsize):x+int(2.5*gsize)]
            
            
            #Check if the cut is usable
            flag = 0
            flag = helper.detectBad(cut)
            if(flag == 1):
                
                if(j == 6000):
                    hdu = fits.PrimaryHDU(cut)
                    hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj1/bad' + str(file), overwrite=True)
                if(j == 10415):
                    hdu = fits.PrimaryHDU(cut)
                    hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj2/bad' + str(file), overwrite=True)
                    
                if(j == 11000):
                    hdu = fits.PrimaryHDU(cut)
                    hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj3/bad' + str(file), overwrite=True)
                    
                
                continue 
# =============================================================================
#             if(flag ==0 ):
#                 sys.exit()
#             if(len(np.where(cut ==0) [0]) == 0):
#                 sys.exit()
#                 
#             print (len(np.where(cut == 0) [0]))
# =============================================================================
            flux, mux, muy, e1, e2, back, size, sigxx, sigyy, sigxy = measure_pythonV.measure(cut)
            if(flux == None or sigxx == None or sigxy == None or mux == None):
                gmux = xFrameLoc[j] - x
                gmuy = yFrameLoc[j] - y
                flux, mux, muy, e1, e2, back, size, sigxx, sigyy, sigxy = measure_pythonV.measureSingleIter(cut, np.sqrt(2*gsizex), np.sqrt(2*gsizey), 2*gsizexy, gmux, gmuy)
                if(flux == None or sigxx == None or sigxy == None or mux == None):
                    flux, mux, muy, e1, e2, back, size, sigxx, sigyy, sigxy = measure_pythonV.measureSingleIter(cut, np.sqrt(gsizex), np.sqrt(gsizey), gsizexy, gmux, gmuy)
                    
                if(flux == None or sigxx == None or sigxy == None or mux == None):
                    if(j == 6000):
                        hdu = fits.PrimaryHDU(cut)
                        hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj1/neg' + str(file), overwrite=True)
                    if(j == 10415):
                        hdu = fits.PrimaryHDU(cut)
                        hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj2/neg' + str(file), overwrite=True)
                        
                    if(j == 11000):
                        hdu = fits.PrimaryHDU(cut)
                        hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj3/neg' + str(file), overwrite=True)
                        
                    #storei[j,frameCount,:] = gsizex, gsizey, gsizexy, gmux, gmuy, 0, 0,0, 0, 0,0, 0 , 0, 0, fwhm, zp, mjd, airmass, mphase, mAngle, 0,0,0,0,0
                    continue
            
            
            storei[j,frameCount,:] = flux, mux, muy, e1, e2, back, size, sigxx-sigxx_temp, sigyy - sigyy_temp, sigxy-sigxy_temp, raArr[j], decArr[j] , back_h, seeing, fwhm, zp, mjd, airmass, mphase, mAngle, 0,0,0,0,0
            if(j == 6000):
                hdu = fits.PrimaryHDU(cut)
                hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj1/' + str(file), overwrite=True)
            if(j == 10415):
                hdu = fits.PrimaryHDU(cut)
                hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj2/' + str(file), overwrite=True)
                
            if(j == 11000):
                hdu = fits.PrimaryHDU(cut)
                hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj3/' + str(file), overwrite=True)
            
            
            del cut
        frameCount += 1
# =============================================================================
#         if(frameCount>= 1):
#             sys.exit()
# =============================================================================
np.save( '/home/dutta26/allObj_params3z.npy', storei)
#del img                
            