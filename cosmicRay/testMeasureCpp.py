#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 10:06:17 2020

@author: anirban
"""

import numpy as np
from astropy.io import fits
import wrapper
from scipy.ndimage import rotate
import pandas as pd
filename = '/scratch/halstead/d/dutta26/gauss.fits'
outputText = '/scratch/halstead/d/dutta26/results.csv'
# =============================================================================
# outputText = '/home/anirban/results.csv'
# sigmaxArr= np.hstack(np.arange(0.5 , 1.5, 0.1), np.arange(2,5,1))
# sigmayArr = np.hstack(np.arange(0.5 , 1.5, 0.1), np.arange(2,5,1))
# signalArr = np.hstack(np.arange(100,1000, 100), np.arange(1000,10000, 1000))
# rotArr = np.arange(0,45,5)
# =============================================================================
sigmaxArr= [1]
sigmayArr = [1.2]
signalArr = [1000]
rotArr = [0.0]
meanFluxErr =[]
stdFluxErr =[]
meanSigxErr =[]
stdSigxErr =[]
meanSigyErr =[]
stdSigyErr =[]
meanEccenErr =[]
stdEccenErr=[]
realSigxArr=[]
realSigyArr=[]
realFluxArr=[]
realRotArr =[]
failRateArr =[]
for sigma_x in sigmaxArr:
    for sigma_y in sigmayArr:
        for flux in signalArr:
            for rot in rotArr:
                fluxErrArr=[]
                sigmaxErrArr = []
                sigmayErrArr =[]
                eccenErrArr =[]
                rate = 0 
                for j in range(100):
                    mux = 20
                    muy = 20
                    x = np.linspace(0, 49, 50)
                    y = np.linspace(0, 49, 50)
                    x, y = np.meshgrid(x, y)
                    z = (1 * np.exp(-((x-mux)**2/(2*sigma_x**2)+ (y-muy)**2/(2*sigma_y**2))))
                    z=(z*flux)/np.sum(z)
                    z = rotate(z,-rot, reshape=False)
                    noise = np.random.normal(100, 10, (50,50))
                    z=z+ noise
                    hdu = fits.PrimaryHDU(z)  
                    hdu.writeto(filename, clobber=True)
                    #Run detection algorithm 
                    filename_encode = filename.encode('utf-8')
                    a=wrapper.measure([b"abc", filename_encode])
                    #Measurement done. Now extract values
                    f=open('/scratch/halstead/d/dutta26/test.txt')
                    content = f.readlines()
                    f.close()
                    if(len(content) == 0):
                        rate = rate + 1
                        continue
                    sigmaxx = float(content[len(content)-1].split()[7])
                    sigmayy = float(content[len(content)-1].split()[8])
                    sigmaxy = float(content[len(content)-1].split()[9])
                    flux_measure =float(content[len(content)-1].split()[3])
                    if(sigmaxx < 0 or sigmayy< 0 or (4*(sigmaxy**2) + (sigmaxx-sigmayy)**2) < 0 ):
                        rate =rate +1
                        continue
                    #Extraction done
                    #Now validate
                    minSize = ((sigmaxx+sigmayy)*0.5) - 0.5*np.sqrt(4*(sigmaxy**2) + (sigmaxx-sigmayy)**2)
                    maxSize = ((sigmaxx+sigmayy)*0.5) + 0.5*np.sqrt(4*(sigmaxy**2) + (sigmaxx-sigmayy)**2)
                    if(minSize<0 or maxSize<0):
                        rate =rate +1
                        continue
                    sigmamin = np.sqrt(2*minSize)
                    sigmamax = np.sqrt(2*maxSize)
                    
                    eccen_measured = np.sqrt(1-sigmamin/sigmamax)
                    if(sigma_x > sigma_y):
                        eccen_real = np.sqrt(1-sigma_y/sigma_x)
                        eccenErrArr.append(abs(eccen_measured-eccen_real)/eccen_real)
                        sigmaxErrArr.append(abs(sigma_x- sigmamax)/sigma_x )
                        sigmayErrArr.append(abs(sigma_y- sigmamin)/sigma_y )
                        
                    else:
                        eccen_real = np.sqrt(1-sigma_x/sigma_y)
                        eccenErrArr.append(abs(eccen_measured-eccen_real)/eccen_real)
                        sigmaxErrArr.append(abs(sigma_y- sigmamax)/sigma_y )
                        sigmayErrArr.append(abs(sigma_x- sigmamin)/sigma_x )
                        
                    fluxErrArr.append(abs(flux_measure-flux)/flux)
                if(rate > 90):
                    failRateArr.append(1)
                    continue
                failRateArr.append(rate/100.0)
                meanFluxErr.append(np.mean(fluxErrArr))
                stdFluxErr.append(np.std(fluxErrArr))
                meanSigxErr.append(np.mean(sigmaxErrArr))
                stdSigxErr.append(np.std(sigmaxErrArr))
                meanSigyErr.append(np.mean(sigmayErrArr))
                stdSigyErr.append(np.std(sigmayErrArr))
                meanEccenErr.append(np.mean(eccenErrArr))
                stdEccenErr.append(np.std(eccenErrArr))
                realSigxArr.append(max(sigma_x, sigma_y))
                realSigyArr.append(min(sigma_x, sigma_y))
                realFluxArr.append(flux)
                realRotArr.append(rot)
# =============================================================================
# 
# df = pd.DataFrame(np.array(failRateArr, meanFluxErr, stdFluxErr,meanSigxErr ,stdSigxErr, meanSigyErr ,
#        stdSigyErr,meanEccenErr, stdEccenErr,realSigxArr, realSigyArr,realFluxArr,realRotArr),
#     columns=('Fail rate', 'Mean Flux Err', 'Std Dev of Flux Err', 'Mean Sigmax Err', 'Std Dev of Sigmax Err',
#              'Mean Sigmay Err', 'Std Dev of Sigmay Err', 'Mean Eccentricity Err', 'Std Dev of Eccentricity Err',
#              'Real Sigmax', 'Real Simgay', 'Real Flux', 'Real Rot'))
# =============================================================================
    

                     
d = {'Fail rate': failRateArr, 'Mean Flux Err': meanFluxErr, 'Std Dev of Flux Err': stdFluxErr ,
     'Mean Sigmax Err': meanSigxErr, 'Std Dev of Sigmax Err': stdSigxErr,
     'Mean Sigmay Err': meanSigyErr, 'Std Dev of Sigmay Err': stdSigyErr,
     'Mean Eccentricity Err': meanEccenErr,'Std Dev of Eccentricity Err':stdEccenErr,
     'Real Sigmax': realSigxArr,'Real Sigmay': realSigyArr, 'Real Flux': realFluxArr,
     'Real Rot': realRotArr}        

df = pd.DataFrame(data=d)            
df.to_csv(outputText, index=False)                     
                    
                    

                    
                