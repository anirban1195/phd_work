#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 13:44:42 2022

@author: dutta26
"""
from astropy.io import fits
import numpy as np
from astropy.io import fits
import numpy as np
import helper, helper1
import measure_pythonV
from astropy.stats import sigma_clipped_stats

f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/star_field.cat')
cat = f.readlines()
f.close()
raArr =[]
decArr=[]
fluxArr=[]
sex_xx=[]
sex_yy=[]
sex_xy=[]
xList=[]
yList=[]
for j in range(len(cat)):
    if((cat[j].split()[0]) == '#'):
     continue
    
    
    raArr.append(float(cat[j].split()[5])) 
    decArr.append(float(cat[j].split()[6]))
    fluxArr.append(float(cat[j].split()[1]))
    sex_xx.append(float(cat[j].split()[7]))
    sex_yy.append(float(cat[j].split()[8]))
    sex_xy.append(float(cat[j].split()[9]))
    xList.append(float(cat[j].split()[3]))
    yList.append(float(cat[j].split()[4]))

coadd_file='/scratch/halstead/d/dutta26/abell_2390/star_field.fits'   
f=fits.open(coadd_file)
data = np.array(f[0].data)
f.close()
ySize,xSize = np.shape(data)
count = 10
store = np.zeros((len(xList), 50), dtype = np.float32)
#Fist find sources and measure params 
for j in range(len(xList)):
    #print (j)
    x = int(round( xList[j]))
    y = int(round(yList[j]))
    size = np.sqrt(sex_xx[j] +sex_yy[j])
    if(size<4):
        size = 4
    
    if(y-int(5*size) < 0 or x-int(5*size)<0 or y+int(5*size)> ySize or x+int(5*size)> xSize):
        continue
    
    cut = data[y-int(5*size): y+int(5*size), x-int(5*size): x+int(5*size)]
    
    
    flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
    
    #if Measurement failed
    if(flux == None or np.isnan(flux)):
        print (j, flux)  
        count += 1
        continue
    else:
        store[j][3:12] = flux, mux, muy, bkg, sigxx, sigyy, sigxy, x, y
        continue
    
#Find all good stars in the frame
star_arr = store[(np.where(store[:,3] > 1000))[0],  : ]
#Now tuse k sigma clip to find usable stars. Just do for sigxx
mean,median, std = sigma_clipped_stats(star_arr[:, 7])
print (mean,median, std)

star_arr = store[(np.where((store[:,7] >= mean-5*std) &
                                      (store[:,7] <= mean+5*std) & (store[:,3] > 1000) ))[0],  : ]

q,r = np.shape(star_arr)
star_temp = np.zeros(( q , 6)   , dtype = np.float32)
star_temp[:,0] = star_arr[:, 7]
star_temp[:,1] = star_arr[:, 8]
star_temp[:,2] = star_arr[:, 9]
star_temp[:,3] = star_arr[:, 10]
star_temp[:,4] = star_arr[:, 11]       
nStars = 10


for j in range(len(xList)):
    #print (j)
    x = int(round( xList[j]))
    y = int(round(yList[j]))
    
    temp = np.copy(star_temp)
    temp[:,5] = ((temp[:,3]-x)**2 + (temp[:,4]-y)**2)**0.5 
    temp = temp[temp[:,5].argsort()]
    #if(j==1000):
    #    sys.exit()
    #Check if same star. Then delete the entry
    if(temp[0,5]<5):
        print (store[j, 2])
        temp = np.delete(temp, (0), axis = 0)
    
    #Checking for nans to avoid code from crashing
    if(np.isnan(np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
        continue
    if(np.isnan(np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
        continue
    if(np.isnan(np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
        continue
    
    avgSigxx = np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
    avgSigyy = np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
    avgSigxy = np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
    
    #avgSigxx = np.mean( temp[0:nStars, 0])
    #avgSigyy = np.mean( temp[0:nStars, 1])
    #avgSigxy = np.mean( temp[0:nStars, 2])   
                      
                      
    store[j,38] = avgSigxx
    store[j,39] = avgSigyy
    store[j,40] = avgSigxy
    store[j,41] = np.std(temp[0:nStars, 0])
    store[j,42] = np.std(temp[0:nStars, 1])
    store[j,43] = np.std(temp[0:nStars, 2])
    
    del temp
loc= (np.where( store[:,3] > 1000 ))[0]    
print(sigma_clipped_stats(store[loc,7] - store[loc,38]))
print(sigma_clipped_stats(store[loc,7]))
