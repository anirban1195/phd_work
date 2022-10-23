#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  7 11:36:33 2022

@author: dutta26
"""


from astropy.io import fits
import numpy as np
import helper, helper1
import measure_pythonV
from astropy.stats import sigma_clipped_stats


raList_star =[]
decList_star =[]
Raindex=5
Decindex=6
elonIndex = 12
xxIndex = 7
yyIndex = 8
xyIndex = 9

coadd_file = '/scratch/halstead/d/dutta26/abell_2390/sample1.fits'
#store = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data.npy')

f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/test_stars.cat')
content = f.readlines()
f.close()
#Create a list of stars with ra and dec 
for j in range(len(content)):
    #Skip the comment lines
    if((content[j].split()[0]) == '#'):
        continue
    #Reject cosmic rays
    if(float(content[j].split()[elonIndex]) > 1.5):
        continue
    
    raList_star.append(float(content[j].split()[Raindex])) 
    decList_star.append(float(content[j].split()[Decindex])) 
    
    
f_test=fits.open('/scratch/halstead/d/dutta26/abell_2390/test.fits', mode='update')
f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/test_all.cat')
content = f.readlines()
f.close()
raList =[]
decList=[]
sex_xx=[]
sex_yy=[]
sex_xy =[]
elonList=[]
#Create a list of stars with ra and dec 
for j in range(len(content)):
    #Skip the comment lines
    if((content[j].split()[0]) == '#'):
        continue

    raList.append(float(content[j].split()[Raindex])) 
    decList.append(float(content[j].split()[Decindex])) 
    sex_xx.append(float(content[j].split()[xxIndex]))
    sex_yy.append(float(content[j].split()[yyIndex]))
    sex_xy.append(float(content[j].split()[xyIndex]))
    elonList.append(float(content[j].split()[elonIndex]))

xList, yList = helper.convertToXY(raList, decList, coadd_file)    
store = np.zeros((len(raList), 50), dtype = np.float32)
#Designate stars as 1
for j in range(len(raList)):
    print (j)
    if(elonList[j] > 1.5):
        continue
    ra = raList[j]
    dec = decList[j]
    store[j, 0:2] = [ra, dec]
    #minVal = int(j/2) - 10000
    #maxVal = int(j/2) + 10000
    minVal = 0
    maxVal = len(raList_star)
    if(minVal< 0):
        minVal = 0
    if(maxVal > len(raList_star) ):
        maxVal = len(raList_star)
    for k in range(minVal, maxVal):
        ra_star = raList_star[k]
        dec_star = decList_star[k]
        
        if(np.sqrt( (ra-ra_star)**2 + (dec-dec_star)**2 ) < 0.0002 ):
            store[j , 2] = 1
            x= int(xList[j])
            y= int(yList[j])
            if(x<22):
                x=22
            if(y<22):
                y=22
            if(y>9970):
                y=9970
            if(x>9970):
                x=9970
            f_test[0].data[y-20, x-20:x+20] = 1000
            f_test[0].data[y+20, x-20:x+20] = 1000
            f_test[0].data[y-20:y+20, x-20] = 1000
            f_test[0].data[y-20:y+20, x+20] = 1000
            
np.save('/scratch/halstead/d/dutta26/abell_2390/wiyn_data.npy', store)

f_test.flush()

#########################################################################################################

f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/test_all.cat')
content = f.readlines()
f.close()
raList =[]
decList=[]
sex_xx=[]
sex_yy=[]
sex_xy =[]
#Create a list of stars with ra and dec 
for j in range(len(content)):
    #Skip the comment lines
    if((content[j].split()[0]) == '#'):
        continue
    
    raList.append(float(content[j].split()[Raindex])) 
    decList.append(float(content[j].split()[Decindex])) 
    sex_xx.append(float(content[j].split()[xxIndex]))
    sex_yy.append(float(content[j].split()[yyIndex]))
    sex_xy.append(float(content[j].split()[xyIndex]))
    
    
coadd_file = '/scratch/halstead/d/dutta26/abell_2390/sample1.fits'
store = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data.npy')
xList, yList = helper.convertToXY(raList, decList, coadd_file)
f=fits.open(coadd_file)
data = np.array(f[0].data)
f.close()
ySize,xSize = np.shape(data)
count = 0

#Fist find stars and measure params 
for j in range(len(xList)):
    #print (j)
    
    ra = raList[j]
    dec = decList[j]
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
# =============================================================================
#         hdu = fits.PrimaryHDU(cut)
#         hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj1/'+str(j)+'.fits')
# =============================================================================
        continue
    else:
        store[j][3:12] = flux, mux, muy, bkg, sigxx, sigyy, sigxy, x, y
        continue
    

#Now do a second pass to find interpolated values 
#Find all good stars in the frame
star_arr = store[(np.where((store[:,2] == 1) ))[0],  : ]
#Now tuse k sigma clip to find usable stars. Just do for sigxx
mean,median, std = sigma_clipped_stats(star_arr[:, 7])
print (mean,median, std)

star_arr = store[(np.where((store[:,2] == 1) & 
                                      (store[:,7] >= mean-5*std) &
                                      (store[:,7] <= mean+5*std) &
                                      (store[:,12] == 0) & 
                                      (store[:,13] == 0) & 
                                      (store[:,14] == 0)))[0],  : ]

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
    ra = raList[j]
    dec = decList[j]
    x = int(round( xList[j]))
    y = int(round(yList[j]))
    
    temp = np.copy(star_temp)
    temp[:,5] = (temp[:,3]-x)**2 + (temp[:,4]-y)**2
    temp = temp[temp[:,5].argsort()]
    
    #Check if same star. Then delete the entry
    if(temp[0,5]<5):
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
    
    
    
    store[j,38] = avgSigxx
    store[j,39] = avgSigyy
    store[j,40] = avgSigxy
    store[j,41] = np.std(temp[0:nStars, 0])
    store[j,42] = np.std(temp[0:nStars, 1])
    store[j,43] = np.std(temp[0:nStars, 2])
    
    del temp
    
np.save('/scratch/halstead/d/dutta26/abell_2390/wiyn_data_final.npy', store)


# =============================================================================
# store = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data_final.npy')
# loc = np.where((store[:,2] == 0 ) &(store[:,3] > 0 ) )[0]
# f=fits.open('/scratch/halstead/d/dutta26/abell_2390/test.fits', mode='update')
# for index in loc:
#     x = int(store[index, 10])
#     y = int(store[index, 11])
#     if(x>9970 or y>9970):
#         continue
#     if(x<30 or y<30):
#         continue
#     f[0].data[y-20, x-20:x+20] = 1000
#     f[0].data[y+20, x-20:x+20] = 1000
#     f[0].data[y-20:y+20, x-20] = 1000
#     f[0].data[y-20:y+20, x+20] = 1000
#     
# f.flush()
# =============================================================================
    
store = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data_final.npy')
a,b = np.shape(store)
# =============================================================================
# flagArr=[]
# count = 0
# for j in range(a):
#     flag = 0
#     flagSig =1
#     print (j)
#     x = int(store[j, 10])
#     y= int (store[j,11])
#     
#     maxVal = j + 2000
#     minVal = j- 2000
#     if(minVal <0):
#         minVal=0
#     if(maxVal > a ):
#         maxVal = a
#     for k in range(minVal, maxVal):
#         if(k == j):
#             continue
#         x1 = int(store[k, 10])
#         y1= int (store[k,11] )
#         if(np.sqrt((x-x1)**2 + (y-y1)**2) < 40):
#             flag = 1
#             
#     if((store[j, 7] - store[j, 38])< 0 or (store[j, 8] - store[j, 39])< 0):
#         flagSig = 0
#             
#     
#     if(flag ==0 and flagSig == 0 and store[k, 2]==0):
#         count +=1
#         flagArr.append(1)
#     else:
#         flagArr.append(0)
# flagArr = np.array(flagArr)
# np.save('/scratch/halstead/d/dutta26/abell_2390/flagArr.npy', flagArr)
# =============================================================================

flagArr = np.load('/scratch/halstead/d/dutta26/abell_2390/flagArr.npy')
fluxDiff=[]
sigxx_diff=[]
sigyy_diff=[]
sigxy_diff=[]
f=fits.open('/scratch/halstead/d/dutta26/abell_2390/sample_long1.fits')    
data = f[0].data
f.close()
c = 0

f=fits.open('/scratch/halstead/d/dutta26/abell_2390/test1.fits', mode='update')

for j in range(a):
    if(flagArr[j] == 0):
        continue
    x = int(store[j, 10])
    y = int (store[j,11])
    f[0].data[y-20, x-20:x+20] = 1000
    f[0].data[y+20, x-20:x+20] = 1000
    f[0].data[y-20:y+20, x-20] = 1000
    f[0].data[y-20:y+20, x+20] = 1000
    
f.flush()
# =============================================================================
#     cut = data[y-20:y+20, x-20:x+20]
#     #hdu = fits.PrimaryHDU(cut)
#     #hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj1/'+str(j)+'.fits', overwrite=True) 
#     flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
#     if(flux == None or np.isnan(flux)):
#         c+= 1
#         continue
#     sigxx_diff.append(sigxx )
#     sigyy_diff.append(sigyy )
#     fluxDiff.append(flux/10 - store[j, 3])
# =============================================================================


