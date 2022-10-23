#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 10:20:18 2022

@author: dutta26
"""


from astropy.io import fits
import numpy as np
import helper, helper1
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import os,shutil,sys

raList_star =[]
decList_star =[]
Raindex=5
Decindex=6
elonIndex = 12
xxIndex = 7
yyIndex = 8
xyIndex = 9
coadd_file = '/scratch/halstead/d/dutta26/abell_2390/sample1_sf.fits'
coadd_df = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data.npy')
for chipFile in os.listdir('/scratch/halstead/d/dutta26/abell_2390/chip_cuts/'):
        os.remove('/scratch/halstead/d/dutta26/abell_2390/chip_cuts/'+chipFile)


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
f=fits.open(coadd_file)
data = np.array(f[0].data)
f.close()

raList = np.array(raList)
decList = np.array(decList)


min_dist_arr = np.load('/scratch/halstead/d/dutta26/abell_2390/min_dist.npy')
count = 0
store[:, 13] = coadd_df[:,13]
store[:, 2] = coadd_df[:, 2]

        
        
            
#np.save('/scratch/halstead/d/dutta26/abell_2390/wiyn_data.npy', store)

#Now measure 
#coadd_file = '/scratch/halstead/d/dutta26/abell_2390/sample1.fits'
#store = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data.npy')
xList, yList = helper.convertToXY(raList, decList, coadd_file)
f=fits.open(coadd_file)
data = np.array(f[0].data)
f.close()
ySize,xSize = np.shape(data)
count = 10
chipCount = 0
fileLoc='/scratch/halstead/d/dutta26/abell_2390/odi_img/'
for files in os.listdir(fileLoc):
    f=fits.open(fileLoc+files)
    data = f[0].data
    f.close()
    chipCount += 1
    xList, yList = helper.convertToXY(raList, decList, fileLoc+files)
    xMax, yMax = np.shape(data)
    #Fist find sources and measure params 
    for j in range(len(xList)):
        #print (j)
        
        ra = raList[j]
        dec = decList[j]
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        if(x<28 or y<28 or x>xMax-28 or y>yMax-28):
            continue
        size = np.sqrt(sex_xx[j] +sex_yy[j])
        if(size<5):
            size = 5
        
        
        cut = data[y-int(5*size): y+int(5*size), x-int(5*size): x+int(5*size)]
        
        zeroLoc= np.where(cut == 0)[0]
        if(len(zeroLoc)> 0):
            continue
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
        
        #if Measurement failed
        if(flux == None or np.isnan(flux)):
            #print (j, flux)  
            count += 1
            continue
        else:
            store[j][3:12] = flux, mux, muy, bkg, sigxx, sigyy, sigxy, x, y
            store[j][44] = chipCount
            continue
    

    #Now do a second pass to find interpolated values 
    #Find all good stars in the frame
    loc = (np.where((store[:,2] == 1)& (store[:,3] > 100)& (store[:,3] < 5000) & (store[:,13] == 1)& (store[:,44] == chipCount) ))[0]
    if(len(loc) ==0):
        continue
    star_arr = store[loc,  : ]
    #Now tuse k sigma clip to find usable stars. Just do for sigxx
    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    print (mean,median, std)
    loc =(np.where((store[:,2] == 1) & 
                                          (store[:,7] >= mean-5*std) &
                                          (store[:,7] <= mean+5*std) &
                                          (store[:,3] > 50000)  &
                                          (store[:,3] < 100000) & (store[:,13] == 1)& (store[:,44] == chipCount) ))[0]
    if(len(loc) == 0):
        continue
    star_arr = store[loc,  : ]
    
    
        
        
    q,r = np.shape(star_arr)
    star_temp = np.zeros(( q , 6)   , dtype = np.float32)
    star_temp[:,0] = star_arr[:, 7]
    star_temp[:,1] = star_arr[:, 8]
    star_temp[:,2] = star_arr[:, 9]
    star_temp[:,3] = star_arr[:, 10]
    star_temp[:,4] = star_arr[:, 11]       
    nStars = 10
    
    
    
    shutil.copy(fileLoc+files,'/scratch/halstead/d/dutta26/abell_2390/chip_cuts/'+str(chipCount)+'.fits')
    f_test=fits.open('/scratch/halstead/d/dutta26/abell_2390/chip_cuts/'+str(chipCount)+'.fits', mode='update')

    for j in range(len(star_arr)):
        x=int(star_arr[j,10])
        y=int(star_arr[j,11])
        f_test[0].data[y-20, x-20:x+20] = 1000
        f_test[0].data[y+20, x-20:x+20] = 1000
        f_test[0].data[y-20:y+20, x-20] = 1000
        f_test[0].data[y-20:y+20, x+20] = 1000
        
    f_test.flush()    
    
    print (max(star_arr[:, 10]), max(star_arr[:, 11]))
    for j in range(len(xList)):
        #print (j)
        ra = raList[j]
        dec = decList[j]
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        
        temp = np.copy(star_temp)
        temp[:,5] = ((temp[:,3]-x)**2 + (temp[:,4]-y)**2)**0.5 
        temp = temp[temp[:,5].argsort()]
        #if(j==1000):
        #    sys.exit()
        #Check if same star. Then delete the entry
        if(temp[0,5]<5):
            #print (store[j, 2])
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

#sys.exit()    
np.save('/scratch/halstead/d/dutta26/abell_2390/wiyn_data_chip50k100k.npy', store)
loc= (np.where((store[:,2] == 1) & (store[:,3] > 2000)  &(store[:,3] < 40000) & (store[:,13] == 1) ))[0]
print(sigma_clipped_stats(store[loc,7] - store[loc,38]))
print(sigma_clipped_stats(store[loc,7]))
print (np.shape(star_arr))




# =============================================================================
# import numpy as np
# from scipy.interpolate import griddata
# 
# 
# 
# points=[]
# values =[] 
# grid_y, grid_x = np.meshgrid(np.linspace(0, 1000, 1000), np.linspace(0, 1000, 1000))
# 
# for index in loc:
#     x = int(round(store[index,10]))
#     y = int(round(store[index,11]))
#     
#     points.append([y,x])
#     values.append( store[index,7])
# 
# 
# grid_z1 = griddata(points, values, (grid_y, grid_x), method='cubic')    
# grid_z1 = np.array(grid_z1, dtype = np.float32).T
# hdu = fits.PrimaryHDU(grid_z1) 
# hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/residue.fits', clobber=True)  
# 
# 
# x = np.linspace(1, 1000, 1000)
# y = np.linspace(0, 1000, 1000)
# X, Y = np.meshgrid(x, y, copy=False)
# Z = grid_z1
# 
# X = X.flatten()
# Y = Y.flatten()
# 
# A = np.array([X*0+1, X, Y, X**2, X**2*Y, X**2*Y**2, Y**2, X*Y**2, X*Y]).T
# B = Z.flatten()
# 
# coeff, r, rank, s = np.linalg.lstsq(A, B)
# 
# 
# =============================================================================





   
# =============================================================================
# store= np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data_final.npy')
# a = store[:,7]-store[:,38]
# b = store[:,8]-store[:,39]
# e1 = (store[:,7]-store[:,8])/(store[:,7]+store[:,8])
# e2 = (2*store[:,9])/(store[:,7]+store[:,8])
# ellip = np.sqrt(e1**2 + e2**2)
# 
# loc = np.where((store[:,2]==0 ) & ((a<0)|(b<0)) & (store[:,13]==0 ) & (store[:,3]<= 2000 ))[0]
# loc_new= []
# for index in loc:
#     print (index)
#     x = int(store[index,10])
#     y = int(store[index,11])
#     flag =0 
#     
#     if(min_dist_arr[index]< 40):
#         flag = 1
# # =============================================================================
# #     maxVal = index+3000
# #     minVal = index-3000
# #     if(maxVal > len(a)):
# #         maxVal = len(a)
# #     if(minVal < 0):
# #         minVal = 0
# #     
# #     for j in range(minVal, maxVal):
# #         if(j == index):
# #             continue
# #         x1 = store[j, 10]
# #         y1 = store[j, 11]
# #         if(np.sqrt((x-x1)**2 +(y-y1)**2) < 40):
# #             flag = 1
# # =============================================================================
#     if(flag==0):
#         loc_new.append(index)
#        
# 
# np.save('/scratch/halstead/d/dutta26/abell_2390/loc_new.npy', loc_new)
# =============================================================================

# =============================================================================
# loc_new = np.load('/scratch/halstead/d/dutta26/abell_2390/loc_new.npy')
# store= np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data_final.npy')
# f=fits.open('/scratch/halstead/d/dutta26/abell_2390/sample_long1.fits')
# long_data = f[0].data
# f.close()
# 
# e1_diff =[]
# e2_diff=[]
# fluxArr=[]
# xList, yList = helper.convertToXY(store[loc_new,0], store[loc_new,1], '/scratch/halstead/d/dutta26/abell_2390/sample_long1.fits') 
# #f_test=fits.open('/scratch/halstead/d/dutta26/abell_2390/test1.fits', mode='update')
# j=0
# for index in loc_new:
#     #ra = store[index,0]
#     #dec = store[index,1]
#     #x = int(store[index,10])
#     #y = int(store[index,11])
#     print (index)
#        
#     x=int(xList[j])
#     y=int(yList[j])
#     j=j+1
#     cut = long_data[y-20:y+20, x-20:x+20]
#     flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
#     if(flux == None or np.isnan(flux)):
#         continue
#     #xx_diff.append(sigxx - store[index,7])
#     e1 = (store[index,7] - store[index,8])/(store[index,7] + store[index,8])
#     e2 = 2* store[index,9]/(store[index,7] + store[index,8])
#     
#     e1_long = (sigxx-sigyy)/(sigxx+sigyy)
#     e2_long = 2*sigxy / (sigxx+sigyy)
#     e1_diff.append(e1-e1_long)
#     e2_diff.append(e2-e2_long)
#     fluxArr.append(store[index,3])
#     #f_test[0].data[y-20, x-20:x+20] = 1000
#     #f_test[0].data[y+20, x-20:x+20] = 1000
#     #f_test[0].data[y-20:y+20, x-20] = 1000
#     #f_test[0].data[y-20:y+20, x+20] = 1000
#     # 
# #f_test.flush()
# 
#     
#     
#     
# # =============================================================================
# # f_test=fits.open('/scratch/halstead/d/dutta26/abell_2390/test1.fits', mode='update')
# # for index in loc_new:
# #     x = int(store[index,10])
# #     y = int(store[index,11])
# #     f_test[0].data[y-20, x-20:x+20] = 1000
# #     f_test[0].data[y+20, x-20:x+20] = 1000
# #     f_test[0].data[y-20:y+20, x-20] = 1000
# #     f_test[0].data[y-20:y+20, x+20] = 1000
# #      
# # f_test.flush()
# # =============================================================================
# 
# 
# =============================================================================
