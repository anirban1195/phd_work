#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 22:11:19 2021

@author: dutta26
"""


import astropy.units as u
import helper
import matplotlib.pyplot as plt
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
from astroquery.gaia import Gaia
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord   
from astropy.coordinates import FK5 
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import sys,shutil
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
lowLtArr =[ [2782,17237], [2346, 16796], [10192, 16350], [10906, 16070], [10995, 19202], 
           [4103, 16424], [2909, 16090], [2398, 16104], [1974, 15839], [687, 14304], 
           [1204,14336], [1269,14755], [1702, 14331], [2204, 14378], [2697, 14372],
           [3219,14414], [3736, 14411], [3076, 13808], [16290, 12053], [13565, 12717], 
           [12983, 8748], [16453, 9099], [12986, 6917], [18688, 1333], [20149, 498], 
           [2692, 15596], [7451, 19267], [2025, 15543], [3559, 15555], [5075, 15555],
           [5583, 15529], [16303, 11273], [16471, 11012], [18689, 2819], [20507, 0]]
upLtArr = [ [3045,17476], [2591, 16879], [10434, 16413], [11195, 16103], [11141, 19471],
           [4498, 16600], [2993, 16233], [2488, 16293], [2157, 16030], [780, 14758],
           [1269,14755], [1490,14847], [1795, 14775], [2297,14698], [2795,14716],
           [3327, 14608], [3811, 14650], [3219, 13933], [16387, 12243], [13759, 12914],
           [13180, 8891], [16682, 9242], [13127, 7136], [18799,2483], [20260, 1752],
           [2877, 15672], [7515, 19392], [2100,15793], [3620, 15704], [5143, 15697],
           [5655, 15701], [16401, 11437], [16531, 11091], [18785, 3052], [23500, 1426]]


rejectLowLtArr =[ [3305,23439], [1677,22593], [1830,22744], [1744, 22551], [497,22425],
                 [130,22484], [128,21952], [1428, 22254], [8883,22543], [9577, 22331], 
                 [8121,21446], [9090, 19997], [20205, 18340], [5215, 13553],[14081, 11060],
                 [14487, 11540],[16835, 10025], [16887, 9771], [16790, 9607], [16417, 9594],
                 [14328, 9890], [16423, 10047], [2694, 10109], [14417, 8496], [13954, 7623],
                 [14663, 7324], [16253, 8173], [16418, 8203], [16685, 6853], [13383, 6472],
                 [13189, 6455], [14969, 5520], [15140, 5540], [5210, 22194], [4695, 22196],
                 [4187, 22199], [12749, 22974]]

rejectUpLtArr= [[3317,23497] , [1696, 22665], [1851, 22819], [1768, 22630], [520, 22504], 
                [155,22576], [144, 22020], [1457, 22331], [8895,22968], [9597, 22655],
                [8131,21764], [9100, 20261], [20245, 18480], [5233, 13650], [14161, 11080],
                [14564,11560], [16909, 10092], [16909,9972], [16902, 9763], [16563, 9681], 
                [14410, 9918], [16468, 10104], [2714, 10166], [14631, 8551], [14170,7680],
                [14865, 7381], [16331, 8195], [16492, 8222], [16778, 6878], [13453, 6495], 
                [13264, 6477], [15043, 5540], [15223, 5557], [5240, 22303], [4727, 22301],
                [4212,22299 ] , [12789, 23088]] 






plotLoc = '/scratch/bell/dutta26/abell_2390/plot/'
Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" # Select early Data Release 3
# =============================================================================
# catalog = str(sys.argv[1])
# outFile = str(sys.argv[2])
# imgFile = str(sys.argv[3])
# =============================================================================
catalog_low = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/abell_2390_sextract_low.cat'
catalog_high = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/abell_2390_sextract_high.cat'

outFile = '/home/dutta26/codes/source_list.pk1'
imgFile = '/scratch/bell/dutta26/abell_2390/abell_ir_coadd_crop.fits'
wted_img_file = '/scratch/bell/dutta26/abell_2390/abell_ir_coadd.fits'
fluxMin = 10
psfMin = 2.5
psfMax = 3.4
totScale = 34244.70
totBkg = 75517.35

f=fits.open(imgFile)
img = f[0].data
sizey, sizex = np.shape(img)

#Find the max, min and mid Ra Dec
w=wcs.WCS(f[0].header)
f.close()
[[startRa, startDec]]  = w.wcs_pix2world([[0,0]], 0)
[[endRa, endDec]] = w.wcs_pix2world([[sizex, sizey]], 0)
midRa = 0.5*(startRa +endRa )
midDec = 0.5* (startDec + endDec)
#Make an image of same size as input and number code it. 0= Normal, 1= Use high arr, 2= Reject
coded_img = np.zeros((sizey,sizex), dtype = np.int16)
for j in range(len(lowLtArr )):
    [xLow , yLow ]= lowLtArr[j]
    [xUp , yUp ]= upLtArr[j]
    coded_img[yLow:yUp, xLow:xUp] = 1
for j in range(len(rejectLowLtArr )):
    [xLow , yLow ]= rejectLowLtArr[j]
    [xUp , yUp ]= rejectUpLtArr[j]
    coded_img[yLow:yUp, xLow:xUp] = 2


#Get stars from GAIA within a degree of mid ra dec
coord = SkyCoord(ra=midRa , dec=midDec, unit=(u.degree, u.degree), frame='icrs')
width = u.Quantity(np.abs(startRa - endRa ), u.deg)
height = u.Quantity(np.abs(startDec - endDec), u.deg)
Gaia.ROW_LIMIT = -1
r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
r.pprint()

df = r.to_pandas()

#sys.exit()

#df.drop(df[df['parallax_over_error'] < 0.5].index, inplace = True)
df.drop(df[df['phot_g_mean_mag'] < 14].index, inplace = True)

#df1 = df[['ra', 'dec', 'parallax_over_error']]
df2 = df[df['parallax_over_error'].notna()]

#Read the catalog 
with open(catalog_low) as f:
    content_low = f.readlines()
    
with open(catalog_high) as f:
    content_high = f.readlines()

xIndex = 3
yIndex = 4
Raindex=5
Decindex=6
xxIndex = 7
yyIndex = 8
xyIndex = 9

#Make an array containing x and y indices 
raList =[]
decList =[]
sex_xx_list = []
sex_yy_list = []
sex_xy_list = []
star_list =[]
for j in range(len(content_low)):
    #Skip the comment lines
    if((content_low[j].split()[0]) == '#'):
        continue
    x = float(content_low[j].split()[xIndex])-1
    y = float(content_low[j].split()[yIndex])-1
    x=int(round(x))
    y=int(round(y))
    if(coded_img[y,x] == 0):
        raList.append(float(content_low[j].split()[Raindex])) 
        decList.append(float(content_low[j].split()[Decindex])) 
        sex_xx_list.append(float(content_low[j].split()[xxIndex]))
        sex_yy_list.append(float(content_low[j].split()[yyIndex]))
        sex_xy_list.append(float(content_low[j].split()[xyIndex]))
        star_list.append(0)
        
for j in range(len(content_high)):
    #Skip the comment lines
    if((content_high[j].split()[0]) == '#'):
        continue
    x = float(content_high[j].split()[xIndex])-1
    y = float(content_high[j].split()[yIndex])-1
    x=int(round(x))
    y=int(round(y))
    if(coded_img[y,x] == 1):
        raList.append(float(content_high[j].split()[Raindex])) 
        decList.append(float(content_high[j].split()[Decindex])) 
        sex_xx_list.append(float(content_high[j].split()[xxIndex]))
        sex_yy_list.append(float(content_high[j].split()[yyIndex]))
        sex_xy_list.append(float(content_high[j].split()[xyIndex]))
        star_list.append(0)

#Make a pandas dataframe for the sources . 1 for star and 0 for galaxy
a=np.array([raList, decList, sex_xx_list, sex_yy_list, sex_xy_list, star_list])
b=np.swapaxes(a, 0,1)
df_source = pd.DataFrame(b,  columns = ['ra', 'dec', 'sex_xx', 'sex_yy', 'sex_xy', 'star_bool'])

    




finalMatchIndexArr =[]
for j in range(len(df2)):
    
#for j in [2947, 2954]:
    ra = float( df2['ra'][j:j+1])
    dec = float( df2['dec'][j:j+1])
    
    fluxMatch=0
    matchFlag = 0 
    matchIndexArr=[]
    distArr =[]
    #print ('$$$$$$$$$$$$$$$$$$')
    for k in range(len(raList)):
        if(np.abs(raList[k] - ra) < 0.003 and np.abs(decList[k] - dec) < 0.003):
            matchIndexArr.append(k)
            distArr.append(np.sqrt((raList[k] - ra)**2 + (decList[k] - dec)**2))
    
    
    #If no match is found       
    if(len(matchIndexArr) == 0):
        continue
        
    #If only 1 match is found and is close 
    if(len(matchIndexArr) == 1 and distArr[0]<= 0.0002):
        #print ('aaaaaaa')
        df_source['star_bool'][matchIndexArr[0]] = 1
        continue
    #If multiple matches are found 
    temp_dist = np.sort(distArr)
    #print (temp_dist)
    #if(temp_dist[0]< 0.0002 and temp_dist[1]> 0.0008):
    if(temp_dist[0]< 0.0002 and temp_dist[1]> 0.0016):
        #print ('bbbbbbbbb')
        leastDistIndex = np.where(distArr == temp_dist[0])
        df_source['star_bool'][matchIndexArr[leastDistIndex[0][0]]] = 1



#Find additional stars not in GAIA catalog and mark them as 2
f=fits.open(wted_img_file)
img = f[0].data
f.close()
sizey, sizex= np.shape(img)
sizeArr =[]
fluxArr =[]
boolArr=[]
xList, yList = helper.convertToXY(raList, decList, wted_img_file)
for j in range(len(xList)):
    
    print (j)
    x = int(round( xList[j]))
    y = int(round(yList[j]))
    
    if(x<30 or y<30 or x>(sizex-30) or y>(sizey-30)):
        fluxArr.append(0)
        sizeArr.append(0)
        continue
    cut = img[y-25: y+25, x-25: x+25]
    flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2)
    print (flux, psf)
    if(flux == None or sigxx== None or np.isnan(psf)):
        fluxArr.append(0)
        sizeArr.append(0)
        boolArr.append(99)
        continue
    fluxArr.append(flux)
    sizeArr.append(psf)
    if(df_source['star_bool'][j] == 1):
        boolArr.append(1)
        continue
    if(np.log10(flux)>np.log10(fluxMin) and psf>psfMin and psf<psfMax):
        df_source['star_bool'][j] = 2
    if(psf>12):
        df_source['star_bool'][j] = 3
    if(flux> 1e4):
        df_source['star_bool'][j] = 4
    size_err = np.sqrt( 0.5 *( (psfMin)**2/(flux*totScale)  + 4*3.14*totBkg*(psfMin**4)/(flux*totScale)**2 ))
    if(psf < psfMin-3*size_err):
        df_source['star_bool'][j] = 5
    boolArr.append(0)
    
    
    
df_source.to_pickle(outFile)       
fluxArr = np.array(fluxArr) 
sizeArr = np.array(sizeArr)
boolArr = np.array(boolArr)
plt.subplot()
#plt.yscale('log')

loc = np.array(df_source.index[df_source['star_bool'] != 1].tolist())
loc = np.where((boolArr != 1)  & (fluxArr >0))[0]
plt.figure(figsize=(13,10))
plt.plot(sizeArr[loc], 25-2.5*np.log10(fluxArr[loc]),'b.', markersize= 2)  
loc = np.array(df_source.index[df_source['star_bool'] == 1].tolist())
loc = np.where((boolArr == 1)  & (fluxArr >0))[0]
plt.plot(sizeArr[loc], 25-2.5*np.log10(fluxArr[loc]),'r.', markersize= 2)  
plt.xlabel('Size')
plt.ylabel('Magnitude')
plt.plot([psfMax, psfMin], [25-2.5*np.log10(1e5), 25-2.5*np.log10(1e5)], 'k-',markersize= 2)
plt.plot([psfMax, psfMin], [25-2.5*np.log10(fluxMin), 25-2.5*np.log10(fluxMin)], 'k-',markersize= 2)
plt.plot([psfMax, psfMax], [25-2.5*np.log10(1e5), 25-2.5*np.log10(fluxMin)], 'k-',markersize= 2)
plt.plot([psfMin, psfMin], [25-2.5*np.log10(1e5), 25-2.5*np.log10(fluxMin)], 'k-',markersize= 2)

flux_arange = np.arange(0.2,1000, 0.1)
size_err = np.sqrt( 0.5 *( (psfMin)**2/(flux_arange*totScale)  + 4*3.14*totBkg*(psfMin**4)/(flux_arange*totScale)**2 ))
plt.plot(psfMin-3*size_err, 25-2.5*np.log10(flux_arange), 'y-')
plt.plot([0,12], [25-2.5*np.log10(1e4), 25-2.5*np.log10(1e4)], 'y-')
plt.plot([12,12], [25-2.5*np.log10(0.2), 25-2.5*np.log10(1e4)], 'y-')
plt.ylim(29,11)
plt.savefig(plotLoc+'flux_vs_size.png')
plt.close()


# =============================================================================
# xList, yList = helper.convertToXY(raList, decList, imgFile)
# shutil.copyfile('/scratch/bell/dutta26/abell_2390/abell_ir_coadd_crop.fits', '/scratch/bell/dutta26/abell_2390/abell_ir_coadd_crop1.fits')
# f=fits.open('/scratch/bell/dutta26/abell_2390/abell_ir_coadd_crop1.fits', mode = 'update')
# data = f[0].data
# sizey, sizex = np.shape(data)
# for j in range(len(xList)):
#     x= int(round(xList[j]))
#     y= int(round(yList[j]))
#     if(x<20 or x>sizex-20 or y<20 or y>sizey-20):
#         continue
#     data[y-19:y+19, x-19] = 1000
#     data[y-19:y+19, x+19] = 1000
#     data[y-19, x-19:x+19] = 1000
#     data[y+19, x-19:x+19] = 1000
#     
# f.flush()
# =============================================================================
