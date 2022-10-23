#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 24 09:53:15 2021

@author: anirban
"""

import astropy.units as u
import helper
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
import sys

Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" # Select early Data Release 3
finalStarFile = str(sys.argv[1])
f=fits.open(str(sys.argv[2]))

img = np.array(f[0].data)
sizey, sizex = np.shape(img)

#Find the max, min and mid Ra Dec
w=wcs.WCS(f[0].header)
f.close()
[[startRa, startDec]]  = w.wcs_pix2world([[0,0]], 0)
[[endRa, endDec]] = w.wcs_pix2world([[sizex, sizey]], 0)
midRa = 0.5*(startRa +endRa )
midDec = 0.5* (startDec + endDec)

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
df.drop(df[df['phot_g_mean_mag'] < 16].index, inplace = True)

#df1 = df[['ra', 'dec', 'parallax_over_error']]
df2 = df[df['parallax_over_error'].notna()]

#Read the catalog 
catalog = str(sys.argv[3])
raArr =[] 
decArr = []
fluxArr = []
xPosArr =[]
yPosArr =[]
with open(catalog) as f:
    content = f.readlines()

#Store catalog elemnts in array
for j in range(len(content)):
    if((content[j].split())[0] == '#'):
       continue
    raArr.append(float(  (content[j].split())[5]) )
    decArr.append(float(  (content[j].split())[6]) )
    fluxArr.append(float(  (content[j].split())[1]) )
    xPosArr.append(float(  (content[j].split())[3]) )
    yPosArr.append(float(  (content[j].split())[4]) )

c_icrs = SkyCoord(ra=raArr*u.degree, dec=decArr*u.degree, frame='icrs')
temp = c_icrs.transform_to(FK5(equinox='J2000.0'))
raArr_2000 = np.array(temp.ra)
decArr_2000 = np.array(temp.dec)

count1 = count2 = count3 = count4 = 0    
finalMatchIndexArr =[]
for j in range(len(df2)):
    ra = float( df2['ra'][j:j+1])
    dec = float( df2['dec'][j:j+1])
    
    fluxMatch=0
    matchFlag = 0 
    matchIndexArr=[]
    distArr =[]
    #print ('$$$$$$$$$$$$$$$$$$')
    for k in range(len(raArr_2000)):
        if(np.abs(raArr_2000[k] - ra) < 0.003 and np.abs(decArr_2000[k] - dec) < 0.003):
            #print ('**********************')
            matchIndexArr.append(k)
            distArr.append(np.sqrt((raArr_2000[k] - ra)**2 + (decArr_2000[k] - dec)**2))
    
    #If no match is found       
    if(len(matchIndexArr) == 0):
        finalMatchIndexArr.append(None)
        count1 += 1
        continue
    #If only 1 match is found and is close 
    if(len(matchIndexArr) == 1 and distArr[0]<= 0.0002):
        finalMatchIndexArr.append(matchIndexArr[0])
        count2 += 1
        continue
    #If multiple matches are found 
    temp_dist = np.sort(distArr)
    if(temp_dist[0]< 0.0002 and temp_dist[1]> 0.002):
        leastDistIndex = np.where(distArr == temp_dist[0])
        finalMatchIndexArr.append(matchIndexArr[leastDistIndex[0][0]])
        count3 += 1
        continue
    else:
        finalMatchIndexArr.append(None)
        count4 += 1
        continue
    
store = np.zeros((10000, 11), dtype =np.float32)

#Now loop through the stars and find7 flux and psf etc 
for j in range(len(finalMatchIndexArr)):
    index = finalMatchIndexArr[j]
    if(index == None):
        continue
    y= int(yPosArr[index])
    x = int( xPosArr[index])
    flux_gaia =  float( df2['phot_g_mean_mag'][j:j+1])
    
    print (x,y)
    if((y - 25)< 0 or (x-25)<0 or (x+26)>= sizex-1 or (y+26)>= sizey-1):
        continue
    
    cut = img[y - 25: y+26, x-25:x+26]
    center = img[y - 5: y+6, x-5:x+6]
    flag = 0
    if(len(np.where(cut >300)[0]) > 0 ):
        flag = helper.vert_stripe(cut)  
    
    #Reject if any central pixels are 0    
    if(len(np.where(center < 0 )[0]) > 0 )
        flag = 1
        
    if(flag == 0):
        flux, mux, muy, e1, e2, back, size, sigxx, sigyy, sigxy = measure_pythonV.measure(cut)
        
        if(flux == None):
            
            store[j,:] = None, None, None,None, None, None, None, None, None, None, None
            continue
       
        store[j,:] = flux, mux, muy, e1, e2, back, size, sigxx, sigyy, sigxy, flux_gaia
    else:
        store[j,:] = None, None, None,None, None, None, None, None, None, None, None
        
store1 = store[store[:,0] != 0]
size = store1[:,6]
e1 = store1[:,3]
e2 = store1[:,4]
sys.exit()

mean_size, median_size, stddev_size = sigma_clipped_stats(size, cenfunc='median' )
mean_e1, median_e1, stddev_e1 = sigma_clipped_stats(e1, cenfunc='median' )
mean_e2, median_e2, stddev_e2 = sigma_clipped_stats(e2, cenfunc='median' )

maxSize, minSize = median_size + 3*stddev_size,  median_size - 3*stddev_size
maxe1, mine1 = median_e1 + 3*stddev_e1 , median_e1 - 3*stddev_e1
maxe2, mine2 = median_e2 + 3*stddev_e2 , median_e2 - 3*stddev_e2


#Now find the final ra dec of these images in the 3 sigma range 
finalRaList =[]
finalDecList =[]
f=open(finalStarFile, 'w+')
for j in range(len(finalMatchIndexArr)):
    size = store[j,6]
    e1 = store[j,3]
    e2 = store[j,4]
    if(size>minSize and size<maxSize and e1>mine1 and e1<maxe1 and e2>mine2 and e2<maxe2):
        finalRaList.append(float( df2['ra'][j:j+1]))
        finalDecList.append(float( df2['dec'][j:j+1]))
        f.write(str(float( df2['ra'][j:j+1])) + ', '+str(float( df2['dec'][j:j+1])) + '\n')
        
        
f.close()
        







