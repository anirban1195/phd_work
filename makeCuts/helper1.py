#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 20 08:17:44 2021

@author: dutta26
"""
import numpy as np
from astropy.io import fits


def vert_stripe(img):
    sizey, sizex = np.shape(img)
    #If cut too small
    if(sizex< 30 or sizey<30):
        return None
    flag = 0
    midx = int(round(sizex/2))
    midy = int(round(sizey/2))
    #Check 2 pix above centre 
    avg = np.mean(img[midy+2 , midx-1:midx+1 ])
    edge = (img[midy+2 , midx-6 ]+img[midy+2 , midx+6 ])/2
    #print (avg,edge, img[midy+2 , midx-1:midx+2 ] , img[midy+2 , midx-4 ], img[midy+2 , midx+4 ])
    if(avg < edge):
        flag = 1
    #Check 3 pix above center 
    avg = np.mean(img[midy+3 , midx-1:midx+1 ])
    edge = (img[midy+3 , midx-7 ] + img[midy+3 , midx+7 ])/2
    #print (avg,edge)
    
    if(avg < edge):
        flag = 1
    
    if(flag == 0):
        return 0
    else:
        return 1
    
    

    
#First read teh ra dec list of detected objects 
#Convert them to pixels of the single frame 
#Throw away any objects that are less than 20 pix near 
#Throw away any object that has lot of nans or zeros 
#The do ellipticity measurements 

def detectBad(cutout, size=3):
    cutout = np.array(cutout, dtype = np.float32)
    sizey, sizex = np.shape(cutout)
    flag = flagC = flagE = flagSize = flagTot = 0
# =============================================================================
#     cutout[np.isnan(cutout)] = 0
#     flag = flagC = flagE = flagSize = flagTot = 0
#     
#     #First check if any central 20x20 region is 
#     a= np.where(cutout[int(sizey/2)-5:int(sizey/2)+5 , int(sizex/2)-5:int(sizex/2)+5] <= 0)
#     if(len(a[0]) > 0):
#         flagC = 1
#         
#     
#         
#     #Next check how many total zero in cutout
#     b= np.where(cutout <= 0)
#     if(len(b[0]) > (0.2*sizex*sizey)):
#         flagTot = 1
#     
#     #Check the edges    
#     subCut = cutout[4:sizey-4, 4:sizex-4]
#     c= np.where(subCut <= 0)
#     if(len(c[0]) > 0):
#         flagE = 1
# =============================================================================
    
    #Check if zeros within 3 sigma of size 
    if(sizex> (6*size) and sizey>(6*size)):
        #print ('aa')
        centerx = int(round(sizex/2))
        centery = int(round(sizey/2))
        half = int(round(3*size))
        subCut = cutout[centery-half:centery+half, centerx-half:centerx+half]
        c= np.where(subCut <= 0)
        #print(len(c[0]))
        if(len(c[0]) > 0):
            flagC= 1
    else:
        #print ('bb')
        c= np.where(cutout <= 0)
        if(len(c[0]) > 0):
            flagC= 1
    
    #FlagC says there are 0s in central areas 
    #Flag tot sas there are too many zeros in the image to be useful
    #FlagE says there are zeros in places other than edges 
    if((flagE == 0 and flagC == 0 and flagTot == 0) ):
        flag = 0
    else:
        flag = 1
    #del data        
    return flag 






#This function checks each cutout and returns a list of good indices
#filename ,  xposition, yposition, and sizeList are inputs
#Size is total size of the cutout

def detectBadListV(filename,  xList, yList, sizeList):
    #print ('aaa')
    f=fits.open(filename)
    data = np.array(f[0].data)
    data[np.isnan(data)] = 0   #Make all nan values 0
    sizex, sizey = np.shape(data)
    f.close()
    indexList=[]
    for j in range(len(xList)):
        #print (j)
        flagC = 0 #Flag for central pizels 
        flagTot = 0 #Flag for total pixels
        flagE = 0 #Flag for edges
        sigma = int(round(1.5*sizeList[j]))
        
        if(xList[j]<0 or xList[j]> sizex):
            flagC = flagTot = 1
            continue
        if(yList[j]<0 or yList[j]> sizey):
            flagC = flagTot = 1
            continue
        
        if(sigma< 20):
            sigma = 20
        if(sigma>50):
            sigma=50
        cutout = data[int(round(yList[j]-sigma)) : int(round(yList[j]+sigma)) , int(round(xList[j]-sigma)) : int(round(xList[j]+sigma))]
        #cutout = data[100:200, 100:200]
        midPt = int(round(sigma))
        a = np.where(cutout[midPt-5:midPt+5, midPt-5:midPt+5] <= 0)
        #a = np.where(cutout[0:50, 0:50] <= 0)
        if(len(a[0]) > 0):
            flagC = 1
        
            
        #Next check how many total zero in cutout. If >100 then reject
        b= np.where(cutout <= 0)
        if(len(b[0]) > 0.2*sigma*sigma):
            flagTot = 1
            
        #Check the edges    
        subCut = cutout[4:sizey-4, 4:sizex-4]
        c= np.where(subCut <= 0)
        if(len(c[0]) > 0):
            flagE = 1
        
        #print (flagC, flagTot)
        if( flagC == 0 and flagTot == 0 and flagE == 0):
            indexList.append(j)
    del data        
    return indexList