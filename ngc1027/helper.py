#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 10:53:36 2020

@author: anirban
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 10:53:36 2020

@author: anirban
"""

from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import sys


#First read teh ra dec list of detected objects 
#Convert them to pixels of the single frame 
#Throw away any objects that are less than 20 pix near 
#Throw away any object that has lot of nans or zeros 
#The do ellipticity measurements 

def detectBad(filename, xList, yList):
    f=fits.open(filename)
    data = np.array(f[0].data)
    print (np.shape(data))
    f.close()
    indexList=[]
    for j in range(len(xList)):
        
        flagC = 0 #Flag for entral pizels 
        flagE = 0 #Flag edges
        flagTot = 0 
        cutout = data[int(round(yList[j]-25)) : int(round(yList[j]+26)) , int(round(xList[j]-25)) : int(round(xList[j]+26))]
        
        
        #First check if any central 10x10 region is 
        a= np.where(cutout[20:31, 20:31] <= 0)
        if(len(a[0]) > 0):
            flagC = 1
            
        #Next check how many total zero in cutout. If >100 then reject
        b= np.where(cutout <= 0)
        if(len(b[0]) > 100):
            flagTot = 1
        
        #Check the edges    
        subCut = cutout[5:45, 5:45]
        c= np.where(subCut <= 0)
        if(len(c[0]) > 10):
            flagE = 1
        
        if(j==5):
            print (cutout, len(b[0]), flagTot)
        #FlagC says there are 0s in central areas 
        #Flag tot sas there are too many zeros in the image to be useful
        #FlagE says there are zeros in places other than edges 
        if((flagE == 0 and flagC == 0 and flagTot == 0) ):
            indexList.append(j)
    del data        
    return indexList

#Gets a list of ra dec and threshold 
#Return a list of x and y coordinates that are good(thresh)
def detNearby(raList, decList, thresh , filename):
    f=fits.open(filename)
    w = wcs.WCS(f[0].header)
    f.close()
    xList= []
    yList =[] 
    print (len(raList))
    #First convert ra dec to x and y
    tempList = np.zeros((len(raList), 2), dtype = np.float32)
    print (np.shape(tempList))
    tempList[:,0] = raList
    tempList[:,1] = decList
    temp = w.wcs_world2pix(tempList, 0)
    xList = temp[:,0]
    yList = temp[:,1]
    print (raList[1:5],xList[1:5])
    #Now check x and y lists for nearby objects
    #If nearby then throw both objects away 
    indexList=[]
    for j in range(len(xList)):
        x=xList[j]
        y=yList[j]
        flag =0
        #Find the range of seach 
        #Saerch neighbors only to save time 
        minVal = j-2000
        maxVal = j+2000
        if(minVal < 0):
            minVal = 0
        if(maxVal > len(xList)):
            maxVal = len(xList)
        
        for k in range(minVal, maxVal):
            if(k == j):
                continue
            x1 = xList[k]
            y1 = yList[k]
            if(((x-x1)**2 + (y-y1)**2) < thresh**2):
                flag = 1
        #If flag=0 then it is a good object
        if(flag ==0 ):
            indexList.append(j)
    return indexList, xList, yList
    
    
    
def createNewList(xList, yList, indexList):
    newxList=[]
    newyList =[]
    for j in range(len(indexList)):
        newxList.append(xList[indexList[j]])
        newyList.append(yList[indexList[j]])
    return newxList, newyList
        
    
def detectBad_combined(filename, xList, yList):
    f=fits.open(filename)
    data = np.array(f[0].data)
    print (np.shape(data))
    f.close()
    indexList=[]
    for j in range(len(xList)):
        
        flagC = 0 #Flag for entral pizels 
        flagE = 0 #Flag edges
        flagTot = 0 
        cutout = data[int(round(yList[j]-25)) : int(round(yList[j]+26)) , int(round(xList[j]-25)) : int(round(xList[j]+26))]
        
        
        #First check if any central 10x10 region is 
        a= np.where(cutout[20:31, 20:31] <= -1)
        if(len(a[0]) > 0):
            flagC = 1
            
        #Next check how many total zero in cutout. If >100 then reject
        b= np.where(cutout <= -1)
        if(len(b[0]) > 100):
            flagTot = 1
        
        #Check the edges    
        subCut = cutout[5:45, 5:45]
        c= np.where(subCut <= -1)
        if(len(c[0]) > 10):
            flagE = 1
        
        if(j==5):
            print (cutout, len(b[0]), flagTot)
        #FlagC says there are 0s in central areas 
        #Flag tot sas there are too many zeros in the image to be useful
        #FlagE says there are zeros in places other than edges 
        if((flagE == 0 and flagC == 0 and flagTot == 0) ):
            indexList.append(j)
    del data        
    return indexList