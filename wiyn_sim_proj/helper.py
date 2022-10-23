#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 18:58:26 2021

@author: dutta26
"""


from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import sys
#import measure_pythonV


#This function checks each cutout and returns a list of good indices
#filename ,  xposition, yposition, and sizeList are inputs
#Size is total size of the cutout

def detectBad(filename,  xList, yList, sizeList):
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
        a = np.where(cutout[midPt-5:midPt+5, midPt-5:midPt+5] == 0)
        #a = np.where(cutout[0:50, 0:50] <= 0)
        if(len(a[0]) > 10):
            flagC = 1
        
            
        #Next check how many total zero in cutout. If >100 then reject
        b= np.where(cutout == 0)
        if(len(b[0]) > 0.5*sigma*sigma):
            flagTot = 1
        
        #print (flagC, flagTot)
        if( flagC == 0 and flagTot == 0):
            indexList.append(j)
    del data        
    return indexList


#Takes list of ra dec objects and return a list of x and y objects
def convertToXY(raList, decList, filename):
    f=fits.open(filename)
    w = wcs.WCS(f[0].header)
    f.close()
    xList= []
    yList =[] 
    #First convert ra dec to x and y
    tempList = np.zeros((len(raList), 2), dtype = np.float32)
    tempList[:,0] = raList
    tempList[:,1] = decList
    temp = w.wcs_world2pix(tempList, 0)
    xList = temp[:,0]
    yList = temp[:,1]
    return xList, yList

#Return new x and y list from only for indices present in indexList
def createNewList(xList, yList, sizeList, indexList):
    newxList=[]
    newyList =[]
    newSizeList=[]
    for j in range(len(indexList)):
        newxList.append(xList[indexList[j]])
        newyList.append(yList[indexList[j]])
        newSizeList.append(sizeList[indexList[j]])
    return newxList, newyList, newSizeList


def runMeasure(filename , xList, yList, indexList, sizeList):
    f=fits.open(filename)
    data = np.array(f[0].data)
    data[np.isnan(data)] = 0   #Make all nan values 0
    sizex, sizey = np.shape(data)
    f.close()
    fluxList =np.zeros(len(xList), dtype=np.float32)
    psfList = np.zeros(len(xList), dtype=np.float32)
    bkgList = np.zeros(len(xList), dtype=np.float32)
    e1List = np.zeros(len(xList), dtype=np.float32)
    e2List = np.zeros(len(xList), dtype=np.float32)
    muXList = np.zeros(len(xList), dtype=np.float32)
    muYList = np.zeros(len(xList), dtype=np.float32)
    
    for j in range(len(indexList)):
        #print (j)
        size = sizeList[indexList[j]]
        x = xList[indexList[j]]
        y = yList[indexList[j]]
        sigma = int(round(1.5*size))
        
        if(sigma< 20):
            sigma = 20
        cutout = data[int(round(y-sigma)) : int(round(y+sigma)) , int(round(x-sigma)) : int(round(x+sigma))]
        
        flux, mux, muy, e1, e2, bkg, psf = measure_pythonV.measure(cutout) 
        fluxList[indexList[j]] = flux
        psfList[indexList[j]] = psf
        bkgList[indexList[j]] = bkg
        e1List[indexList[j]]= e1
        e2List[indexList[j]] = e2
        muXList[indexList[j]] = mux
        muYList[indexList[j]] = muy
        del cutout
        
       
    return fluxList, psfList, bkgList, e1List, e2List, muXList, muYList







#First read teh ra dec list of detected objects 
#Convert them to pixels of the single frame 
#Throw away any objects that are less than 20 pix near 
#Throw away any object that has lot of nans or zeros 
#The do ellipticity measurements 

def detectBad1(filename, xList, yList):
    f=fits.open(filename)
    data = np.array(f[0].data)
    print (np.shape(data))
    f.close()
    indexList=[]
    for j in range(len(xList)):
        
        flagC = 0 #Flag for entral pizels 
        flagE = 0 #Flag edges
        flagTot = 0 
        cutout = data[int(round(yList[j]-25)) : int(round(yList[j]+25)) , int(round(xList[j]-25)) : int(round(xList[j]+25))]
        
        
        #First check if any central 20x20 region is 
        a= np.where(cutout[15:35, 15:35] <= 0)
        if(len(a[0]) > 0):
            flagC = 1
            
        #Next check how many total zero in cutout
        b= np.where(cutout <= 0)
        if(len(b[0]) > 100):
            flagTot = 1
        
        #Check the edges    
        subCut = cutout[3:47, 3:47]
        c= np.where(subCut <= 0)
        if(len(c[0]) > 0):
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
def detNearby1(raList, decList, thresh , filename):
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
    
    
def correct(sigmaxx_m , sigmayy_m, sigmaxy_m, sigmaxx_p, sigmayy_p, sigmaxy_p, exx_m, eyy_m, exy_m, exx_p, eyy_p, exy_p):
    g1= np.random.normal(0, 1, 10000)
    g2= np.random.normal(0, 1, 10000)
    g3= np.random.normal(0, 1, 10000)
    g4= np.random.normal(0, 1, 10000)
    g5= np.random.normal(0, 1, 10000)
    g6= np.random.normal(0, 1, 10000)
    
    sigmaxx_c_arr = sigmaxx_m + exx_m * g1 - sigmaxx_p - exx_p*g2
    sigmayy_c_arr = sigmayy_m + eyy_m * g3 - sigmayy_p - eyy_p*g4
    sigmaxy_c_arr = sigmaxy_m + exy_m * g5 - sigmaxy_p - exy_p*g6
    
    temp_arr = sigmaxx_c_arr + sigmayy_c_arr - 2*np.abs(sigmaxy_c_arr)
    
    loc = np.where((sigmaxx_c_arr >0) & (sigmayy_c_arr >0) & (temp_arr>0))[0]
    #print (len(loc))
    if(len(loc) == 0):
        return 0,0,0,0  
    sigmaxx_c = np.mean(sigmaxx_c_arr[loc])
    sigmayy_c = np.mean(sigmayy_c_arr[loc])
    sigmaxy_c = np.mean(sigmaxy_c_arr[loc])
    
    return sigmaxx_c, sigmayy_c, sigmaxy_c, len(loc)  
