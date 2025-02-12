#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 18:58:26 2021

@author: dutta26
"""


from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import os
import sys
from astropy import wcs 
from astropy.nddata.utils import Cutout2D
from astropy.stats import sigma_clipped_stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
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
    tempList = np.zeros((len(raList), 2), dtype = np.float64)
    tempList[:,0] = raList
    tempList[:,1] = decList
    temp = w.all_world2pix(tempList, 0)
    #temp = w.wcs_world2pix(tempList, 0)
    #temp = w.world_to_pixel(tempList)
    xList = temp[:,0]
    yList = temp[:,1]
    return xList, yList


#Takes list of ra dec objects and return a list of x and y objects
def convertToRaDec(xList, yList, filename):
    f=fits.open(filename)
    w = wcs.WCS(f[0].header)
    f.close()
    raList= []
    decList =[] 
    #First convert ra dec to x and y
    tempList = np.zeros((len(xList), 2), dtype = np.float64)
    tempList[:,0] = xList
    tempList[:,1] = yList
    temp = w.all_pix2world(tempList, 0)
    #temp = w.wcs_world2pix(tempList, 0)
    #temp = w.world_to_pixel(tempList)
    raList = temp[:,0]
    decList = temp[:,1]
    return raList, decList

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
    
    
def getLoc (raList, decList , filename):
    f=fits.open(filename)
    
    raStartArr =np.zeros((31 ,8 ,8), dtype = np.float32)
    raEndArr =np.zeros((31 ,8 ,8), dtype = np.float32)
    decStartArr=np.zeros((31 ,8 ,8), dtype = np.float32)
    decEndArr = np.zeros((31 ,8 ,8), dtype = np.float32)
    for j in range(1,31):
         w = wcs.WCS(f[j].header)
         for y in range(8):
             for x in range(8):
                 xStart = 5+(x*508)
                 xEnd = xStart + 470
                 yStart = 5+ (y*505)
                 yEnd = yStart + 484
                 [[raStartArr[j,y,x], decStartArr[j,y,x]]] = w.wcs_pix2world([[yStart, xStart]], 0)
                 [[raEndArr[j,y,x], decEndArr[j,y,x]]] = w.wcs_pix2world([[yEnd, xEnd]], 0)
                 
    
    segArr=[]
    chip_x =[]
    chip_y =[]
    
    for j in range(len(raList)):
        loc = np.where((raStartArr<raList[j]) & (raEndArr>raList[j]) & (decStartArr<decList[j]) & (decEndArr>decList[j]))
        #print (loc[0])
        if(len(loc[0] )== 0):
            segArr.append(-1)
            chip_x.append(-1)
            chip_y.append(-1)
        else:
            segArr.append(loc[0][0])
            chip_x.append(loc[1][0])
            chip_y.append(loc[2][0])
    return segArr, chip_x, chip_y    




#Takes in measured sigmas (so sigmas should be theoretically > psf but are not)
#returns sigmas that are >psf
def correct(sigmaxx_m , sigmayy_m, sigmaxy_m, sigmaxx_p, sigmayy_p, sigmaxy_p, exx_m, eyy_m, exy_m, exx_p, eyy_p, exy_p, niter = 30000):
    g1= np.random.normal(0, 1, niter)
    g2= np.random.normal(0, 1, niter)
    g3= np.random.normal(0, 1, niter)
    g4= np.random.normal(0, 1, niter)
    g5= np.random.normal(0, 1, niter)
    g6= np.random.normal(0, 1, niter)
    
    sigmaxx_c_arr = sigmaxx_m + exx_m * g1 - sigmaxx_p - exx_p*g2
    sigmayy_c_arr = sigmayy_m + eyy_m * g3 - sigmayy_p - eyy_p*g4
    sigmaxy_c_arr = sigmaxy_m + exy_m * g5 - sigmaxy_p - exy_p*g6
    size_arr = sigmaxx_c_arr+sigmayy_c_arr
    temp_arr = sigmaxx_c_arr + sigmayy_c_arr - 2*np.abs(sigmaxy_c_arr)
    
    #loc = np.where((sigmaxx_c_arr >0) & (sigmayy_c_arr >0) & (temp_arr>0))[0]
    loc = np.where((size_arr >0)  )[0]
    #print (len(loc))
    if(len(loc) == 0):
        return None,None,None,0  
    sigmaxx_c = np.median(sigmaxx_c_arr[loc])
    sigmayy_c = np.median(sigmayy_c_arr[loc])
    sigmaxy_c = np.median(sigmaxy_c_arr[loc])
    
    return sigmaxx_c, sigmayy_c, sigmaxy_c, len(loc)  




    
def findEffSigma(sigma_arr):
    sigma_arr = 2* np.array(sigma_arr)
    tot = 1/np.sum(1/sigma_arr)
    return tot*len(sigma_arr)/2


def findEffSigmaxy( sigmaxx_arr, sigmayy_arr, sigmaxy_arr):
    sigmaxx_arr = np.array(sigmaxx_arr)
    sigmayy_arr = np.array(sigmayy_arr)
    sigmaxy_arr = 2*np.array(sigmaxy_arr)
    rho = sigmaxy_arr/ (np.sqrt(2*sigmaxx_arr) * np.sqrt(2*sigmayy_arr) )
    #print (sigmaxy_arr, (np.sqrt(2*sigmaxx_arr) * np.sqrt(2*sigmayy_arr) ))
    rhoEff = np.mean(rho)
    sigxx_eff = findEffSigma(sigmaxx_arr)
    sigyy_eff = findEffSigma(sigmayy_arr)
    sigxy_eff = (rhoEff * np.sqrt(2*sigxx_eff) * np.sqrt(2*sigyy_eff))/2
    return sigxx_eff, sigyy_eff, sigxy_eff


def make2Dsubcuts( df_slice,  seqNo, writeLoc, sfLoc, band, j):
    xPos = df_slice[10]
    yPos = df_slice[11]
    
    size = int(round(np.sqrt(df_slice[38] +  df_slice[39])))
    for files in os.listdir(sfLoc):
        
        if(str(df_slice[44]) in files):
            
            file1 = sfLoc +files
            hdu = fits.open(file1)[0]
            wcs = WCS(hdu.header)
            
            
            if(size<6):
                size = 6
            size = int(8*size)
            cutout = Cutout2D(hdu.data,(xPos, yPos), (size, size), wcs=wcs)

           
            hdu.data = cutout.data
            hdu.header.update(cutout.wcs.to_header())

            # Write the cutout to a new FITS file
            
            cutout_filename = writeLoc + str(seqNo)+'_' +str(band)+'_'+str(j)+ '.fits'
            hdu.writeto(cutout_filename, overwrite=True)
            

    
def measure_new(img ,lut1, lut2, guessFlux = 100, guessmux = 0, guessmuy=0 , guessAlphax =3, guessAlphay =3, guessAlphaxy =0 , counter_target = 100, distort = 0, fixBkg=0, back_calc=0):
    #Shape the image properly
    img = np.array(img)
    sizey, sizex = np.shape(img)
    
    #Resize the cutout if too big
    if(sizex > 100):
        midX = int(round(sizex/2.0))
        img = img[: , midX-50: midX+50]
    if(sizey > 100):
        midY = int(round(sizey/2.0))
        img = img[midY-50: midY+50, :]
    if(sizex< 20 or sizey<20):
        return None, None, None,None, None, None, None, None, None, None
    
    #If guess sigmas it too small or large
    if(guessAlphax< 0.9):
        guessAlphax = 0.9
    if(guessAlphax> 10):
        guessAlphax = 10
        
    if(guessAlphay< 0.9):
        guessAlphay = 0.9
    if(guessAlphay>10):
        guessAlphay = 10
        
    if(abs(guessAlphaxy)> 100):
        guessAlphaxy = 0
    
    
    #Now define meshgrid for calculation
    sizey, sizex = np.shape(img)
    x = np.linspace(0, sizex-1, sizex)
    y = np.linspace(0, sizey-1, sizey)
    x= x -sizex/2.0 + 0.5 
    y= y -sizey/2.0 + 0.5 
    x, y = np.meshgrid(x, y)
    
    #Initialize variables
    delSigxx = 9999
    delSigyy = 9999

    prevSigxx = 9999
    prevSigyy = 9999
    
    alphax = guessAlphax
    alphay = guessAlphay
    alphaxy = guessAlphaxy
    sigxx_calc = sigyy_calc = sigxy_calc = 0
    mux_calc = guessmux
    muy_calc= guessmuy
    flux_calc = 0 
    e1 = e2 = 0.0
    
    med = np.median(img)
    if(fixBkg==0):
        back_calc = np.median(img)
    
    total = np.sum(img)
    counter = 0
    
    
    
    #Loop until convergence
    
    while( (abs(delSigxx)>0.001 or abs(delSigyy)> 0.001) and counter<counter_target):
        #print (alphax, alphay,alphaxy)
        #Correct alphxy to avoid nans
        while( (alphaxy/(alphax*alphay))**2  >= 1):
            #print (alphax, alphay, alphaxy)
            if(alphaxy > 0 ):
                alphaxy = alphaxy - 0.1
            if (alphaxy < 0):
                alphaxy = alphaxy + 0.1
            if (abs(alphaxy) > 5000):
                counter = counter_target+999
                break
        #If size too large exit
        if(abs(alphaxy)> sizex/2  or abs(alphax) > sizex/2 or abs(alphay)> sizey/2):
            counter = counter_target+999
            break
        
        arbConst = 2*(1- (alphaxy/(alphax*alphay))**2 )  
        A =1/(2*np.pi*alphax*alphay*np.sqrt(1- (alphaxy/(alphax*alphay))**2 ))
        k=(A * np.exp(-((x-mux_calc)**2/(arbConst*alphax**2)+ (y-muy_calc)**2/(arbConst*alphay**2) - 2*alphaxy*(y-muy_calc)*(x-mux_calc)/(arbConst* alphax**2 * alphay**2 ) )))
        
        #Apply distortion correction
        if(distort == 1):
            sqrtImg = np.sqrt(img)
            q= np.abs(img-back_calc)/sqrtImg
            q2 = q*q
            q3 = q2*q
            temp1 = img - back_calc + 1.41421* sqrtImg*(0.477/ np.exp(0.928*q +0.247*q2 + 0.04*q3))
            temp2 = 1.41421* sqrtImg*(0.477/ np.exp(0.551*q - 0.06*q2 + 0.003*q3))
            img1 = img - back_calc
            img1[img1>=0] = temp1[img1>=0]
            img1[img1<0] = temp2[img1<0]
        else:
            img1 = img - back_calc
        
        t1= np.sum(x*y* (img1) * k)
        t2 = np.sum(k*(img1))
        t3 = np.sum(x*(img1)*k)
        t4 = np.sum(y*(img1)*k)
        t5 = np.sum(x * x * (img1) * k)
        t6 = np.sum(y * y * (img1) * k)
        t7 = np.sum(k**2)
        
        mux_calc = t3/t2
        muy_calc = t4/t2
        flux_calc = t2/t7
        
        
        sigxy_calc = (t1/t2) - (t3*t4)/(t2*t2)
        sigxx_calc = (t5/t2) - (t3*t3) / (t2*t2)
        sigyy_calc = (t6/t2) - (t4*t4) / (t2*t2)
        
        
        #If negative sigmas, break
        if(sigxx_calc<0 or sigyy_calc<0):
            counter = counter_target+999
            break
        
        e1 = (sigxx_calc - sigyy_calc)/(sigxx_calc + sigyy_calc)
        e2 = 2*sigxy_calc/(sigxx_calc + sigyy_calc)
        #ellip = np.sqrt(e1*e1 + e2*e2)
        
        if(med != 0 and fixBkg==0):
            back_calc = (total - flux_calc)/ (sizex*sizey)
            
        delSigxx = prevSigxx - sigxx_calc
        delSigyy = prevSigyy - sigyy_calc
        
        prevSigxx = sigxx_calc
        prevSigyy = sigyy_calc
        
        alphax = sigxx_calc*2
        alphay = sigyy_calc*2
        if(alphax <0.9 ):
            alphax = 0.9
        if(alphay <0.9 ):
            alphay = 0.9
        alphax = np.sqrt(alphax)
        alphay = np.sqrt(alphay)
        alphaxy = 2.0*sigxy_calc
       
        
        counter += 1
        
    #print (alphax, alphay, alphaxy)
    #If failed convergence, return None         
    if(counter == (counter_target+999)):
        return None, None, None,None, None, None, None, None, None, None
    else:
        if((guessAlphax**2 *guessAlphay**2 - guessAlphaxy**2)<=0 or guessFlux<=0):
            return None, None, None,None, None, None, None, None, None, None
        
        if((med< 0 or back_calc<0) and distort == 1):
            return None, None, None,None, None, None, None, None, None, None
        #Appy correction if distorted
        sqrtba = np.sqrt(med)*3.14*np.sqrt(guessAlphax**2 *guessAlphay**2 - guessAlphaxy**2)
        fbysqrtba = np.log10(guessFlux/sqrtba)
        if(distort == 1 and counter_target == 1):
            #print ('aa')
            #lut = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
            lut = lut1
            if(fbysqrtba< -3.275):
                index = 0
                sigma_factor = lut[index,2]
                flux_bias= lut[index,1]  * sqrtba
            elif(fbysqrtba> 2.1):
                index = 217
                sigma_factor = lut[index,2]
                flux_bias= lut[index,1]  * sqrtba
            else:
                index = int((fbysqrtba + 3.3)/0.025)
                sigma_factor = np.interp(fbysqrtba,lut[index-1:index+2,0], lut[index-1:index+2,2] )
                flux_bias= np.interp(fbysqrtba,lut[index-1:index+2,0], lut[index-1:index+2,1] ) * sqrtba
            flux_calc -= flux_bias
            sigxx_calc -= sigma_factor*(guessAlphax)**2
            sigyy_calc -= sigma_factor*(guessAlphay)**2
            sigxy_calc -= sigma_factor*(guessAlphaxy)    
                
        elif(distort == 1 and counter_target > 1):
            #print ('bb')
            #lut = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
            lut = lut2
            if(fbysqrtba< 1.225):
                index = 0
                sigma_factor = lut[index,2]
                flux_bias= lut[index,1]  * sqrtba
            elif(fbysqrtba> 2.825):
                index = 34
                sigma_factor = lut[index,2]
                flux_bias= lut[index,1]  * sqrtba
            else:
                index = int((fbysqrtba -1.175)/0.05)
                sigma_factor = np.interp(fbysqrtba,lut[index-1:index+2,0], lut[index-1:index+2,2] )
                flux_bias= np.interp(fbysqrtba,lut[index-1:index+2,0], lut[index-1:index+2,1] ) * sqrtba
            flux_calc -= flux_bias
            sigxx_calc -= sigma_factor*(guessAlphax)**2
            sigyy_calc -= sigma_factor*(guessAlphay)**2
            sigxy_calc -= sigma_factor*(guessAlphaxy)
            
        
        return flux_calc, mux_calc, muy_calc, e1, e2, back_calc, np.sqrt(sigxx_calc + sigyy_calc), sigxx_calc, sigyy_calc, sigxy_calc


def findBoundRadius(store):
    boundArr =[]
    for j in range(len(store)):
        thresh =0.11
        flux = store[j,3]
        fact = 1
        if(flux>1000):
            fact = 2
        elif(flux<200):
            fact = 1
        else:
            fact = (flux-200)/800 +1
        maxVal = max(store[j,7] , store[j,8])
        size = np.sqrt(maxVal*2)
        if(np.isnan(size) or np.isnan(flux) or flux<=0 or size<=0):
            boundArr.append(0)
            continue
        r= fact*np.sqrt( -2* size**2 * np.log(thresh*2*np.pi*size**2/ flux))
        if(r<13):
            r= 13
        boundArr.append(r)
    
    flagArr=[]
    minArr=[]
    angleArr=[]
    for j in range(len(store)):
        print (j)
        flux = store[j,3]
        size = np.sqrt(store[j,7] + store[j,8])
        if(np.isnan(size) or np.isnan(flux) or flux<=0 or size<=0):
            flagArr.append(0)
            minArr.append(0)
            angleArr.append(0)
            continue
        x,y = store[j,10], store[j,11]
        dist = 999999999
        flag = 0
        maxLt = j+1000
        minLt = j-1000
        if(minLt<0):
            minLt = 0
        if(maxLt>len(store)):
            maxLt = len(store)
        minDist =99999999
        angle = 0
        for k in range(minLt, maxLt):
            if(j==k):
                continue
            x1,y1 = store[k,10], store[k,11]
            dist = np.sqrt((x-x1)**2 + (y-y1)**2)
            if(dist>200):
                continue
            crit_dist = boundArr[j] + boundArr[k]
            if(dist<crit_dist):
                flag = 1
            if(dist<minDist):
                minDist =dist
                angle =np.arctan2(y1-y, x1-x)
        flagArr.append(flag)
        minArr.append(minDist)
        angleArr.append(angle)
    return boundArr, flagArr, minArr, angleArr



def bkgFlagCalculate(img, store,  arr_min=6000, arr_max=15000):
    flagArr=[]
    angleArr= np.arange(0, 3.14+0.01, 3.14/4)
    gmean, gmedian, gstd = sigma_clipped_stats(img[arr_min:arr_max, arr_min:arr_max])
    for j in range(len(store)):
        #print (j, 'aa')
        flag =0
        thresh = gstd
        flux = store[j,3]
        fact = 1
# =============================================================================
#         if(flux>1000):
#             fact = 2
#         elif(flux<200):
#             fact = 1
#         else:
#             fact = (flux-200)/800 +1
# =============================================================================
        maxVal = max(store[j,7] , store[j,8])
        size = np.sqrt(maxVal*2)
        if(np.isnan(size) or np.isnan(flux) or flux<=0 or size<=0):
            flagArr.append(99)
            continue
        test = -2* size**2 * np.log(thresh*2*np.pi*size**2/ flux)
        if(test< 0):
            #print (thresh*2*np.pi*size**2/ flux , flux, size)
            r= 13
        else:
            r= fact*np.sqrt( -2* size**2 * np.log(thresh*2*np.pi*size**2/ flux))
        if(r<13):
            r= 13
        if(np.isnan(r) or r==None):
            flag = 99
            flagArr.append(flag)
            continue
        x,y = int(round(store[j,10])), int(round(store[j,11]))
        r1 = int(round((r/1.414)))
        
        med1 = np.median(img[y+r1:y+r1+4 , x+r1:x+r1+4] )
        med2 = np.median(img[y+r1:y+r1+4 , x-r1-4:x-r1] )
        med3 = np.median(img[y-r1-4:y-r1, x-r1-4:x-r1] )
        med4 = np.median(img[y-r1-4:y-r1 , x+r1:x+r1+4] )
        r = int(round(r))
        med5 = np.median(img[y-2:y+2 , x+r:x+r+4] )
        med6 = np.median(img[y-2:y+2 , x-r-4:x-r] )
        med7 = np.median(img[y+r:y+r+4 , x-2:x+2] )
        med8 = np.median(img[y-r-4:y-r , x-2:x+2] )
        
        lt_low = gmedian - 3*gstd
        lt_high = gmedian + 3*gstd
        tempArr= np.array([med1, med2, med3, med4, med5, med6, med7, med8])
        for val in tempArr:
            if(np.isnan(val) or val == None):
                flag = 1
        l1 = l2= 0
        #l1 = len(np.where( (tempArr> (gmedian + 3*gstd) ) | (tempArr< (gmedian - 3*gstd) ) ) [0])
        #l2 = len(np.where( (tempArr> (gmedian + 2*gstd) ) | (tempArr< (gmedian - 2*gstd) ) ) [0])
        #l3 = len(np.where( (tempArr>gmedian + 1*gstd) | (tempArr<gmedian - 1*gstd)))
        if(l1>=1  or l2>=2):
            flag = 1
        
        #Check for self consistency
        locMean, locMed, locStd = sigma_clipped_stats(tempArr, sigma=3)
# =============================================================================
#         if(locMed >  (gmedian + 3*gstd) or locMed <  (gmedian + 3*gstd)):
#             flag = 1
# =============================================================================
        
        locStd = gstd/2
        if(locMed > (gmedian+2*gstd)):
            locStd = locMed/3
        lc1 = len(np.where( (tempArr> (locMed + 3*locStd) ) | (tempArr< (locMed - 3*locStd) ) ) [0])
        if(lc1 >= 1):
            flag = 1
        if(j in [8336, 9046, 8148]):
            print (gmedian, gstd, med1, med2, med3, med4, med5, med6, med7,med8, flag,j, locStd, locMed,r)
        flagArr.append(flag)
        
    return flagArr
        
            
        
def getPSFErr(xxArr, yyArr, xyArr, fluxArr, wtArr, B):
    wtArr = wtArr/np.sum(wtArr) #Normalize the weights
    errArr=[]
    
    size = np.sqrt(xxArr + yyArr)
    #print (size)
    errArr=np.sqrt(size**4/fluxArr + 4*size**6*np.pi * B/(fluxArr**2)  ) 
     
    return np.sqrt(np.sum(wtArr**2*errArr**2))



def gaussian(x, mu, sig, A):
    return (A)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))




def fitGaussian(xxArr, yyArr, xyArr, fluxArr, bkg, plotLoc, band):
    
    
    mean_psf_size = np.median(np.sqrt(xxArr+yyArr))
    err_xx = err_yy = np.sqrt( (mean_psf_size**4/np.median(fluxArr) + (4* np.pi* mean_psf_size**6 * bkg)/(np.median(fluxArr)**2) ) )
    err_xy = err_xx * 0.707

    k_guess_xx = np.std(xxArr)
    k_guess_yy = np.std(yyArr)
    k_guess_xy = np.std(xyArr)
    
    mean_xx = np.median( xxArr)
    mean_yy = np.median( yyArr)
    mean_xy = np.median( xyArr)
    
    
    #Do for xx
    data = xxArr
    counts,bin_edges = np.histogram(data, bins=40)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.0
    err = []
    for j in range(1, len(bin_edges)):
        bin_start = bin_edges[j-1]
        bin_end = bin_edges[j]
        loc1 = np.where((data>=bin_start)& (data<bin_end))[0]
        if(len(loc1) <=5):
            err.append(0)
        else:
            err.append(np.nanstd(data[loc1]))
            
    plt.errorbar(bin_centres, counts, xerr=None, yerr=np.sqrt(counts), fmt='o', markersize = 5)
    xVal =[]
    yVal =[]
    tot = 0
    start = mean_xx - 5*np.sqrt(k_guess_xx**2 + err_xx**2)
    cutoff = mean_xx + 5*np.sqrt(k_guess_xx**2 + err_xx**2)
    #center = 13.6
    #print (mean_xx, k_guess_xx)
    for j in range(len(bin_centres)):
        
        if(bin_centres[j]> cutoff or bin_centres[j]< start):
            continue
        xVal.append(bin_centres[j])
        yVal.append(counts[j])
    
    
    parameters, covariance = curve_fit(gaussian, xVal, yVal, [mean_xx, np.sqrt(k_guess_xx**2 + err_xx**2) , 80])    
    res_err_xx = np.sqrt(parameters[1]**2 - err_xx**2)
    #print (parameters[1], err_xx)
    x_values = np.linspace(mean_xx - 10*np.sqrt(k_guess_xx**2 + err_xx**2),mean_xx + 10*np.sqrt(k_guess_xx**2 + err_xx**2), 25000)
    plt.plot(x_values, gaussian(x_values,parameters[0], parameters[1], parameters[2]), 'r-',linewidth = 2)
    plt.savefig(plotLoc+'xx_fit_'+band+'.png')
    plt.xlabel('Sigma xx')
    plt.close()
    
    
    
    #Do for yy
    data = yyArr
    counts,bin_edges = np.histogram(data, bins=40)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.0
    err = []
    for j in range(1, len(bin_edges)):
        bin_start = bin_edges[j-1]
        bin_end = bin_edges[j]
        loc1 = np.where((data>=bin_start)& (data<bin_end))[0]
        if(len(loc1) <=5):
            err.append(0)
        else:
            err.append(np.nanstd(data[loc1]))
    plt.errorbar(bin_centres, counts, xerr=None, yerr=np.sqrt(counts), fmt='o', markersize = 5)       
    xVal =[]
    yVal =[]
    tot = 0
    start = mean_yy - 5*np.sqrt(k_guess_yy**2 + err_yy**2)
    cutoff = mean_yy + 5*np.sqrt(k_guess_yy**2 + err_yy**2)
    #center = 13.6
    for j in range(len(bin_centres)):
        if(bin_centres[j]> cutoff or bin_centres[j]< start):
            continue
        xVal.append(bin_centres[j])
        yVal.append(counts[j])
    
    
    parameters, covariance = curve_fit(gaussian, xVal, yVal, [mean_yy, np.sqrt(k_guess_yy**2 + err_yy**2) , 80])    
    res_err_yy = np.sqrt(parameters[1]**2 - err_yy**2)
    x_values = np.linspace(mean_yy - 10*np.sqrt(k_guess_yy**2 + err_yy**2),mean_yy + 10*np.sqrt(k_guess_yy**2 + err_yy**2), 25000)
    plt.plot(x_values, gaussian(x_values,parameters[0], parameters[1], parameters[2]), 'r-',linewidth = 2)
    plt.savefig(plotLoc+'yy_fit_'+band+'.png')
    plt.xlabel('Sigma yy')
    plt.close()
    
    #Do for xy
    data = xyArr
    counts,bin_edges = np.histogram(data, bins=40)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.0
    err = []
    for j in range(1, len(bin_edges)):
        bin_start = bin_edges[j-1]
        bin_end = bin_edges[j]
        loc1 = np.where((data>=bin_start)& (data<bin_end))[0]
        if(len(loc1) <=5):
            err.append(0)
        else:
            err.append(np.nanstd(data[loc1]))
    plt.errorbar(bin_centres, counts, xerr=None, yerr=np.sqrt(counts), fmt='o', markersize = 5)

    xVal =[]
    yVal =[]
    tot = 0
    start = mean_xy - 5*np.sqrt(k_guess_xy**2 + err_xy**2)
    cutoff = mean_xy + 5*np.sqrt(k_guess_xy**2 + err_xy**2)
    #center = 13.6
    for j in range(len(bin_centres)):
        if(bin_centres[j]> cutoff or bin_centres[j]< start):
            continue
        xVal.append(bin_centres[j])
        yVal.append(counts[j])
    
    
    parameters, covariance = curve_fit(gaussian, xVal, yVal, [mean_xy, np.sqrt(k_guess_xy**2 + err_xy**2) , 80])    
    res_err_xy = np.sqrt(parameters[1]**2 - err_xy**2)
    x_values = np.linspace(mean_xy - 10*np.sqrt(k_guess_xy**2 + err_xy**2),mean_xy + 10*np.sqrt(k_guess_xy**2 + err_xy**2), 25000)
    plt.plot(x_values, gaussian(x_values,parameters[0], parameters[1], parameters[2]), 'r-',linewidth = 2)
    plt.savefig(plotLoc+'xy_fit_'+band+'.png')
    plt.xlabel('Sigma xy')
    plt.close()
    
    return res_err_xx, res_err_yy, res_err_xy
    
    
    
    
            
        
    