#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 16:02:16 2023

@author: dutta26
"""

import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
def measure(img, guessAlphax =3, guessAlphay =3, guessAlphaxy =0 , counter_target = 100, bkg1= 0):
    #Shape the image properly
    img = np.array(img)
    sizey, sizex = np.shape(img)
    #print (sizex,sizey)
    
    if(sizex > 100):
        midX = int(round(sizex/2.0))
        img = img[: , midX-50: midX+50]
    if(sizey > 100):
        midY = int(round(sizey/2.0))
        img = img[midY-50: midY+50, :]
    if(sizex< 30 or sizey<30):
        return None, None, None,None, None, None, None, None, None, None
    #hdu = fits.PrimaryHDU(img)
    #hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/cpp_temp.fits', overwrite=True)
    #Now commence actual calculation
    sizey, sizex = np.shape(img)
    x = np.linspace(0, sizex-1, sizex)
    y = np.linspace(0, sizey-1, sizey)
    x= x -sizex/2.0 + 0.5 
    y= y -sizey/2.0 + 0.5 
    #mux=muy= np.random.normal(29.5, 0.15)
    x, y = np.meshgrid(x, y)
    
    delSigxx = 999
    delSigyy = 999
    delSigxy = 999
    prevSigxx = 9999
    prevSigyy = 9999
    if(guessAlphax< 0.9):
        guessAlphax = 0.9
    if(guessAlphay< 0.9):
        guessAlphay = 0.9
    if(abs(guessAlphaxy)> 100):
        guessAlphaxy = 0
    alphax = guessAlphax
    alphay = guessAlphay
    alphaxy = guessAlphaxy
    sigxx_calc =0
    sigyy_calc=0
    sigxy_calc = 0
    mux_calc = 0
    muy_calc = 0
    flux_calc = 0 
    e1 = 0.0
    e2 = 0.0
    med = np.median(img)
    if(bkg1 == 0):
        back_calc = med 
    else:
        back_calc = 0
    total = np.sum(img)
    counter = 0
    
    
    
    while(abs(delSigxx)>0.001 and abs(delSigyy)> 0.001 and counter<counter_target):
        #print (counter)
        
        while( (alphaxy/(alphax*alphay))**2  >= 1):
            #print (alphax, alphay, alphaxy)
            if(alphaxy > 0 ):
                alphaxy = alphaxy - 0.101235
            if (alphaxy < 0):
                alphaxy = alphaxy + 0.101235
            if (abs(alphaxy) > 5000):
                counter = 100
                break
        if(abs(alphaxy)> 50  or abs(alphax) > sizex/2 or abs(alphay)> sizey/2):
            counter = counter_target
            break
        arbConst = 2*(1- (alphaxy/(alphax*alphay))**2 )  
        A =1/(2*np.pi*alphax*alphay*np.sqrt(1- (alphaxy/(alphax*alphay))**2 ))
        k=(A * np.exp(-((x-mux_calc)**2/(arbConst*alphax**2)+ (y-muy_calc)**2/(arbConst*alphay**2) - 2*alphaxy*(y-muy_calc)*(x-mux_calc)/(arbConst* alphax**2 * alphay**2 ) )))
        
        #img1 = img-back_calc
        #img1[img1<0] = 0
        
        q= np.abs(img-back_calc)/np.sqrt(img)
        temp1 = img - back_calc + np.sqrt(2*img)*(0.477/ np.exp(0.928*q +0.247*q*q + 0.04*q*q*q))
        temp2 = np.sqrt(2*img)*(0.477/ np.exp(0.551*q - 0.06*q*q + 0.003*q*q*q))
        img1 = img - back_calc
        
        img1[img1>=0] = temp1[img1>=0]
        img1[img1<0] = temp2[img1<0]
        
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
        
        #print (t1, t2, t3, t4, t5, t6)
        #print (sigxy_calc, sigxx_calc, sigyy_calc)
        
        if(sigxx_calc<0 or sigyy_calc<0):
            counter = counter_target
            break
        
        e1 = (sigxx_calc - sigyy_calc)/(sigxx_calc + sigyy_calc)
        e2 = 2*sigxy_calc/(sigxx_calc + sigyy_calc)
        #ellip = np.sqrt(e1*e1 + e2*e2)
        
        if(med != 0 ):
            back_calc = (total - flux_calc)/ (sizex*sizey)
            #total = flux_calc + (sizex*sizey* back_calc)
        #back_calc = bkg
        
        delSigxx = prevSigxx - sigxx_calc
        delSigyy = prevSigyy - sigyy_calc
        
        prevSigxx = sigxx_calc
        prevSigyy = sigyy_calc
        
        alphax = sigxx_calc*2
        alphay = sigyy_calc*2
        if(alphax <0.15 ):
            alphax = 0.15
        if(alphay <0.15 ):
            alphay = 0.15
        alphax = np.sqrt(alphax)
        alphay = np.sqrt(alphay)
        alphaxy = 2.0*sigxy_calc
        #print (flux_calc, alphax, alphay, alphaxy, back_calc)
        
        counter += 1
        #print (counter)
        if (counter >= counter_target):
            break
    #print (counter) 
    #print (t1,t2,t3,t4,t5,t6,t7, alphax)             
    if(counter>= 98):
        return None, None, None,None, None, None, None, None, None, None
    else:
        return flux_calc, mux_calc, muy_calc, e1, e2, back_calc, np.sqrt(sigxx_calc + sigyy_calc), sigxx_calc, sigyy_calc, sigxy_calc


# =============================================================================
# fluxArr =np.array([30, 110, 220, 574])
# sigx = 4
# sigy = 4
# sigxy = 10
# bkg = 200
# sizex=sizey = 80
# 
# snrArr = fluxArr/ np.sqrt(3.14*sigx*sigy*bkg*4)
# measuredArr=[]
# measuredErrArr=[]
# correctedArr =[]
# for flux in fluxArr:
#     
#     tempArr=[]
#     for j in range(200):
#         
#         muArr= [sizex/2.0-0.5, sizey/2.0-0.5]
#         cov = [[(sigx**2), sigxy], [sigxy, (sigy**2)]]
#         const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
#         if(const<1):
#             const = 1
#         x, y = np.random.multivariate_normal(muArr, cov, const).T
#         x = np.int32(np.round(x))
#         y = np.int32(np.round(y))
#         obj = np.zeros((sizex,sizey))
#         np.add.at(obj, (y,x), 1)
#         
#         
#         noise = np.random.normal(bkg, np.sqrt(bkg), (sizex,sizey))
#         tot = np.add(obj,noise)
#         
#         flux_measure, mux_measure, muy_measure, e1_measure, e2, bkg_measure, psf_measure, sigxx_measure,sigyy_measure, sigxy_measure = measure(tot, sigx , sigy, sigxy, 1)
#         tempArr.append(flux_measure)
#         
#         
#     lut = np.load('/scratch/bell/dutta26/abell_2390/calib_combined.npy')
#     sqrtba = np.sqrt(bkg)*3.14*np.sqrt(sigx**2 *sigy**2 - sigxy**2)
#     if(flux<= 10):
#         index = int(round(flux/0.5))
#         
#     elif(flux<100):
#         index = int(round(flux/5)) + 18
#         
#     elif(flux<1000):
#         index = int(round(flux/50)) + 36
#         
#     elif(flux<10000):
#         index = int(round(flux/500)) + 54
#         
#     elif(flux<50000):
#         index = int(round(flux/5000)) + 72
#         
#     if(sqrtba <= 200):
#         baIndex = int(round((sqrtba-50)))*82
#     elif(sqrtba <= 1000):
#         baIndex = int(round((sqrtba-200)/5)*82) + 12300
#     elif(sqrtba <= 4950):
#         baIndex = int(round((sqrtba-1000)/30)*82) + 25420
#         
#     
#     if(lut[index][0]< flux):
#         nextIndex = index +1
#     else:
#         nextIndex = index - 1
#     if(index >  nextIndex):
#         index, nextIndex = nextIndex, index
#         
#     if(lut[baIndex][1]< sqrtba):
#         nextbaIndex = baIndex + 82
#     else:
#         nextbaIndex = baIndex - 82
#     if(baIndex > nextbaIndex):
#         baIndex, nextbaIndex = nextbaIndex, baIndex
#         
#     a1 = np.interp(sqrtba ,[lut[baIndex+index][1], lut[nextbaIndex+index][1]] , [lut[baIndex+index][2], lut[nextbaIndex+index][2]])
#     a2 = np.interp(sqrtba ,[lut[baIndex+nextIndex][1], lut[nextbaIndex+nextIndex][1]] , [lut[baIndex+nextIndex][2], lut[nextbaIndex+nextIndex][2]])
#     a3 = np.interp(flux,[lut[baIndex+index][0], lut[baIndex+nextIndex][0]], [a1,a2] )
#     
#     print (a3, index, nextIndex, baIndex, nextbaIndex)
#     
#     measuredArr.append(np.nanmedian(tempArr))
#     measuredErrArr.append(np.nanstd(tempArr))
#     correctedArr.append(np.nanmedian(tempArr) - a3)
#     
# 
# 
# plt.errorbar(np.log10(snrArr), measuredArr, yerr= measuredErrArr, fmt = '.k', label ='Measured')
# plt.errorbar(np.log10(snrArr), correctedArr, yerr= measuredErrArr, fmt = '.r', label ='Corrected')
# plt.errorbar(np.log10(snrArr), fluxArr, yerr= measuredErrArr, fmt = '.b',label ='True')
# plt.legend()
# plt.xlabel('Log SNR')
# plt.ylabel('Flux')
# plt.title('Elliptical Source')
# 
# =============================================================================

fluxArr =np.array([ 10000, 20000, 50000])
sigx = 3
sigy = 3
sigxy = 0
bkg = 200
sizex=sizey = 80
lut = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')

snrArr = fluxArr/ np.sqrt(3.14*sigx*sigy*bkg*4)
measuredArr_flux=[]
measuredErrArr_flux=[]
correctedArr_flux =[]

measuredArr_sigx=[]
measuredErrArr_sigx=[]
correctedArr_sigx =[]

measuredArr_sigxy=[]
measuredErrArr_sigxy=[]
correctedArr_sigxy =[]

for flux in fluxArr:
    
    tempArr_sigx=[]
    tempArr_sigxy=[]
    tempArr_flux=[]
    for j in range(200):
        
        muArr= [sizex/2.0-0.5, sizey/2.0-0.5]
        cov = [[(sigx**2), sigxy], [sigxy, (sigy**2)]]
        const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
        if(const<1):
            const = 1
        x, y = np.random.multivariate_normal(muArr, cov, const).T
        x = np.int32(np.round(x))
        y = np.int32(np.round(y))
        obj = np.zeros((sizex,sizey))
        np.add.at(obj, (y,x), 1)
        
        
        noise = np.random.normal(bkg, np.sqrt(bkg), (sizex,sizey))
        tot = np.add(obj,noise)
        
        flux_measure, mux_measure, muy_measure, e1_measure, e2, bkg_measure, psf_measure, sigxx_measure,sigyy_measure, sigxy_measure = measure(tot, sigx, sigy, sigxy, 1)
        tempArr_sigx.append(sigxx_measure)
        tempArr_sigxy.append(sigxy_measure)
        tempArr_flux.append(flux_measure)
        
        
    
    sqrtba = np.sqrt(bkg)*3.14*np.sqrt((sigx)**2 *sigy**2 - sigxy**2)
    fbysqrtba = np.log10(flux/sqrtba)
    print (fbysqrtba, 'aaaaa')
    if(fbysqrtba< -3.3):
        index = 1
    elif(fbysqrtba> 2.125):
        index = 216
    else:
        index = int((fbysqrtba + 3.3)/0.025)
    print (fbysqrtba, 'aaaaa', index)    
    a3 = np.interp(fbysqrtba,lut[index-1:index+2,0], lut[index-1:index+2,2] )
    
    print (a3)
    
    measuredArr_sigx.append(np.nanmedian(tempArr_sigx))
    measuredErrArr_sigx.append(np.nanstd(tempArr_sigx))
    correctedArr_sigx.append(np.nanmedian(tempArr_sigx) - a3*(sigx)**2)
    
    
    measuredArr_sigxy.append(np.nanmedian(tempArr_sigxy))
    measuredErrArr_sigxy.append(np.nanstd(tempArr_sigxy))
    correctedArr_sigxy.append(np.nanmedian(tempArr_sigxy) - sigxy*a3)
    
    a3 = np.interp(fbysqrtba,lut[index-1:index+2,0], lut[index-1:index+2,1] ) * sqrtba
    
    
    measuredArr_flux.append(np.nanmedian(tempArr_flux))
    measuredErrArr_flux.append(np.nanstd(tempArr_flux))
    correctedArr_flux.append(np.nanmedian(tempArr_flux) - a3)

e1 = (sigx**2 - sigy**2)/(sigx**2 +sigy**2)
e2 = 2*sigxy/(sigx**2 +sigy**2)


plt.figure(1)
plt.subplot(221)
plt.errorbar(np.log10(snrArr), measuredArr_sigx, yerr= measuredErrArr_sigx, fmt = '.k')
plt.errorbar(np.log10(snrArr), correctedArr_sigx, yerr= measuredErrArr_sigx, fmt = '.r')
plt.errorbar(np.log10(snrArr), sigx**2* np.ones(len(snrArr))/2, yerr= measuredErrArr_sigx, fmt = '.b')
plt.title('e1 , e2 = '+str(e1)[0:5] + ' , '+str(e2)[0:5])
plt.ylabel('Sigmaxx')

plt.subplot(222)
plt.errorbar(np.log10(snrArr), measuredArr_sigxy, yerr= measuredErrArr_sigxy, fmt = '.k')
plt.errorbar(np.log10(snrArr), correctedArr_sigxy, yerr= measuredErrArr_sigxy, fmt = '.r')
plt.errorbar(np.log10(snrArr), sigxy* np.ones(len(snrArr))/2, yerr= measuredErrArr_sigxy, fmt = '.b')
plt.ylabel('Sigmaxy')

plt.subplot(223)
plt.errorbar(np.log10(snrArr), np.log10(measuredArr_flux), yerr= np.log10(measuredErrArr_flux), fmt = '.k')
plt.errorbar(np.log10(snrArr), np.log10(correctedArr_flux), yerr= np.log10(measuredErrArr_flux), fmt = '.r')
plt.errorbar(np.log10(snrArr), np.log10(fluxArr), yerr= np.log10(measuredErrArr_flux), fmt = '.b')
plt.ylabel('Log Flux')




plt.legend()
#plt.xlabel('Log SNR')
#plt.ylabel('Sigxx')





