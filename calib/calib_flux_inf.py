#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  5 08:10:01 2023

@author: dutta26
"""


import sys
from astropy.io import fits
import numpy as np
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
                alphaxy = alphaxy - 0.1
            if (alphaxy < 0):
                alphaxy = alphaxy + 0.1
            if (abs(alphaxy) > 5000):
                counter = counter_target + 1
                break
        if(abs(alphaxy)> 50  or abs(alphax) > sizex/2 or abs(alphay)> sizey/2):
            counter = counter_target + 1
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
            counter = counter_target + 1
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
        
    #print (counter) 
    #print (t1,t2,t3,t4,t5,t6,t7, alphax)             
    if(counter == counter_target + 1):
        return None, None, None,None, None, None, None, None, None, None
    else:
        return flux_calc, mux_calc, muy_calc, e1, e2, back_calc, np.sqrt(sigxx_calc + sigyy_calc), sigxx_calc, sigyy_calc, sigxy_calc






id_num = str(sys.argv[1])

sizex = sizey = 80
sigxy = 0
fluxArr = np.hstack((np.arange(500, 1000, 100), 
                     np.arange(1000, 10000, 1000) , np.arange(10000, 110000, 10000)))
bkgArr= [ 50, 100, 250, 500, 750, 1000, 1500, 2000]
sigxArr= [ 3,  4,  5, 6]
if(id_num == '1004'):
    sigxArr= [ 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7 ]
    fluxArr = np.hstack((np.arange(30000, 110000, 10000)))
    bkgArr= [50, 60, 70]
totRows = len(fluxArr)*len(bkgArr)*len(sigxArr)
store = np.ones((totRows, 4))
store = store*np.nan

count = 0
for sigx in sigxArr:
    if(id_num == '1001'):
        sigy =sigx
        sigxy = 0
    if(id_num == '1002'):
        sigy = 4
        sigxy = 0
    if(id_num == '1003'):
        sigy = 4
        sigxy = 5
    if(id_num == '1004'):
        sigy = 2
        sigxy = 0
        
    for bkg in bkgArr:
        
        print (bkg, sigx, sigy)
        area = 3.14*np.sqrt(sigx**2*sigy**2 - sigxy**2)
        biasArr=[]
        biasArr2=[]
        for flux in fluxArr:
            
            temp_measured =[]
            temp_measured2=[]
            for j in range(500):
                muArr= [sizex/2.0-0.5, sizey/2.0-0.5]
                cov = [[(sigx**2), sigxy], [sigxy, (sigy**2)]]
                const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
                if(const < 1):
                    const = 1
                x, y = np.random.multivariate_normal(muArr, cov, const).T
                x = np.int32(np.round(x))
                y = np.int32(np.round(y))
                obj = np.zeros((sizex,sizey))
                np.add.at(obj, (y,x), 1)
                
                
                noise = np.random.normal(bkg, np.sqrt(bkg), (sizex,sizey))
                noise = np.round(noise)
                obj = obj + noise
                flux_calc, mux_calc, muy_calc, e1_calc, e2_calc, bkg_calc, psf_calc, sigxx_calc, sigyy_calc, sigxy_calc = measure(obj, sigx, sigy, sigxy, 100)
                if(flux_calc == None or sigxx_calc == None or np.isnan(flux_calc) or np.isnan(psf_calc) or psf_calc == None):
                    continue
                if(abs(sigxx_calc - sigx**2/2) > 100 or  abs(sigyy_calc - sigy**2/2) > 100):
                    continue
                if(sigxx_calc < 0 or sigxx_calc > 100 or sigyy_calc<0 or sigyy_calc>100):
                    continue
                
                temp_measured.append(flux_calc)
                temp_measured2.append(sigxx_calc)
            
            if(len(temp_measured) == 0):
                #print ('aa')
                store [count, 0], store [count, 1], store [count, 2], store [count, 3] = None, None, None, None
            else:
                biasArr.append(np.nanmedian(temp_measured) - flux)
                biasArr2.append( (np.nanmedian(temp_measured2) - sigx**2/2)/ sigx**2)
                store [count, 0], store [count, 1], store [count, 2], store [count, 3] = flux, np.sqrt(bkg)*area, np.nanmedian(temp_measured) - flux, (np.nanmedian(temp_measured2) - sigx**2/2)/ sigx**2
                #store2[count, 0], store2[count, 1], store2[count, 2] = flux, bkg*3.14*sigx*sigy, (np.nanmedian(temp_measured2) - sigx**2/2)/ sigx**2
            
            count += 1
        #biasArr= np.array(biasArr)
        #biasArr2= np.array(biasArr2)
        #plt.plot(np.log10(fluxArr), biasArr2, label = str(bkg*3.14*16)[0:5])
            
#plt.ylabel('Sigmaxx Bias')   
#plt.xlabel('Log True Flux')     
#plt.legend()
np.save('/scratch/bell/dutta26/abell_2390/calib_' +str(id_num)+'_inf.npy', store)    
#np.save('/scratch/bell/dutta26/abell_2390/sigma_calib_66.npy', store2)    










# =============================================================================
# id_num = str(sys.argv[1])
# #id_num = '1'
# import matplotlib.pyplot as plt
# from scipy.special import erf
# sizex = sizey= 80
# sigx = 6.25
# sigy = 6.25
# sigxy = 0
# #fluxArr = [10, 20, 50, 100, 250, 500, 750, 1000, 2500, 5000, 10000]
# #bkgArr= [ 50, 100, 200]
# fluxArr = np.hstack((np.arange(0, 10, 0.5), np.arange(10, 100, 5), np.arange(100, 1000, 50), 
#                      np.arange(1000, 10000, 500) , np.arange(10000, 50000, 5000)))
# #fluxArr= [15000]
# #bkgArr =   np.hstack((np.arange(50, 100, 2), np.arange(100, 1000, 20), np.arange(1000, 3200, 200)))
# #sqrtbaArr = np.hstack( (np.arange(50, 200, 1), np.arange(200, 1000, 5), np.arange(1000, 4950, 30)))
# if(id_num == '1'):
#     sqrtbaArr = np.arange(50, 200, 1)
# if(id_num == '2'):
#     sqrtbaArr = np.arange(200, 1000, 5)
# if(id_num == '3'):
#     sqrtbaArr = np.arange(1000, 4950, 30)
# 
# areaArr= np.arange(1.5, 6.23, 0.2)**2 * 3.14
# totRows = len(fluxArr)*len(sqrtbaArr)
# store = np.ones((totRows, 4))
# store = store*np.nan
# 
# count = 0
# for sqrtba in sqrtbaArr:
#     for area in areaArr:
#         #print (ba/area)
#         bkg = (sqrtba/area)**2
#         sigx = sigy = np.sqrt(area/3.14)
#         if(bkg >= 49 and bkg<=2000):
#             break
#         else:
#             continue
#     
#     #print (bkg, sigx, sigy)
#     biasArr=[]
#     biasArr2=[]
#     for flux in fluxArr:
#         
#         temp_measured =[]
#         temp_measured2=[]
#         for j in range(200):
#             muArr= [sizex/2.0-0.5, sizey/2.0-0.5]
#             cov = [[(sigx**2), sigxy], [sigxy, (sigy**2)]]
#             const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
#             if(const < 1):
#                 const = 1
#             x, y = np.random.multivariate_normal(muArr, cov, const).T
#             x = np.int32(np.round(x))
#             y = np.int32(np.round(y))
#             obj = np.zeros((sizex,sizey))
#             np.add.at(obj, (y,x), 1)
#             
#             
#             noise = np.random.normal(bkg, np.sqrt(bkg), (sizex,sizey))
#             noise = np.round(noise)
#             obj = obj + noise
#             flux_calc, mux_calc, muy_calc, e1_calc, e2_calc, bkg_calc, psf_calc, sigxx_calc, sigyy_calc, sigxy_calc = measure(obj, sigx, sigy, sigxy, 100)
#             if(flux_calc == None or sigxx_calc == None or np.isnan(flux_calc) or np.isnan(psf_calc) or psf_calc == None):
#                 continue
#             if(abs(sigxx_calc - sigx**2/2) > 100 or  abs(sigyy_calc - sigy**2/2) > 100):
#                 continue
#             if(sigxx_calc < 0 or sigxx_calc > 100 or sigyy_calc<0 or sigyy_calc>100):
#                 continue
#             
#             temp_measured.append(flux_calc)
#             temp_measured2.append(sigxx_calc)
#         
#         if(len(temp_measured) == 0):
#             store [count, 0], store [count, 1], store [count, 2], store [count, 3] = None, None, None, None
#         else:
#             biasArr.append(np.nanmedian(temp_measured) - flux)
#             biasArr2.append( (np.nanmedian(temp_measured2) - sigx**2/2)/ sigx**2)
#             store [count, 0], store [count, 1], store [count, 2], store [count, 3] = flux, np.sqrt(bkg)*3.14*sigx*sigy, np.nanmedian(temp_measured) - flux, (np.nanmedian(temp_measured2) - sigx**2/2)/ sigx**2
#             #store2[count, 0], store2[count, 1], store2[count, 2] = flux, bkg*3.14*sigx*sigy, (np.nanmedian(temp_measured2) - sigx**2/2)/ sigx**2
#         
#         count += 1
#     #biasArr= np.array(biasArr)
#     #biasArr2= np.array(biasArr2)
#     #plt.plot(np.log10(fluxArr), biasArr2, label = str(bkg*3.14*16)[0:5])
#         
# #plt.ylabel('Sigmaxx Bias')   
# #plt.xlabel('Log True Flux')     
# #plt.legend()
# np.save('/scratch/bell/dutta26/abell_2390/calib_' +str(id_num)+'inf.npy', store)    
# #np.save('/scratch/bell/dutta26/abell_2390/sigma_calib_66.npy', store2)    
# 
# 
# 
# 
# 
# 
# 
# =============================================================================


















# =============================================================================
# import matplotlib.pyplot as plt
# from scipy.special import erf
# sizex = sizey= 80
# sigx = 2
# sigy = 2
# sigxy = 0
# fluxArr = [8000]
# bkgArr= [ 800]
# 
# 
# count = 0
# for bkg in bkgArr:
#     
#     print (bkg, sigx, sigy)
#     biasArr=[]
#     biasArr2=[]
#     for flux in fluxArr:
#         
#         temp_measured =[]
#         temp_measured2=[]
#         for j in range(500):
#             muArr= [sizex/2.0-0.5, sizey/2.0-0.5]
#             cov = [[(sigx**2), sigxy], [sigxy, (sigy**2)]]
#             const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
#             if(const < 1):
#                 const = 1
#             x, y = np.random.multivariate_normal(muArr, cov, const).T
#             x = np.int32(np.round(x))
#             y = np.int32(np.round(y))
#             obj = np.zeros((sizex,sizey))
#             np.add.at(obj, (y,x), 1)
#             
#             
#             noise = np.random.normal(bkg, np.sqrt(bkg), (sizex,sizey))
#             noise = np.round(noise)
#             obj = obj + noise
#             flux_calc, mux_calc, muy_calc, e1_calc, e2_calc, bkg_calc, psf_calc, sigxx_calc, sigyy_calc, sigxy_calc = measure(obj, sigx, sigy, sigxy, 1,0)
#             
#             temp_measured.append(flux_calc)
#             temp_measured2.append(sigxx_calc)
#             
#         biasArr.append((np.nanmedian(temp_measured) - flux)*np.sqrt(bkg))
#         biasArr2.append( (np.nanmedian(temp_measured2) - sigx**2/2)/ sigx**2)
#         
#     #biasArr= np.array(biasArr)
#     #biasArr2= np.array(biasArr2)
#     #plt.plot(np.log10(fluxArr), biasArr2, label = str(bkg*3.14*16)[0:5])
#         
# #plt.ylabel('Sigmaxx Bias')   
# #plt.xlabel('Log True Flux')     
# #plt.legend()
# print (biasArr , biasArr2)
# =============================================================================
