#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 10:07:24 2023

@author: dutta26
"""
import numpy as np
def measure(img ,lut1, lut2,guessFlux = 100, guessmux = 0, guessmuy=0 , guessAlphax =3, guessAlphay =3, guessAlphaxy =0 , counter_target = 100, distort = 0):
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
    if(sizex< 30 or sizey<30):
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
    back_calc = np.median(img)
    total = np.sum(img)
    counter = 0
    
    
    
    #Loop until convergence
    
    while( (abs(delSigxx)>0.001 or abs(delSigyy)> 0.001) and counter<counter_target):
        
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
        
        if(med != 0 ):
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





# =============================================================================
# import measure_pythonV
# import matplotlib.pyplot as plt
# fluxArr =np.array([1700,1800, 19000]) 
# sigx = 5
# sigy = 5
# sigxy = 5
# bkg = 200
# sizex=sizey = 60
# lut = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
# 
# snrArr = fluxArr/ np.sqrt(3.14*bkg*4*np.sqrt(sigx**2*sigy**2 - sigxy**2))
# measuredArr_flux=[]
# measuredErrArr_flux=[]
# correctedArr_flux =[]
# 
# measuredArr_sigx=[]
# measuredErrArr_sigx=[]
# correctedArr_sigx =[]
# 
# measuredArr_sigxy=[]
# measuredErrArr_sigxy=[]
# correctedArr_sigxy =[]
# 
# for flux in fluxArr:
#     
#     tempArr_sigx=[]
#     tempArr_sigxy=[]
#     tempArr_flux=[]
#     for j in range(1000):
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
#         flux_measure, mux_measure, muy_measure, e1_measure, e2, bkg_measure, psf_measure, sigxx_measure,sigyy_measure, sigxy_measure = measure(tot, flux, 0 ,0, sigx, sigy, sigxy, 1, 1)
#         #flux_measure, mux_measure, muy_measure, e1_measure, e2, bkg_measure, psf_measure, sigxx_measure,sigyy_measure, sigxy_measure = measure_pythonV.measure_v2(tot, sigx, sigy, sigxy, 100)
# 
#         tempArr_sigx.append((sigxx_measure- sigyy_measure) / (sigxx_measure+ sigyy_measure))
#         tempArr_sigxy.append(2*sigxy_measure/ (sigxx_measure+ sigyy_measure))
#         tempArr_flux.append(flux_measure)
#         
#         
#    
#     
#     measuredArr_sigx.append(np.nanmean(tempArr_sigx))
#     measuredErrArr_sigx.append(np.nanstd(tempArr_sigx))
#     #correctedArr_sigx.append(np.nanmedian(tempArr_sigx) - a3*(sigx)**2)
#     
#     
#     measuredArr_sigxy.append(np.nanmean(tempArr_sigxy))
#     measuredErrArr_sigxy.append(np.nanstd(tempArr_sigxy))
#     #correctedArr_sigxy.append(np.nanmedian(tempArr_sigxy) - sigxy*a3)
#     
#     
#     
#     measuredArr_flux.append(np.nanmedian(tempArr_flux))
#     measuredErrArr_flux.append(np.nanstd(tempArr_flux))
#     #correctedArr_flux.append(np.nanmedian(tempArr_flux) - a3)
# 
# e1 = (sigx**2 - sigy**2)/(sigx**2 +sigy**2)
# e2 = 2*sigxy/(sigx**2 +sigy**2)
# 
# 
# plt.figure(1)
# plt.subplot(221)
# plt.errorbar(np.log10(snrArr), measuredArr_sigx, yerr= measuredErrArr_sigx, fmt = '.k')
# #plt.errorbar(np.log10(snrArr), correctedArr_sigx, yerr= measuredErrArr_sigx, fmt = '.r')
# #plt.errorbar(np.log10(snrArr), sigx**2* np.ones(len(snrArr))/2, yerr= measuredErrArr_sigx, fmt = '.b')
# plt.errorbar(np.log10(snrArr), e1* np.ones(len(snrArr)), yerr= measuredErrArr_sigx, fmt = '.b')
# plt.title('e1 , e2 = '+str(e1)[0:5] + ' , '+str(e2)[0:5])
# plt.ylabel('Sigmaxx')
# 
# plt.subplot(222)
# plt.errorbar(np.log10(snrArr), measuredArr_sigxy, yerr= measuredErrArr_sigxy, fmt = '.k')
# #plt.errorbar(np.log10(snrArr), correctedArr_sigxy, yerr= measuredErrArr_sigxy, fmt = '.r')
# #plt.errorbar(np.log10(snrArr), sigxy* np.ones(len(snrArr))/2, yerr= measuredErrArr_sigxy, fmt = '.b')
# plt.errorbar(np.log10(snrArr), e2* np.ones(len(snrArr)), yerr= measuredErrArr_sigxy, fmt = '.b')
# 
# plt.ylabel('Sigmaxy')
# 
# plt.subplot(223)
# plt.errorbar(np.log10(snrArr), np.log10(measuredArr_flux), yerr= np.log10(measuredErrArr_flux), fmt = '.k')
# #plt.errorbar(np.log10(snrArr), np.log10(correctedArr_flux), yerr= np.log10(measuredErrArr_flux), fmt = '.r')
# plt.errorbar(np.log10(snrArr), np.log10(fluxArr), yerr= np.log10(measuredErrArr_flux), fmt = '.b')
# plt.ylabel('Log Flux')
# 
# 
# 
# 
# plt.legend()
# #plt.xlabel('Log SNR')
# #plt.ylabel('Sigxx')
# 
# =============================================================================

