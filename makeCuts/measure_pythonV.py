#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 21:24:54 2021

@author: dutta26
"""

import numpy as np
from astropy.io import fits
def measure(img):
    #Shape the image properly
    img = np.array(img)
    sizey, sizex = np.shape(img)
    #print (sizex,sizey)
    
    if(sizex > 100):
        midX = int(round(sizex/2.0))
        img = img[ :, midX-50: midX+50]
    if(sizey > 100):
        midY = int(round(sizey/2.0))
        img = img[midY-50: midY+50, :]
    if(sizex< 30 or sizey<30):
        return None, None, None,None, None, None, None
    hdu = fits.PrimaryHDU(img)
    hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/cpp_temp.fits', overwrite=True)
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
    alphax = 6
    alphay = 6
    alphaxy =0
    sigxx_calc =0
    sigyy_calc=0
    sigxy_calc = 0
    mux_calc = 0
    muy_calc =0
    flux_calc = 0 
    e1 = 0.0
    e2 = 0.0
    med = np.median(img)
    back_calc = med
    total = np.sum(img)
    counter = 0
    while(abs(delSigxx)>0.001 and abs(delSigyy)> 0.001 and counter<100):
        #print (counter)
        while( (alphaxy/(alphax*alphay))**2  >= 1):
            #print (alphaxy, (alphaxy/(alphax*alphay))**2 )
            if(alphaxy > 0 ):
                alphaxy = alphaxy - 0.101235
            if (alphaxy < 0):
                alphaxy = alphaxy + 0.101235
            if (abs(alphaxy) > 5000):
                counter = 100
                break
        if(abs(alphaxy)> 50  or abs(alphax) > sizex/4 or abs(alphay)> sizey/4):
            counter = 100
            break
        arbConst = 2*(1- (alphaxy/(alphax*alphay))**2 )  
        A =1/(2*np.pi*alphax*alphay*np.sqrt(1- (alphaxy/(alphax*alphay))**2 ))
        k=(A * np.exp(-((x-mux_calc)**2/(arbConst*alphax**2)+ (y-muy_calc)**2/(arbConst*alphay**2) - 2*alphaxy*(y-muy_calc)*(x-mux_calc)/(arbConst* alphax**2 * alphay**2 ) )))
        
        
        t1= np.sum(x*y* (img - back_calc) * k)
        t2 = np.sum(k*(img - back_calc))
        t3 = np.sum(x*(img- back_calc)*k)
        t4 = np.sum(y*(img- back_calc)*k)
        t5 = np.sum(x * x * (img - back_calc) * k)
        t6 = np.sum(y * y * (img - back_calc) * k)
        t7 = np.sum(k**2)
        
        mux_calc = t3/t2
        muy_calc = t4/t2
        flux_calc = t2/t7
        
        
        sigxy_calc = (t1/t2) - (t3*t4)/(t2*t2)
        sigxx_calc = (t5/t2) - (t3*t3) / (t2*t2)
        sigyy_calc = (t6/t2) - (t4*t4) / (t2*t2)
        
        if(sigxx_calc<0 or sigyy_calc<0):
            counter = 100
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
    if(counter>= 98):
        return None, None, None,None, None, None, None
    else:
        return flux_calc, mux_calc, muy_calc, e1, e2, back_calc, np.sqrt(sigxx_calc + sigyy_calc)


def measure_v2(img, guessAlphax =3, guessAlphay =3, guessAlphaxy =0 , counter_target = 100):
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
    back_calc = med 
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
        
        
        t1= np.sum(x*y* (img - back_calc) * k)
        t2 = np.sum(k*(img - back_calc))
        t3 = np.sum(x*(img- back_calc)*k)
        t4 = np.sum(y*(img- back_calc)*k)
        t5 = np.sum(x * x * (img - back_calc) * k)
        t6 = np.sum(y * y * (img - back_calc) * k)
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
    if(counter>= 98):
        return None, None, None,None, None, None, None, None, None, None
    else:
        return flux_calc, mux_calc, muy_calc, e1, e2, back_calc, np.sqrt(sigxx_calc + sigyy_calc), sigxx_calc, sigyy_calc, sigxy_calc


def measureReturnGuess(img):
    #Shape the image properly
    img = np.array(img)
    sizey, sizex = np.shape(img)
    #print (sizex,sizey)
    
    if(sizex > 100):
        midX = int(round(sizex/2.0))
        img = img[: , midX-50: midX+50]
    if(sizey > 100):
        midY = int(round(sizey/2.0))
        img = img[midY-50: midY+50 , :]
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
    alphax = 6
    alphay = 6
    alphaxy =0
    sigxx_calc =0
    sigyy_calc=0
    sigxy_calc = 0
    mux_calc = 0
    muy_calc =0
    flux_calc = 0 
    e1 = 0.0
    e2 = 0.0
    med = np.median(img)
    back_calc = med
    total = np.sum(img)
    counter = 0
    while(abs(delSigxx)>0.001 and abs(delSigyy)> 0.001 and counter<100):
        #print (counter)
        while( (alphaxy/(alphax*alphay))**2  >= 1):
            #print (alphaxy, (alphaxy/(alphax*alphay))**2 )
            if(alphaxy > 0 ):
                alphaxy = alphaxy - 0.101235
            if (alphaxy < 0):
                alphaxy = alphaxy + 0.101235
            if (abs(alphaxy) > 5000):
                counter = 100
                break
        if(abs(alphaxy)> 50  or abs(alphax) > sizex/4 or abs(alphay)> sizey/4):
            counter = 100
            break
        arbConst = 2*(1- (alphaxy/(alphax*alphay))**2 )  
        A =1/(2*np.pi*alphax*alphay*np.sqrt(1- (alphaxy/(alphax*alphay))**2 ))
        k=(A * np.exp(-((x-mux_calc)**2/(arbConst*alphax**2)+ (y-muy_calc)**2/(arbConst*alphay**2) - 2*alphaxy*(y-muy_calc)*(x-mux_calc)/(arbConst* alphax**2 * alphay**2 ) )))
        
        t1= np.sum(x*y* (img - back_calc) * k)
        t2 = np.sum(k*(img - back_calc))
        t3 = np.sum(x*(img- back_calc)*k)
        t4 = np.sum(y*(img- back_calc)*k)
        t5 = np.sum(x * x * (img - back_calc) * k)
        t6 = np.sum(y * y * (img - back_calc) * k)
        t7 = np.sum(k**2)
        
        mux_calc = t3/t2
        muy_calc = t4/t2
        flux_calc = t2/t7
        
        
        sigxy_calc = (t1/t2) - (t3*t4)/(t2*t2)
        sigxx_calc = (t5/t2) - (t3*t3) / (t2*t2)
        sigyy_calc = (t6/t2) - (t4*t4) / (t2*t2)
        #print (t1,t2,t3,t4)
        #print (sigxx_calc, sigyy_calc, sigxy_calc)
        if(sigxx_calc<0 or sigyy_calc<0):
            counter = 100
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
        #if (counter == 1):
        #    break
    #print (counter)        
    if(counter>= 98):
        return None, None, None,None, None, None, None, None, None, None
    else:
        return flux_calc, mux_calc, muy_calc, e1, e2, back_calc, np.sqrt(sigxx_calc + sigyy_calc), sigxx_calc, sigyy_calc, sigxy_calc



def measure_tver(img, guessAlphax =3, guessAlphay =3, guessAlphaxy =0 , counter_target = 1):
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
    back_calc = med 
    total = np.sum(img)
    counter = 0
    t1=t2=t3=t4=t5=t6=t7=0
    
    
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
        
        back_calc = np.median(img)
        t1= np.sum(x*y* (img - back_calc) * k)
        t2 = np.sum(k*(img - back_calc))
        t3 = np.sum(x*(img- back_calc)*k)
        t4 = np.sum(y*(img- back_calc)*k)
        t5 = np.sum(x * x * (img - back_calc) * k)
        t6 = np.sum(y * y * (img - back_calc) * k)
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
    
    return t1,t2,t3,t4,t5,t6,t7


def measure_v2T(img, guessAlphax =3, guessAlphay =3, guessAlphaxy =0 , counter_target = 100):
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
    back_calc = med 
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
        
        t1= np.sum(x*y* (img - back_calc) * k)
        t2 = np.sum(k*(img - back_calc))
        t3 = np.sum(x*(img- back_calc)*k)
        t4 = np.sum(y*(img- back_calc)*k)
        t5 = np.sum(x * x * (img - back_calc) * k)
        t6 = np.sum(y * y * (img - back_calc) * k)
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
    
    return t1,t2,t3,t4,t5,t6,t7

# =============================================================================
# #Testing the code
# bkg = 10
# num =0 
# sigx=6
# sigy=6
# mux=muy= 0
# sizex=sizey=180
# 
# aList=[]
# flux = 5000.0
# delBkg =0
# delmu = 0
# delSig =0
# muArr= [sizex/2.0-0.5+np.random.normal(0,delmu), sizey/2.0-0.5+np.random.normal(0,delmu)]
# cov = [[(sigx**2)+np.random.normal(0,delSig), 0], [0, (sigy**2)+np.random.normal(0,delSig)]]
# const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
# print (const)
# # =============================================================================
# # x, y = np.random.multivariate_normal(muArr, cov, const).T
# # x = np.int32(np.round(x))
# # y = np.int32(np.round(y))
# # obj = np.zeros((sizex,sizey))
# # np.add.at(obj, (x,y), 1)
# # =============================================================================
# 
# new_sigxSq = (sigx**2)+np.random.normal(0,delSig)
# new_sigySq = (sigy**2)+np.random.normal(0,delSig)
# x = np.linspace(0, sizex-1, sizex)
# y = np.linspace(0, sizey-1, sizey)
# x, y = np.meshgrid(x, y)
# obj = np.exp(-((x-sizex/2.0+0.5+np.random.normal(0,delmu))**2/(2*new_sigxSq) +(y-sizex/2.0+0.5+np.random.normal(0,delmu))**2/(2*new_sigySq)))
# obj = obj/np.sum(obj)
# obj = obj*const
# 
# 
# newBkg=np.random.normal(bkg,delBkg)
# noise = np.random.normal(newBkg, np.sqrt(newBkg), (sizex,sizey))
# noise = np.round(noise)
# 
# obj = np.add(obj,noise)
# obj = obj[40:140, 50:120]
# print (measure(obj))
# =============================================================================
