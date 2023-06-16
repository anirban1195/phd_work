#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 17:14:06 2023

@author: anirban
"""
import numpy as np
def measure(img ,lut1, lut2, guessFlux = 100, guessmux = 0, guessmuy=0 , guessAlphax =3, guessAlphay =3, guessAlphaxy =0 , counter_target = 100, distort = 0):
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
        #print (alphax, alphay, alphaxy)
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
        k = 1
        img1 = img
        
        t1= np.sum(x*y* (img1) * k)
        t2 = np.sum(k*(img1))
        t3 = np.sum(x*(img1)*k)
        t4 = np.sum(y*(img1)*k)
        t5 = np.sum(x * x * (img1) * k)
        t6 = np.sum(y * y * (img1) * k)
        t7 = np.sum(k**2)
        
        t8 = np.sum(x*x*x* (img1) * k)
        t9 = np.sum(y*x*x* (img1) * k)
        t10 = np.sum(y*y*x* (img1) * k)
        t11 = np.sum(y*y*y* (img1) * k)
        
        t12 = np.sum(x*x*x*x* (img1) * k)
        t13 = np.sum(y*x*x*x* (img1) * k)
        t14 = np.sum(y*y*x*x* (img1) * k)
        t15 = np.sum(y*y*y*x* (img1) * k)
        t16 = np.sum(y*y*y*y* (img1) * k)
        #print (t1,t2,t3,t4,t5,t6, t7)
        mux_calc = t3/t2
        muy_calc = t4/t2
        flux_calc = t2/t7
        
        
        sigxy_calc = (t1/t2) - (t3*t4)/(t2*t2)
        sigxx_calc = (t5/t2) - (t3*t3) / (t2*t2)
        sigyy_calc = (t6/t2) - (t4*t4) / (t2*t2)
        
        sigx4_calc = t12/t2 - 4*(mux_calc)*(t8/t2) + 6*(mux_calc)**2*(t5/t2) -3*(mux_calc)**4
        sigy4_calc = t16/t2 - 4*(muy_calc)*(t11/t2) + 6*(muy_calc)**2*(t6/t2) -3*(muy_calc)**4
        sigx2y2_calc = t14/t2 - 2*(muy_calc)*(t9/t2) - 2*(mux_calc)*(t10/t2) + (muy_calc)**2 * (t5/t2) + (mux_calc)**2*(t6/t2) + 4*mux_calc*muy_calc*(t1/t2) - 3*mux_calc**2*muy_calc**2
        sigx3y_calc = t13/t2 - (muy_calc)*(t8/t2) - 3*(mux_calc)*(t9/t2) + 3*mux_calc*muy_calc*(t5/t2) + 3*(mux_calc)**2 * (t1/t2) - 3*mux_calc**3 *muy_calc
        sigxy3_calc = t15/t2 - (mux_calc)*(t11/t2) - 3*(muy_calc)*(t10/t2) + 3*muy_calc*mux_calc*(t6/t2) + 3*(muy_calc)**2 * (t1/t2) - 3*muy_calc**3*mux_calc
        
        #If negative sigmas, break
        if(sigxx_calc<0 or sigyy_calc<0):
            counter = counter_target+999
            break
        
      
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
        break
    
    #print (t9/t2,t5/t2, sigxx_calc)
    #avg = (t9/t2)/((t5/t2)**2 )
    #avg = ((t9/t2)+ 4*(t10/t2) + 6*(t11/t2) + 4*(t12/t2) +  (t13/t2) )/  ((t5/t2)+ (t6/t2) +2*(t1/t2))**2
    #avg = ((t9/t2)+  2*(t11/t2)  +  (t13/t2) )/  ((t5/t2)+ (t6/t2) )**2
    #print (sigx3y_calc, sigxy3_calc)
    k1 = (sigx4_calc + 4*sigx3y_calc + 6*sigx2y2_calc + 4*sigxy3_calc + sigy4_calc)/ (sigxx_calc + sigyy_calc + 2*sigxy_calc)**2
    k2 = (sigx4_calc + 2*sigx2y2_calc + sigy4_calc)/ (sigxx_calc + sigyy_calc)**2
    return k1,k2
    

