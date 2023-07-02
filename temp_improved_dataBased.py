#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 07:39:45 2023

@author: dutta26
"""

import numpy as np
# =============================================================================
# sigxArr = [2.5,3,4,5,6]
# sigyArr = [2.5,3,4,5,6]
# =============================================================================

sigxArr = [7]
sigyArr = [2.9167,3.5, 7]

modeArr = [1,2]
sizex =sizey = 100
iterations = 1000
store = np.zeros ((300, 30), dtype = np.float32 )
cnt = 0
for mode in modeArr:
    if(mode == 1):
        fluxArr=[50, 500, 10000, 100000, 2e6]
    if(mode == 2):
        fluxArr=[6000, 12000, 15000, 100000, 2e6]
        
        
    for sigx in sigxArr:
        print (sigx)
        for sigy in sigyArr:
            if(sigx< sigy):
                continue 
            
            
            if(sigx == sigy):
                thetaArr=[0]
            else:
                thetaArr= [0, 45, 90]
            for theta in thetaArr:
                #print (sigx, sigy, theta)
                sigxx = ((sigx * np.cos(theta*np.pi/180))**2  + (sigy * np.sin(theta*np.pi/180))**2)
                sigyy = ((sigx * np.sin(theta*np.pi/180))**2  + (sigy * np.cos(theta*np.pi/180))**2)
                sigxy = (sigx**2 - sigy**2)*np.sin(theta*np.pi/180)*np.cos(theta*np.pi/180)
                #print (sigxx, sigyy, sigxy)
                for flux in fluxArr:
                    if(mode == 1):
                        bkgArr=[0]
                    else:
                        bkgArr= [50, 200, 400]
                    for bkg in bkgArr:
                        
                        #Intrinsic component
                        if(mode == 1):
                            
                            perfectShape = 0
                        #Bkg component
                        else:
                            
                            perfectShape = 1
                            
                        fixSigmas = 1 #1 for yes and 0 for no
                        breakAt1 = 1 #1 for yes and 0 for no
                        
                        #Error arrays
                        fluxErrArr=[]
                        sizeErrArr=[]
                        muErrArr=[]
                        xxErrArr=[]
                        yyErrArr=[]
                        xyErrArr=[]
                        e1ErrArr=[]
                        e2ErrArr=[]
                        ellipErrArr=[]
                        
                        #Error arrays
                        fluxErrArr2=[]
                        sizeErrArr2=[]
                        muErrArr2=[]
                        xxErrArr2=[]
                        yyErrArr2=[]
                        xyErrArr2=[]
                        e1ErrArr2=[]
                        e2ErrArr2=[]
                        ellipErrArr2=[]
                        
                        for trials in range(5):
                            for j in range(iterations):
                                num=den=0
                                if(perfectShape == 0):
                                    #Make real object by sampling. 
                                    
                                    muArr= [sizex/2.0-0.5, sizey/2.0-0.5]
                                    cov = [[sigxx, sigxy], [sigxy, sigyy]]
                                    const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
                                    x, y = np.random.multivariate_normal(muArr, cov, const).T
                                    x = np.int32(np.round(x))
                                    y = np.int32(np.round(y))
                                    obj = np.zeros((sizex,sizey))
                                    np.add.at(obj, (y,x), 1)
                                
                        
                                elif(perfectShape == 1):
                                    
                                    x = np.linspace(0, sizex-1, sizex)
                                    y = np.linspace(0, sizey-1, sizey)
                                    x = x-sizex/2.0+0.5
                                    y = y-sizey/2.0+0.5
                                    x, y = np.meshgrid(x, y)
                                    
                                    alphax1 = np.sqrt(sigxx)
                                    alphay1 = np.sqrt(sigyy)
                                    alphaxy1 = sigxy
                                    A1 =1/(2*np.pi*alphax1*alphay1*np.sqrt(1- (alphaxy1/(alphax1*alphay1))**2 ))
                                    arb_const1 = 2*(1 - (alphaxy1/ (alphax1*alphay1))**2 )
                                    k1=(A1 * np.exp(-((x)**2/(arb_const1*alphax1**2)+ (y)**2/(arb_const1*alphay1**2)  - 2*alphaxy1*(y)*(x)/(arb_const1*alphax1**2 * alphay1**2 ) )))
                                    
                                    obj = k1/np.sum(k1)
                                    const = int(round(np.random.normal(flux, np.sqrt(flux))))
                                    obj = obj*const
                                
                                
                                
                                #Use integer background values
                                noise = np.random.normal(bkg, np.sqrt(bkg), (sizex,sizey))
                                
                                #End make real object and background
                                
                                x = np.linspace(0, sizex-1, sizex)
                                y = np.linspace(0, sizey-1, sizey)
                                x= x -sizex/2.0 + 0.5 
                                y= y -sizey/2.0 + 0.5 
                                #mux=muy= np.random.normal(29.5, 0.15)
                                x, y = np.meshgrid(x, y)
                                #Add noise and image 
                                tot = np.add(obj,noise)
                                
                                
                                #z = k*(np.random.normal(0, np.sqrt(bkg), (sizex,sizey)) )
                                delSigxx = 999
                                delSigyy = 999
                                delSigxy = 999
                                prevSigxx = 9999
                                prevSigyy = 9999
                                alphax = 3
                                alphay = 3
                                alphaxy =0
                                sixx_calc =0
                                sigyy_calc=0
                                sigxy_calc = 0
                                mux_calc = 0
                                muy_calc =0
                                flux_calc = 0 
                                e1 = 0.0
                                e2 = 0.0
                                med = np.median(tot)
                                back_calc = med
                                #back_calc = bkg
                                total = np.sum(tot)
                                counter = 0
                                #while(abs(delSigxx)>0.0001 or abs(delSigyy)> 0.0001  ):
                                while(counter < 100):
                                    if( fixSigmas == 1):
                                        alphax=np.sqrt(sigxx)
                                        alphay=np.sqrt(sigyy)
                                        alphaxy = sigxy
                                    
                                    A =1/(2*np.pi*alphax*alphay*np.sqrt(1- (alphaxy/(alphax*alphay))**2 ))
                                    arb_const = 2*(1 - (alphaxy/ (alphax*alphay))**2 )
                                    k=(A * np.exp(-((x-mux_calc)**2/(arb_const*alphax**2)+ (y-muy_calc)**2/(arb_const*alphay**2)  - 2*alphaxy*(y-muy_calc)*(x-mux_calc)/(arb_const*alphax**2 * alphay**2 ) )))
                                    
                                    
                                    t1= np.sum(x*y* (tot - back_calc) * k)
                                    t2 = np.sum(k*(tot - back_calc))
                                    t3 = np.sum(x*(tot- back_calc)*k)
                                    t4 = np.sum(y*(tot- back_calc)*k)
                                    t5 = np.sum(x * x * (tot - back_calc) * k)
                                    t6 = np.sum(y * y * (tot - back_calc) * k)
                                    t7 = np.sum(k**2)
                                    
                                    mux_calc = t3/t2
                                    muy_calc = t4/t2
                                    flux_calc = t2/t7
                                    
                                    
                                    sigxy_calc = (t1/t2) - (t3*t4)/(t2*t2)
                                    sigxx_calc = (t5/t2) - (t3*t3) / (t2*t2)
                                    sigyy_calc = (t6/t2) - (t4*t4) / (t2*t2)
                                    
                                    
                                    e1 = (sigxx_calc - sigyy_calc)/(sigxx_calc + sigyy_calc)
                                    e2 = 2*sigxy_calc/(sigxx_calc + sigyy_calc)
                                    ellip = np.sqrt(e1*e1 + e2*e2)
                                    
                                    if(med != 0 ):
                                        back_calc = (total - flux_calc)/ (sizex*sizey)
                                        #total = flux_calc + (sizex*sizey* back_calc)
                                    #back_calc = bkg
                                    
                                    delSigxx = prevSigxx - sigxx_calc
                                    delSigyy = prevSigyy - sigyy_calc
                                    
                                    prevSigxx = sigxx_calc
                                    prevSigyy = sigyy_calc
                                    
                                    alphax = sigxx_calc*2.0
                                    alphay = sigyy_calc*2.0
                                    if(alphax <0.15 ):
                                        alphax = 0.15
                                    if(alphay <0.15 ):
                                        alphay = 0.15
                                    alphax = np.sqrt(alphax)
                                    alphay = np.sqrt(alphay)
                                    alphaxy = 2.0*sigxy_calc
                                    counter += 1
                                    
                                    if(np.abs(delSigxx) <0.001 and np.abs(delSigyy)<0.001):
                                        break
                                    
                                    if(breakAt1 == 1):
                                        break
                                sigxx_th = 0.5*sigxx
                                sigyy_th = 0.5*sigyy
                                sigxy_th = 0.5*sigxy
                                e1_th = (sigxx_th - sigyy_th)/(sigxx_th + sigyy_th)
                                e2_th = 2*sigxy_th/(sigxx_th + sigyy_th)
                                ellip_th = np.sqrt(e1_th*e1_th + e2_th*e2_th)
                                fluxErrArr2.append(flux-flux_calc)
                                sizeErrArr2.append(np.sqrt(sigxx_th+sigyy_th) - np.sqrt(sigxx_calc+sigyy_calc) )
                                muErrArr2.append(np.sqrt( (mux_calc)**2 +(muy_calc)**2  ))
                                xxErrArr2.append(sigxx_th - sigxx_calc)
                                yyErrArr2.append(sigyy_th - sigyy_calc)
                                xyErrArr2.append(sigxy_th - sigxy_calc)
                                e1ErrArr2.append(e1- e1_th)
                                e2ErrArr2.append(e2- e2_th)
                                ellipErrArr2.append(ellip- ellip_th)
                                
                            fluxErrArr.append(np.std(fluxErrArr2))
                            sizeErrArr.append(np.std(sizeErrArr2))
                            muErrArr.append(np.std(muErrArr2))
                            xxErrArr.append(np.std(xxErrArr2))
                            yyErrArr.append(np.std(yyErrArr2))
                            xyErrArr.append(np.std(xyErrArr2))
                            e1ErrArr.append(np.std(e1ErrArr2))
                            e2ErrArr.append(np.std(e2ErrArr2))
                            ellipErrArr.append(np.std(ellipErrArr2))
                        
# =============================================================================
#                         print(flux, bkg, sigx, sigy, theta, sigxx, sigyy, sigxy, mode)    
#                         print (np.mean(fluxErrArr), np.mean(sizeErrArr), np.mean(muErrArr))
#                         print (np.mean(xxErrArr), np.mean(yyErrArr), np.mean(xyErrArr))
#                         print (np.mean(e1ErrArr), np.mean(e2ErrArr), np.mean(ellipErrArr))
#                         print (bkg, mode)
#                         print ('************************')
# =============================================================================
                        store[cnt,0:9] = flux, bkg, sigx, sigy, theta, sigxx, sigyy, sigxy, mode
                        
                        #Store errors
                        store[cnt, 9:12] = np.mean(fluxErrArr), np.mean(sizeErrArr), np.mean(muErrArr)
                        store[cnt, 12:15] = np.mean(xxErrArr), np.mean(yyErrArr), np.mean(xyErrArr)
                        store[cnt, 15:18] = np.mean(e1ErrArr), np.mean(e2ErrArr), np.mean(ellipErrArr)
                        
                        #Store error of errors
                        store[cnt, 18:21] = np.std(fluxErrArr), np.std(sizeErrArr), np.std(muErrArr)
                        store[cnt, 21:24] = np.std(xxErrArr), np.std(yyErrArr), np.std(xyErrArr)
                        store[cnt, 24:27] = np.std(e1ErrArr), np.std(e2ErrArr), np.std(ellipErrArr)
                        
                        cnt += 1
                        
np.save('/scratch/bell/dutta26/abell_2390/measure_test_data2m_supp.npy', store)
                        
                            
                            
                            