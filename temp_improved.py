#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 12:50:25 2020

@author: anirban
"""

import numpy as np

#bkgArr = np.arange(20,30, 10)3
bkgArr = [0]
num =0 
area_gamma = 36
sigx=5
sigy=3.5
sigxy = 0
mux=muy= 0
sizex=sizey=100

aList=[]
flux = 2000
delBkg =0
delmu = 0
delSig =0
mux_calc_arr=[]

d11 =[]
d22=[]

fixSigmas = 1 #1 for yes and 0 for no
breakAt1 = 1 #1 for yes and 0 for no
perfectShape = 0#1 for yes and 0 for no

for bkg in bkgArr:
    a=[]
    for j in range(1000):
        num=den=0
        if(perfectShape == 0):
            #Make real object by sampling. 
            muRandX = np.random.normal(0, delmu)
            muRandY = np.random.normal(0, delmu)
            muArr= [sizex/2.0-0.5+ muRandX, sizey/2.0-0.5+muRandY]
            cov = [[(sigx**2)+np.random.normal(0,delSig), sigxy], [sigxy, (sigy**2)+np.random.normal(0,delSig)]]
            const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
            x, y = np.random.multivariate_normal(muArr, cov, const).T
            x = np.int32(np.round(x))
            y = np.int32(np.round(y))
            obj = np.zeros((sizex,sizey))
            np.add.at(obj, (y,x), 1)
        
# =============================================================================
#         obj = np.zeros((sizex,sizey))
#         x = np.random.normal(sizex/2.0-0.5, sigx, const )
#         y = np.random.normal(sizex/2.0-0.5, sigy, const )
#         x = np.int32(np.round(x))
#         y = np.int32(np.round(y))
#         np.add.at(obj, (y,x), 1)
# 
# =============================================================================
        elif(perfectShape == 1):
            new_sigxSq = (sigx**2)+np.random.normal(0,delSig)
            new_sigySq = (sigy**2)+np.random.normal(0,delSig)
            x = np.linspace(0, sizex-1, sizex)
            y = np.linspace(0, sizey-1, sizey)
            x = x-sizex/2.0+0.5
            y = y-sizey/2.0+0.5
            x, y = np.meshgrid(x, y)
            
            alphax1 = sigx
            alphay1 = sigy
            alphaxy1 = sigxy
            A1 =1/(2*np.pi*alphax1*alphay1*np.sqrt(1- (alphaxy1/(alphax1*alphay1))**2 ))
            arb_const1 = 2*(1 - (alphaxy1/ (alphax1*alphay1))**2 )
            k1=(A1 * np.exp(-((x)**2/(arb_const1*alphax1**2)+ (y)**2/(arb_const1*alphay1**2)  - 2*alphaxy1*(y)*(x)/(arb_const1*alphax1**2 * alphay1**2 ) )))
            
            obj = k1/np.sum(k1)
            const = int(round(np.random.normal(flux, np.sqrt(flux))))
            obj = obj*const
        
        
        
        #Use integer background values
        #noise = np.random.normal(bkg, np.sqrt(bkg), (sizex,sizey))
        newBkg=np.random.normal(bkg,delBkg)
        noise = np.random.normal(newBkg, np.sqrt(newBkg), (sizex,sizey))
        noise = np.round(noise)
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
                alphax=sigx
                alphay=sigy
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
            #mux_calc = 0.0
            #muy_calc = 0.0
            #flux_calc = flux
            
            sigxy_calc = (t1/t2) - (t3*t4)/(t2*t2)
            sigxx_calc = (t5/t2) - (t3*t3) / (t2*t2)
            sigyy_calc = (t6/t2) - (t4*t4) / (t2*t2)
            
            #sigxy_calc = 0.0
            #sigxx_calc = 0.5*(sigx**2)
            #sigyy_calc = 0.5*(sigy**2)
            
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
            
        
        #print(counter,sigxx_calc ,sigyy_calc )
        sigxx_th = 0.5*(sigx**2)
        sigyy_th = 0.5*(sigy**2)
        sigxy_th = 0.5*sigxy
        e1_th = (sigxx_th - sigyy_th)/(sigxx_th + sigyy_th)
        e2_th = 2*sigxy_th/(sigxx_th + sigyy_th)
        ellip_th = np.sqrt(e1_th*e1_th + e2_th*e2_th)
        if(abs(mux_calc)>8 or abs(muy_calc)>8 ):
            
            continue
        if((abs(flux_calc-flux)/np.sqrt(flux)) > 10):
            continue
        if((sigxx_calc + sigyy_calc) <0 or counter>= 97):
            continue
        if(np.isnan(np.sqrt(sigxx_calc + sigyy_calc))):
            continue
        #del x,y,k,obj,noise, tot
        #a.append(flux_calc)
        #mux_calc_arr.append(mux_calc)
        #a.append(np.sqrt( (mux_calc)**2 +(muy_calc)**2  ))
        #a.append(np.sqrt(sigxx_th+sigyy_th) - np.sqrt(sigxx_calc+sigyy_calc))
        #a.append(sigxx_th - sigxx_calc)
        #a.append( np.sqrt(sigxx_calc + sigyy_calc))
        #a.append(mux-mux_calc)
        #a.append((e2-e2_th)  )
        #d11.append(sigxx_th - sigxx_calc)
        #d22.append(sigyy_th - sigyy_calc)
        
        
    angle = 0.5*np.arctan2(e2_th, e1_th)
    area = 2*np.pi*np.sqrt(sigxx_th * sigyy_th - (sigxy/2)**2)
    print (np.nanmean(a), np.nanstd(a), len(a), area)
    
    #print (np.nanstd(a)**2 * (sigxx_th + sigyy_th)**4 *flux / (sigxx_th**2 * sigyy_th**2* np.sqrt(1-abs(e2_th))  )) #Intrinsic Error for e1 in all circumstance
    #print (np.nanstd(a)**2 * (sigxx_th + sigyy_th)**4 *flux*np.sqrt(1-abs(e1_th)) / (sigxx_th**2 * sigyy_th**2 * (1-abs(e2_th))  ))  # = 8 Intrinsic Error for e2 in all circumstance
    
    #print (np.nanstd(a)**2 * (sigxx_th + sigyy_th)**4 *flux**2 / (sigxx_th**2 * sigyy_th**2 * area*bkg*4*np.sqrt(1-abs(e2_th))))  #4AB Error for e1 in all circumstance
    #print (np.nanstd(a)**2 * (sigxx_th + sigyy_th)**4 *flux**2*np.sqrt(1-abs(e1_th)) / (sigxx_th**2 * sigyy_th**2 * area*bkg*4* (1-abs(e2_th)))) # = 164Ab Error for e2 in all circumstance
    #temp = np.nanstd(a) - 1.16*sigx*sigx/np.sqrt(flux)
    
    #print (np.nanstd(a)**2 *flux**2* (1-ellip_th**2)/ (4*area*area* bkg/3.1416), ellip_th) #= 0.25 #4AB component for centroid and size
    print (np.nanstd(a) **2 *flux* (1-ellip_th**2)*np.pi**2.5/ (area), ellip_th) #Intrinsic component for centroid and size
    
    #print ( np.nanstd(a) **2 * (flux /sigxx_th**2 ) ) #=1.3 for sigxx intrinsic comp
    #print ( np.nanstd(a) **2 * (flux**2 /(sigxx_th**2* 4*area*  bkg ) )) #= 2.2 for sigxx bkg comp
    
    #print ( np.nanstd(a) **2 * (flux /sigyy_th**2 ) ) #=1.3 for sigyy intrinsic comp
    #print ( np.nanstd(a) **2 * (flux**2 /(sigyy_th**2* 4*area*  bkg) )) #= 2.2 for sigyy bkg comp
    
    #print ( np.nanstd(a) **2 * (flux /(sigyy_th*sigxx_th*0.5 *(1+abs(e2_th)**2) ) ) ) #=1.3 for sigxy intrinsic comp
    #print ( np.nanstd(a) **2 * (flux**2 /(sigyy_th*sigxx_th*0.5* 4*area*  bkg*(1+abs(e2_th)**2) ) )) #= 2.2 for sigxy bkg comp
    
    #print ( np.nanstd(a) **2 * (sigxx_th + sigyy_th)*flux / ((sigxx_th**2 + sigyy_th**2 ) ) )
    #print (np.sqrt(flux))
    
    #For flux
    #print (np.nanstd(a) **2/ flux)
    #print (np.nanstd(a) **2/ (4*area*bkg))
    #print( np.sqrt( (area/ (np.pi * flux) + 4*area**2*bkg/ (flux**2 * np.pi))   ) )
    errSize  =  np.sqrt( (area/3.1416)/(4*flux) )
    size = np.sqrt(2*sigxx_th )
    print (np.nanstd(a) /(size*errSize))
    #print ( np.nanstd(a) **2 * (flux**2 /(sigyy_th**2* 4*area*  bkg) )) #= 2.2 for sigyy bkg comp
#print ('Final aa') 
#print (np.std(aList), np.mean(aList), ellip_th, area)   

