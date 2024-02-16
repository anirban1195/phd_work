#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 15:46:49 2023

@author: dutta26
"""
import numpy as np

def integrate(z):
    arr= np.arange(0, z, 0.01)
    #a= np.sum( 0.01/np.sqrt(0.3*(1+arr)**3 + 0.7))
    #return a
    return np.sum( 0.01/np.sqrt(0.3*(1+arr)**3 + 0.7))*(1/(1+z))

def getF(x):
    if(x<1):
        #b = (1/np.sqrt( 1 - x**2 ))* np.log(x+ np.sqrt(x**2 +1))
        b = (1/np.sqrt( 1 - x**2 ))*  np.arccosh(1/x)
        a = (1/(x**2 - 1)) * (1 - b)
    elif(x == 1):
        a= 0.3333
    elif(x>1):
        b = (1/np.sqrt(x**2 - 1))* np.arccos(1/x)
        a = (1/(x**2 - 1)) * (1 - b)
    else:
        a= None
    return a
        
       
    
def getg(x):
    if(x<1):
       #b = (1/np.sqrt( 1 - x**2 ))* np.log(x+ np.sqrt(x**2 +1))
       b = (1/np.sqrt( 1 - x**2 ))* np.arccosh(1/x)
       a = np.log(x/2) + b
       
    elif(x == 1):
        a = 1+np.log(1/2)
    elif(x>1):
        b = (1/np.sqrt(x**2 - 1))* np.arccos(1/x)
        a = np.log(x/2) + b
    else:
        a=None
    return a

def getGamma(centra, centdec, zLens, ra, dec, zs, rs):
    deldec = dec - centdec
    delra = ra - centra 
    dl = 4285.71* integrate(zLens)
    ds = 4285.71* integrate(zs)
    dls = (ds - dl)
    #print (dl, ds, dls)
    ks = 9.65e-4 * (dl*dls)/ds
    thetas = (rs/dl)*(180/np.pi) #In degrees
    R = np.sqrt(deldec**2 + delra**2) #In degrees
    x = R/thetas
    print (x)
    #Limit to weak lensing regime
    if(x<1.25):
        x=1.25
    
    gammaTot = 2*ks* ( 2*getg(x)/x**2 - getF(x))
    kappaTot = 2*ks*getF(x)
    
    angle = np.arctan2(deldec, delra)
    gamma1 = np.cos(2*angle)*gammaTot
    gamma2 = np.sin(2*angle)*gammaTot
    return (-1*gamma1, -1*gamma2, kappaTot) 
    

# =============================================================================
# x= np.arange(0.5, 2, 0.03)
# y=[]
# for xVal in x:
#     y.append(getF(xVal))
# import matplotlib.pyplot as plt    
# plt.plot(x,y, 'b.')
# =============================================================================
