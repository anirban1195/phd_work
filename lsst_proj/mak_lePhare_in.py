#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 18:25:40 2022

@author: dutta26
"""
import pandas as pd
import numpy as np

zpArr=[29.21, 31.23, 31.37, 31.18, 30.84, 30.35 ]
f=open('/scratch/halstead/d/dutta26/lsst/magList.in', 'w+')

df0 = pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df1_0.pk1')
df0=  np.array(df0)

df1 = pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df1_1.pk1')
df1 = np.array(df1)

df2 = pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df1_2.pk1')
df2 = np.array(df2)

df3 = pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df1_3.pk1')
df3 = np.array(df3)

df4 = pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df1_4.pk1')
df4 = np.array(df4)

df5 = pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df1_5.pk1')
df5 = np.array(df5)

for j in range(len(df1)):
    mag0= err0=mag1=err1=mag2=err2=mag3=err3=mag4=err4=mag5=err5=-99
    
    if(df0[j,3] == None or np.isnan(df0[j,3]) or df0[j,3]<=0.5 or np.isinf(df0[j,3])):
        pass
    else:
        mag0 = zpArr[0] - 2.5*np.log10(df0[j,3])
        err0 = np.sqrt(df0[j,3])/df0[j,3]
    
    if(df1[j,3] == None or np.isnan(df1[j,3]) or df1[j,3]<=0.5 or np.isinf(df1[j,3])):
        pass
    else:
        mag1 = zpArr[1] - 2.5*np.log10(df1[j,3])
        err1 = np.sqrt(df1[j,3])/df1[j,3]
        
    if(df2[j,3] == None or np.isnan(df2[j,3]) or df2[j,3]<=0.5 or np.isinf(df2[j,3])):
        pass
    else:
        mag2 = zpArr[2] - 2.5*np.log10(df2[j,3])
        err2 = np.sqrt(df2[j,3])/df2[j,3]
        
    if(df3[j,3] == None or np.isnan(df3[j,3]) or df3[j,3]<=0.5 or np.isinf(df3[j,3])):
        pass
    else:
        mag3 = zpArr[3] - 2.5*np.log10(df3[j,3])
        err3 = np.sqrt(df3[j,3])/df3[j,3]
        
    if(df4[j,3] == None or np.isnan(df4[j,3]) or df4[j,3]<=0.5 or np.isinf(df4[j,3])):
        pass
    else:
        mag4 = zpArr[4] - 2.5*np.log10(df4[j,3])
        err4 = np.sqrt(df4[j,3])/df4[j,3]
        
    if(df5[j,3] == None or np.isnan(df5[j,3]) or df5[j,3]<=0.5 or np.isinf(df5[j,3])):
        pass
    else:
        mag5 = zpArr[5] - 2.5*np.log10(df5[j,3])
        err5 = np.sqrt(df5[j,3])/df5[j,3]
        
    f.write(str(j)+' '+str(mag0)[0:5]+ ' '+str(err0)[0:5]+' '+str(mag1)[0:5]+ ' '+str(err1)[0:5]+' '+str(mag2)[0:5]+ ' '+str(err2)[0:5] +' '+str(mag3)[0:5]+ ' '+str(err3)[0:5]+
            ' '+str(mag4)[0:5]+ ' '+str(err4)[0:5] + ' '+str(mag5)[0:5]+ ' '+str(err5)[0:5]+ '\n')
    
f.close()
        
    
        
    
    
    