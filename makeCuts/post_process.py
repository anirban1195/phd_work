#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 20 20:28:00 2021

@author: dutta26
"""


import pandas as pd
import numpy as np 
import sys

#bandList = ['u', 'g', 'r', 'i', 'z']
bandList = ['r']
#u_frame = np.load('/home/dutta26/codes/singleFrame_u.npy')
#g_frame = np.load('/home/dutta26/codes/singleFrame_g.npy')
#r_frame = np.load('/home/dutta26/codes/singleFrame_r.npy')
#i_frame = np.load('/home/dutta26/codes/singleFrame_i.npy')
#z_frame = np.load('/home/dutta26/codes/singleFrame_z.npy')

#u_coadd = pd.read_pickle('/home/dutta26/codes/coaddSc_u_.pk1')
#g_coadd = pd.read_pickle('/home/dutta26/codes/coaddSc_g_.pk1')
#r_coadd = pd.read_pickle('/home/dutta26/codes/coaddSc_r_.pk1')
#i_coadd = pd.read_pickle('/home/dutta26/codes/coaddSc_i_.pk1')
#z_coadd = pd.read_pickle('/home/dutta26/codes/coaddSc_z_.pk1')

#ir_coadd = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')

#First process the single frames 
for band in bandList:
    frame_data = np.load('/home/dutta26/codes/singleFrame_'+band+'.npy')
    coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_'+band+'_.pk1')
    ir_coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')
    
    ir_coadd_data = ir_coadd_data.set_axis(['ra', 'dec', 'star_flag', 'flux', 'mux', 'muy', 'e1', 'e2', 'bkg',
       'size', 'xx', 'yy', 'xy', 'x', 'y', 'force_flag', 'vert_flag',
       'bad_flag', 'back_sky', 'seeing', 'zp', 'fwhm', 'mjd', 'airmass',
       'mphase', 'mAngle', 'expTime', 'focus', 'zp_n','skymag', 'depth', 'mRa',
       'mDec', 'mag', 'corrxx', 'corryy', 'corrxy', 'corr_flag', 'Empty', 'Empty',
       'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty',
       'Empty', 'Empty',], axis=1, inplace=False)
    
    ir_coadd_data.loc[:,'corrxx'] = ir_coadd_data['xx']
    ir_coadd_data.loc[:,'corryy'] = ir_coadd_data['yy']
    ir_coadd_data.loc[:,'corrxy'] = ir_coadd_data['xy']
    
    coadd_data = coadd_data.set_axis(['ra', 'dec', 'star_flag', 'flux', 'mux', 'muy', 'e1', 'e2', 'bkg',
       'size', 'xx', 'yy', 'xy', 'x', 'y', 'force_flag', 'vert_flag',
       'bad_flag', 'back_sky', 'seeing', 'zp', 'fwhm', 'mjd', 'airmass',
       'mphase', 'mAngle', 'expTime', 'focus', 'zp_n','skymag', 'depth', 'mRa',
       'mDec', 'mag', 'corrxx', 'corryy', 'corrxy', 'corr_flag', 'Empty', 'Empty',
       'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty',
       'Empty', 'Empty',], axis=1, inplace=False)
    coadd_data.loc[:,'corrxx'] = coadd_data['xx']
    coadd_data.loc[:,'corryy'] = coadd_data['yy']
    coadd_data.loc[:,'corrxy'] = coadd_data['xy']
    
    
    a,b,c = np.shape(frame_data)
    store = np.zeros((a+2, b ,c), dtype = np.float32)
    
    #Process frame data for magnitudes
    for j in range(a):
        temp = frame_data[j,:,:]
        indexList = np.where(temp[:,3] != 0)[0]
        frame_data[j, indexList,33 ] = frame_data[j, indexList, 20 ] - 2.5*np.log10(frame_data[j, indexList, 3 ]/frame_data[j, indexList, 26 ])
        del temp, indexList
    frame_data[:, :,34 ] = frame_data[:, :,10 ]
    frame_data[:, :,35 ] = frame_data[:, :,11 ]
    frame_data[:, :,36 ] = frame_data[:, :,12 ]
    
    #First process IR coadd 
    star_arr = ir_coadd_data.loc[(ir_coadd_data['star_flag'] == 1) & (ir_coadd_data['vert_flag'] == 0) & (ir_coadd_data['force_flag'] == 0)]
    star_temp = np.zeros(( len(star_arr['xx']) , 6)   , dtype = np.float32)
    star_temp[:,0] = star_arr['xx']
    star_temp[:,1] = star_arr['yy']
    star_temp[:,2] = star_arr['xy']
    star_temp[:,3] = star_arr['x']
    star_temp[:,4] = star_arr['y']
    
    for j in range(b):
        #print (j)
        if(ir_coadd_data['flux'][j] == 0 or ir_coadd_data['e1'][j] == 0 or np.isnan(ir_coadd_data['e1'][j])):
            continue
        #For each object find the nearest 10 neighbors
        x = ir_coadd_data['x'][j]
        y = ir_coadd_data['y'][j]
        temp = np.copy(star_temp)
        temp[:,5] = np.sqrt((temp[:,3]-x)**2 + (temp[:,4]-y)**2)
        temp = temp[temp[:,5].argsort()]
        #If same star then continue 
        
        if(temp[0, 5]< 5):
            temp = np.delete(temp, 0, 0)
        if(np.isnan(np.sum( temp[0:10, 0] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
            sys.exit()
        if(np.isnan(np.sum( temp[0:10, 1] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
            sys.exit()
        if(np.isnan(np.sum( temp[0:10, 2] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
            sys.exit()
        avgSigxx = np.sum( temp[0:10, 0] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
        avgSigyy = np.sum( temp[0:10, 1] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
        avgSigxy = np.sum( temp[0:10, 2] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
        
        if((ir_coadd_data.loc[j].at['xx'] - avgSigxx + 1)< 0 or (ir_coadd_data.loc[j].at['yy'] - avgSigyy + 1)< 0 ):
            continue
        
        ir_coadd_data.loc[j].at['corrxx'] = ir_coadd_data.loc[j].at['xx'] - avgSigxx + 1
        ir_coadd_data.loc[j].at['corryy'] = ir_coadd_data.loc[j].at['yy'] - avgSigyy + 1
        ir_coadd_data.loc[j].at['corrxy'] = ir_coadd_data.loc[j].at['xy'] - avgSigxy 
        ir_coadd_data.loc[j].at['corr_flag'] = 1
        del temp
    del star_arr, star_temp    
    #Now process the band coadd 
    star_arr = coadd_data.loc[(coadd_data['star_flag'] == 1) & (coadd_data['vert_flag'] == 0) & (coadd_data['force_flag'] == 0)]
    star_temp = np.zeros(( len(star_arr['xx']) , 6)   , dtype = np.float32)
    star_temp[:,0] = star_arr['xx']
    star_temp[:,1] = star_arr['yy']
    star_temp[:,2] = star_arr['xy']
    star_temp[:,3] = star_arr['x']
    star_temp[:,4] = star_arr['y']
    
    for j in range(b):
        if(coadd_data['flux'][j] == 0 or coadd_data['e1'][j] == 0 or np.isnan(coadd_data['e1'][j])):
            continue
        #print (j)
        #For each object find the nearest 10 neighbors
        x = coadd_data['x'][j]
        y = coadd_data['y'][j]
        temp = np.copy(star_temp)
        temp[:,5] = np.sqrt((temp[:,3]-x)**2 + (temp[:,4]-y)**2)
        temp = temp[temp[:,5].argsort()]
        #If same star then continue 
        if(temp[0, 5]< 5):
            temp = np.delete(temp, 0, 0)
        if(np.isnan(np.sum( temp[0:10, 0] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
            sys.exit()
        if(np.isnan(np.sum( temp[0:10, 1] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
            sys.exit()
        if(np.isnan(np.sum( temp[0:10, 2] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
            sys.exit()
        avgSigxx = np.sum( temp[0:10, 0] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
        avgSigyy = np.sum( temp[0:10, 1] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
        avgSigxy = np.sum( temp[0:10, 2] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
        
        if((coadd_data.loc[j].at['xx'] - avgSigxx+ 1)< 0 or (coadd_data.loc[j].at['yy'] - avgSigyy+ 1)< 0 ):
            continue
        
        coadd_data.loc[j].at['corrxx'] = coadd_data.loc[j].at['xx'] - avgSigxx+ 1
        coadd_data.loc[j].at['corryy'] = coadd_data.loc[j].at['yy'] - avgSigyy+ 1
        coadd_data.loc[j].at['corrxy'] = coadd_data.loc[j].at['xy'] - avgSigxy
        coadd_data.loc[j].at['corr_flag'] = 1
        del temp
    del star_arr, star_temp
        
    print ('vvvvvvvvvvv')    
    #Now process indiviual frames
    for k in range(a):
        print (k)
        #Star = 1 , vert flag =0 , bad_flag = 0 and force = 0
        star_arr = frame_data[k, (np.where((frame_data[k,:,2] == 1) & (frame_data[k,:,16] == 0) & (frame_data[k,:,17] == 0) & (frame_data[k,:,15] == 0)))[0],: ]
        q,r = np.shape(star_arr)
        star_temp = np.zeros(( q , 6)   , dtype = np.float32)
        star_temp[:,0] = star_arr[:, 10]
        star_temp[:,1] = star_arr[:, 11]
        star_temp[:,2] = star_arr[:, 12]
        star_temp[:,3] = star_arr[:, 13]
        star_temp[:,4] = star_arr[:, 14]
        
        for j in range(b):
            
            if(frame_data[k,j,3] == 0 or frame_data[k,j,6] == 0 or np.isnan(frame_data[k,j,6])):
                continue
            #For each object find the nearest 10 neighbors
            x = frame_data[k,j, 13]
            y = frame_data[k,j, 14]
            temp = np.copy(star_temp)
            temp[:,5] = np.sqrt((temp[:,3]-x)**2 + (temp[:,4]-y)**2)
            temp = temp[temp[:,5].argsort()]
            
            #If same star then continue 
            if(temp[0, 5]< 5):
                temp = np.delete(temp, 0, 0)
            
            if(np.isnan(np.sum( temp[0:10, 0] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
                sys.exit()
            if(np.isnan(np.sum( temp[0:10, 1] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
                sys.exit()
            if(np.isnan(np.sum( temp[0:10, 2] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
                sys.exit()
               
            avgSigxx = np.sum( temp[0:10, 0] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
            avgSigyy = np.sum( temp[0:10, 1] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
            avgSigxy = np.sum( temp[0:10, 2] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
            
            if((frame_data[k,j, 10] - avgSigxx+ 1)< 0 or (frame_data[k,j, 11] - avgSigyy+ 1)< 0 ):
                continue
            
            frame_data[k,j, 34] = frame_data[k,j, 10] - avgSigxx+ 1
            frame_data[k,j, 35] = frame_data[k,j, 11] - avgSigyy+ 1
            frame_data[k,j, 36] = frame_data[k,j, 12] - avgSigxy
            frame_data[k,j, 37] = 1
            del temp
        del star_arr, star_temp
        
    
    #Now rewrite all data to a single frame
    store[0,:,:] = ir_coadd_data.to_numpy(dtype = np.float32)
    store[1,:,:] = coadd_data.to_numpy(dtype = np.float32)
    store[2:a+2, :, :] = frame_data
    
    np.save('/scratch/halstead/d/dutta26/abell_2390/processed_' +str(band)+'.npy', store)
    
    
    
    
    
    
    
    
    
    
    
    
    