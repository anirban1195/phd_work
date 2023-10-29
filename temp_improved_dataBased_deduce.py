#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 12:12:23 2023

@author: dutta26
"""
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
store = np.vstack((np.load('/scratch/bell/dutta26/abell_2390/measure_test_data2m.npy'),np.load('/scratch/bell/dutta26/abell_2390/measure_test_data2m_supp.npy')))
elong = store[:,2]/store[:,3]
area = np.pi*store[:,2]*store[:,3]
e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
e2 = 2* store[:,7]/(store[:,5] + store[:,6])
ellip = np.sqrt(e1**2 + e2**2)

# =============================================================================
# #FLUX ##########################
# loc = np.where((store[:,0]>0) & (store[:,8]==1)& (elong>0.3) & (elong<20) )[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# ellip = np.sqrt(e1**2 + e2**2)
# def ef_flux_intrinsic(X , a):
#     N,A,B,e1,e2, elong = X
#     temp = 1-elong
#     ellip = np.sqrt(e1**2 + e2**2)
#     return np.sqrt(a*N)
# popt, pcov = curve_fit(ef_flux_intrinsic, (store[:,0],  area, store[:,1], e1, e2, 1/elong), store[:,9])
# print (popt, pcov)
# #plt.plot( np.arange(len(store[:,9]))  , store[:,9], 'b.')
# #plt.plot(  np.arange(len(store[:,9])) , ef_flux_intrinsic((store[:,0],  area, store[:,1], e1, e2, 1/elong), popt[0]), 'r+')
# =============================================================================

# =============================================================================
# loc = np.where((store[:,0]>0) & (store[:,8]==2)& (elong<10.8))[0]
# store = store[loc,:]
# area = np.pi*store[:,2]*store[:,3]
# 
# def ef_flux_bkg(X, a):
#     N,A,B = X
#     return np.sqrt(1.0*N + a*A*B)
# popt, pcov = curve_fit(ef_flux_bkg, (store[:,0],  area, store[:,1]), store[:,9],[4])
# print (popt, pcov)
# plt.plot( np.arange(len(store[:,9]))  , np.log10(store[:,9]), 'b.')
# plt.plot(  np.arange(len(store[:,9])) , np.log10(ef_flux_bkg((store[:,0],  area, store[:,1]), popt[0])), 'r+')
# 
# =============================================================================



# =============================================================================
# #SIZE ###############
# loc = np.where((store[:,0]>0) & (store[:,8]==1)& (ellip == 0.6) )[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# def ef_size_intrinsic(X, a):
#     N,A,B,e1,e2, elong = X
#     temp = 1-elong
#     ellip = np.sqrt(e1**2 + e2**2)
#     return np.sqrt(a*A/(np.pi*N* (1-ellip**2)))
# popt, pcov = curve_fit(ef_size_intrinsic, (store[:,0],  area, store[:,1], e1, e2, 1/elong), store[:,10],[1])
# print (popt, pcov)
# #plt.plot( np.arange(len(store[:,10]))  , np.log10(store[:,10]), 'b.')
# #plt.plot( np.arange(len(store[:,10]))  , np.log10(ef_size_intrinsic((store[:,0],  area, store[:,1], e1, e2, 1/elong), 0.18)), 'r+')
# 
# 
# =============================================================================



# =============================================================================
# loc = np.where((store[:,0]>0) & (store[:,8]==2)& (ellip ==0.0)  & (store[:,1] == 50))[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# def ef_size_bkg(X, a):
#     N,A,B,e1,e2, elong = X
#     ellip = np.sqrt(e1**2 + e2**2)
#     temp = 1/elong**2
#     return np.sqrt(a*A**2 *B/(np.pi* N**2 * (1-ellip**2)))##################################
# popt, pcov = curve_fit(ef_size_bkg, (store[:,0],  area, store[:,1], e1, e2, 1/elong), store[:,10],[1])
# print (popt, pcov)
# #plt.plot( np.arange(len(store[:,10]))  , np.log10(store[:,10]), 'b.')
# #plt.plot( np.arange(len(store[:,10]))  , np.log10(ef_size_bkg((store[:,0],  area, store[:,1], e1, e2, 1/elong), popt[0])), 'r+')
# 
# =============================================================================




#CENTROID#################################

# =============================================================================
# loc = np.where((store[:,0]>0) & (store[:,8]==1)& (ellip==0.6))[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# def ef_cent_intrinsic(X, a):
#     N,A,B,e1,e2, elong = X
#     temp = 1-elong
#     ellip = np.sqrt(e1**2 + e2**2)
#     return np.sqrt(a*A/(np.pi*N* (1-ellip**2)))
# popt, pcov = curve_fit(ef_cent_intrinsic, (store[:,0],  area, store[:,1], e1, e2, 1/elong), store[:,11],[1])
# print (popt, pcov)
# #plt.plot( np.arange(len(store[:,11]))  , np.log10(store[:,11]), 'b.')
# #plt.plot( np.arange(len(store[:,11]))  , np.log10(ef_cent_intrinsic((store[:,0],  area, store[:,1], e1, e2, 1/elong), 0.18)), 'r+')
# 
# 
# =============================================================================
# =============================================================================
# loc = np.where((store[:,0]>0) & (store[:,8]==2)& (ellip==0.0) & (store[:,1] == 400))[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# def ef_cent_bkg(X, a):
#     N,A,B,e1,e2, elong = X
#     ellip = np.sqrt(e1**2 + e2**2)
#     temp = 1/elong**2
#     return np.sqrt(a*A**2 *B/(np.pi* N**2 * (1-ellip**2)))##################################
# popt, pcov = curve_fit(ef_cent_bkg, (store[:,0],  area, store[:,1], e1, e2, 1/elong), store[:,11],[1])
# print (popt, pcov)
# #plt.plot( np.arange(len(store[:,11]))  , np.log10(store[:,11]), 'b.')
# #plt.plot( np.arange(len(store[:,11]))  , np.log10(ef_cent_bkg((store[:,0],  area, store[:,1], e1, e2, 1/elong), popt[0])), 'r+')
# 
# 
# =============================================================================






#SIGMAXX####################

# =============================================================================
# loc = np.where((store[:,0]>0) & (store[:,8]==1)& (elong<10.0985)  )[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# def ef_xx_intrinsic(X, lam, mu):
#     N,A,B,e1,e2, xx , yy, xy, elong = X
#     ellip = np.sqrt(e1**2 + e2**2)
#     size = np.sqrt((xx + yy)/2)
#     temp = (1-elong)**0.5
#     #print (temp)
#     return np.sqrt(lam * (size)**4 * ((1+e1)**2 +mu ) / (N) )
#     #return np.sqrt( lam * (A)**2   / (N*(1-temp)**2) )
#     #return np.sqrt( lam * (xx)**2  / (N) )
# popt, pcov = curve_fit(ef_xx_intrinsic, (store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong), store[:,12],[1, 1])
# print (popt, pcov)
# plt.plot( np.arange(len(store[:,12]))  , np.log10(store[:,12]), 'b.')
# plt.plot( np.arange(len(store[:,12]))  , np.log10(ef_xx_intrinsic((store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong), popt[0], popt[1])), 'r+')
# 
# =============================================================================


# =============================================================================
# loc = np.where((store[:,0]>0) & (store[:,8]==2)& (elong<10.505)  )[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# def ef_xx_bkg(X, lam,mu):
#     N,A,B,e1,e2, xx , yy, xy, elong = X
#     ellip = np.sqrt(e1**2 + e2**2)
#     size = np.sqrt((xx + yy)/2)
#     temp = 1/elong**2
#     #print(temp)
#     return np.sqrt(lam * (size)**4 * ((1+e1)**2 +mu )*A*B / (N**2) )
#     #return np.sqrt(lam * (A)**2 *A*B*temp / (N**2) ) #WOW 
# popt, pcov = curve_fit(ef_xx_bkg, (store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong ), store[:,12],[1,0])
# print (popt, pcov)
# plt.plot( np.arange(len(store[:,12]))  , np.log10(store[:,12]), 'b.')
# plt.plot( np.arange(len(store[:,12]))  , np.log10(ef_xx_bkg((store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong), popt[0], popt[1])), 'r+')
# 
# 
# =============================================================================


#SIGMAYY############################

# =============================================================================
# loc = np.where((store[:,0]>0) & (store[:,8]==1)& (elong<10.0985)  )[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# def ef_yy_intrinsic(X, lam, mu):
#     N,A,B,e1,e2, xx , yy, xy, elong = X
#     ellip = np.sqrt(e1**2 + e2**2)
#     size = np.sqrt((xx + yy)/2)
#     temp = (1-elong)**0.5
#     #print (temp)
#     return np.sqrt(lam * (size)**4 * ((1-e1)**2 +mu ) / (N) )
#     #return np.sqrt( lam * (A)**2   / (N*(1-temp)**2) )
#     #return np.sqrt( lam * (xx/2)**2  / (N) )
# popt, pcov = curve_fit(ef_yy_intrinsic, (store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong), store[:,13],[1, 1])
# print (popt, pcov)
# plt.plot( np.arange(len(store[:,13]))  , np.log10(store[:,13]), 'b.')
# plt.plot( np.arange(len(store[:,13]))  , np.log10(ef_yy_intrinsic((store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong), popt[0], popt[1])), 'r+')
# 
# =============================================================================

# =============================================================================
# loc = np.where((store[:,0]>0) & (store[:,8]==2)& (elong<10.505)  )[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# def ef_yy_bkg(X, lam,mu):
#     N,A,B,e1,e2, xx , yy, xy, elong = X
#     ellip = np.sqrt(e1**2 + e2**2)
#     size = np.sqrt((xx + yy)/2)
#     temp = 1/elong**2
#     #print(temp)
#     return np.sqrt(lam * (size)**4 * ((1-e1)**2 +mu )*A*B / (N**2) )
#     #return np.sqrt(lam * (A)**2 *A*B*temp / (N**2) ) #WOW 
# popt, pcov = curve_fit(ef_yy_bkg, (store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong ), store[:,13],[1,0])
# print (popt, pcov)
# plt.plot( np.arange(len(store[:,13]))  , np.log10(store[:,13]), 'b.')
# plt.plot( np.arange(len(store[:,13]))  , np.log10(ef_yy_bkg((store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong), popt[0], popt[1])), 'r+')
# 
# 
# 
# =============================================================================




#SIGMAXY #################
# =============================================================================
# loc = np.where((store[:,0]>0) & (store[:,8]==1)& (elong<10.0985)  )[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# def ef_xy_intrinsic(X, lam, mu):
#     N,A,B,e1,e2, xx , yy, xy, elong = X
#     ellip = np.sqrt(e1**2 + e2**2)
#     size = np.sqrt((xx + yy)/2)
#     temp = (1-elong)**0.5
#     #print (temp)
#     #return np.sqrt(lam * (size)**4 * ((e2)**2 +mu ) / (N) )
#     return np.sqrt( lam * xx*yy*0.5*(1+e2**2) / N )
#     #return np.sqrt( lam * (xx/2)**2  / (N) )
# popt, pcov = curve_fit(ef_xy_intrinsic, (store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong), store[:,14],[1, 1])
# print (popt, pcov)
# plt.plot( np.arange(len(store[:,14]))  , store[:,14], 'b.')
# plt.plot( np.arange(len(store[:,14]))  , ef_xy_intrinsic((store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong), popt[0], popt[1]), 'r+')
# 
# 
# =============================================================================
# =============================================================================
# loc = np.where((store[:,0]>0) & (store[:,8]==2)& (elong<1.505)  )[0]
# store = store[loc,:]
# elong = elong[loc]
# area = np.pi*store[:,2]*store[:,3]
# e1 = (store[:,5] - store[:,6])/(store[:,5] + store[:,6])
# e2 = 2* store[:,7]/(store[:,5] + store[:,6])
# def ef_xy_bkg(X, lam,mu):
#     N,A,B,e1,e2, xx , yy, xy, elong = X
#     ellip = np.sqrt(e1**2 + e2**2)
#     size = np.sqrt((xx + yy)/2)
#     temp = 1/elong**2
#     #print(temp)
#     #return np.sqrt(lam * (size)**4 * ((e2)**2 +mu )*A*B / (N**2) )
#     return np.sqrt( lam * xx*yy*0.5*(1+e2**2) * A*B / N**2 )
# popt, pcov = curve_fit(ef_xy_bkg, (store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong ), store[:,14],[1,0])
# print (popt, pcov)
# plt.plot( np.arange(len(store[:,14]))  , np.log10(store[:,14]), 'b.')
# plt.plot( np.arange(len(store[:,14]))  , np.log10(ef_xy_bkg((store[:,0],  area, store[:,1], e1, e2, store[:,5], 
#                                                           store[:,6], store[:,7], 1/elong), popt[0], popt[1])), 'r+')
# 
# =============================================================================


# =============================================================================
# colors = area
# ellip = np.sqrt(e1**2+e2**2)
# plt.scatter(store[:,0], store[:,11], c=colors, cmap = 'viridis', marker = '.', s= 0)
# plt.clim(20,160)
# 
# # =============================================================================
# # plt.scatter(np.log10(store[:,0]),  np.log10(ef_cent_bkg((store[:,0],  area, store[:,1], e1, e2, 1/elong), 0.9)), c=colors, cmap = 'viridis', marker = '+', s= (100+ellip))
# # plt.clim(30,70)
# # 
# # plt.scatter(np.log10(1.1*store[:,0]),  np.log10(ef_cent_bkg((store[:,0],  area, store[:,1], e1, e2, 1/elong), 1)), c=colors, cmap = 'viridis', marker = 'v', s= (100+ellip))
# # plt.clim(30,70)
# # 
# # =============================================================================
# 
# cmap =plt.cm.viridis
# norm= Normalize(vmin=20,vmax = 160)
# for j in range(int(len(store)/5)):
#     plt.errorbar(store[j*5: (j+1)*5,0], store[j*5: (j+1)*5,11], yerr= store[j*5: (j+1)*5,20], fmt='k.',  color = matplotlib.colors.rgb2hex(cmap(norm(area[j*5]))), markersize=10)
#     plt.plot(store[j*5: (j+1)*5,0],  ef_cent_bkg((store[j*5: (j+1)*5,0],  area[j*5: (j+1)*5], store[j*5: (j+1)*5,1], e1[j*5: (j+1)*5], e2[j*5: (j+1)*5], 1/elong[j*5: (j+1)*5]), 0.9), '--', color = matplotlib.colors.rgb2hex(cmap(norm(area[j*5]))))
#     plt.plot(store[j*5: (j+1)*5,0],  ef_cent_bkg((store[j*5: (j+1)*5,0],  area[j*5: (j+1)*5], store[j*5: (j+1)*5,1], e1[j*5: (j+1)*5], e2[j*5: (j+1)*5], 1/elong[j*5: (j+1)*5]), 1), '-', color = matplotlib.colors.rgb2hex(cmap(norm(area[j*5]))))
# 
# plt.yscale('log')  
# plt.xscale('log')       
# plt.title (' Background Centroid error when e=0.0 and B=400  \n Coefficients: Dashed=0.9. Solid = 1')
# plt.colorbar(orientation='vertical', label= 'Area')
# plt.xlabel('Flux')
# plt.ylabel('Uncertainty')
# =============================================================================

#For xx, yy and xy
colors = area
ellip = np.sqrt(e1**2+e2**2)
plt.scatter(store[:,0], store[:,11], c=colors, cmap = 'viridis', marker = '.', s= 0)
plt.clim(20,160)

# =============================================================================
# plt.scatter(np.log10(store[:,0]),  np.log10(ef_cent_bkg((store[:,0],  area, store[:,1], e1, e2, 1/elong), 0.9)), c=colors, cmap = 'viridis', marker = '+', s= (100+ellip))
# plt.clim(30,70)
# 
# plt.scatter(np.log10(1.1*store[:,0]),  np.log10(ef_cent_bkg((store[:,0],  area, store[:,1], e1, e2, 1/elong), 1)), c=colors, cmap = 'viridis', marker = 'v', s= (100+ellip))
# plt.clim(30,70)
# 
# =============================================================================

cmap =plt.cm.viridis
norm= Normalize(vmin=20,vmax = 160)
for j in range(int(len(store)/5)):
    plt.errorbar(store[j*5: (j+1)*5,0], store[j*5: (j+1)*5,11], yerr= store[j*5: (j+1)*5,20], fmt='k.',  color = matplotlib.colors.rgb2hex(cmap(norm(area[j*5]))), markersize=10)
    plt.plot(store[j*5: (j+1)*5,0],  ef_cent_bkg((store[j*5: (j+1)*5,0],  area[j*5: (j+1)*5], store[j*5: (j+1)*5,1], e1[j*5: (j+1)*5], e2[j*5: (j+1)*5], 1/elong[j*5: (j+1)*5]), 0.9), '--', color = matplotlib.colors.rgb2hex(cmap(norm(area[j*5]))))

plt.yscale('log')  
plt.xscale('log')       
plt.title (' Background Centroid error when e=0.0 and B=400  \n Coefficients: Dashed=0.9. Solid = 1')
plt.colorbar(orientation='vertical', label= 'Area')
plt.xlabel('Flux')
plt.ylabel('Uncertainty')