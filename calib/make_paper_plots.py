#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 19:05:35 2023

@author: dutta26
"""

import numpy as np
import matplotlib.pyplot as plt
data = np.load('/scratch/bell/dutta26/abell_2390/paper_plot_data_sf_vs_combined_15pert.npy')

snr = np.log10(data[:,7]/np.sqrt(4*data[:,8]*data[:,11]))

fig, axs = plt.subplots(2, 2, sharex = True )
axs[0, 0].plot(snr, data[:,6], 'b.', markersize = 1)
axs[0, 0].set( ylabel = '% Failure')
axs[0, 1].plot(snr, data[:,0]*100, 'b.', markersize = 1)
#axs[0, 1].plot(snr, data[:,3]*100, 'r.', markersize = 1)
axs[0, 1].set( ylabel = '% Error in Flux')
axs[1, 0].plot(snr, data[:,1]*100, 'b.', markersize = 1)
#axs[1, 0].plot(snr, data[:,4]*100, 'r.', markersize = 1)
axs[1, 0].set( xlabel = r'Log( SNR )', ylabel = r' % Error in $\sigma_{xx}$')
axs[1, 1].plot(snr, data[:,2]*100, 'b.', markersize = 1)
#axs[1, 1].plot(snr, data[:,5]*100, 'r.', markersize = 1)
axs[1, 1].set( xlabel = r'Log( SNR )', ylabel = r'% Error in $\sigma_{xy}$')

# =============================================================================
# plt.figure(1)
# plt.subplot(221)
# plt.plot(snr, data[:,6], 'b.')
# 
# plt.subplot(222)
# plt.plot(snr, data[:,0]*100, 'b.', markersize = 1)
# plt.plot(snr, data[:,3]*100, 'g.', markersize = 1)
# 
# 
# 
# plt.subplot(223)
# plt.plot(snr, data[:,1]*100, 'b.', markersize = 1)
# plt.plot(snr, data[:,4]*100, 'g.', markersize = 1)
# 
# plt.subplot(224)
# plt.plot(snr, data[:,2]*100, 'b.', markersize = 1)
# plt.plot(snr, data[:,5]*100, 'g.', markersize = 1)
# # 
# # plt.subplot(222)
# # plt.errorbar(np.log10(snrArr), measuredArr_sigxy, yerr= measuredErrArr_sigxy, fmt = '.k')
# # #plt.errorbar(np.log10(snrArr), correctedArr_sigxy, yerr= measuredErrArr_sigxy, fmt = '.r')
# # #plt.errorbar(np.log10(snrArr), sigxy* np.ones(len(snrArr))/2, yerr= measuredErrArr_sigxy, fmt = '.b')
# # plt.errorbar(np.log10(snrArr), e2* np.ones(len(snrArr)), yerr= measuredErrArr_sigxy, fmt = '.b')
# # 
# # plt.ylabel('Sigmaxy')
# # 
# # plt.subplot(223)
# # plt.errorbar(np.log10(snrArr), np.log10(measuredArr_flux), yerr= np.log10(measuredErrArr_flux), fmt = '.k')
# # #plt.errorbar(np.log10(snrArr), np.log10(correctedArr_flux), yerr= np.log10(measuredErrArr_flux), fmt = '.r')
# # plt.errorbar(np.log10(snrArr), np.log10(fluxArr), yerr= np.log10(measuredErrArr_flux), fmt = '.b')
# # plt.ylabel('Log Flux')
# # 
# =============================================================================
