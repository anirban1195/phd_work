#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 17:58:46 2024

@author: dutta26
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats

dataSet_i = np.load('/scratch/bell/dutta26/test_data_i.npy')

dataSet_r = np.load('/scratch/bell/dutta26/test_data_r.npy')

dataSet = np.vstack((dataSet_i, dataSet_r))

snrArr = 10**np.linspace(-0.5,1.5, 10)

fig, axs = plt.subplots(2, 2, sharex = True )
axs[0,0].set_xscale('log')
tot = 0
# =============================================================================
# axs[0,1].set_yscale('log')
# axs[1,0].set_yscale('log')
# axs[1,1].set_yscale('log')
# =============================================================================
for j in range(len(snrArr) - 1):
    snr_low = snrArr[j]
    snr_high = snrArr[j+1]
    snr_mid= 0.5*(snr_high+snr_low)
    
    snr_cut_loc = np.where((dataSet[:,:, 0]>snr_low) & (dataSet[:,:, 0]<snr_high))
    master_arr = dataSet[snr_cut_loc[0], snr_cut_loc[1],:]
    
    
    #For convergence
    loc_c_fail = np.where((master_arr[:,1] == -99))
    loc_c_success = np.where((master_arr[:,1] != 0) & (master_arr[:,1] != -99) )
    
    print((len(loc_c_fail[0]) + len(loc_c_success[0])))
    tot += (len(loc_c_fail[0]) + len(loc_c_success[0]))
    
    fail_rate_c = len(loc_c_fail[0])/(len(loc_c_fail[0]) + len(loc_c_success[0]))
    axs[0, 0].plot(snr_mid, fail_rate_c , 'b.', markersize = 8)
    
    temp = (master_arr[loc_c_success[0], 1] - master_arr[loc_c_success[0], 13])/ master_arr[loc_c_success[0], 13]
    mean, frac_flux_err_c_med, frac_flux_err_c_std = sigma_clipped_stats(temp.flatten())
    axs[0, 1].errorbar(snr_mid, frac_flux_err_c_med, yerr = frac_flux_err_c_std, fmt = 'b.', markersize = 8)
    print (frac_flux_err_c_med, frac_flux_err_c_std)
    
    temp = (master_arr[loc_c_success[0], 2] - master_arr[loc_c_success[0], 14])/ master_arr[loc_c_success[0], 14]
    mean, frac_sigxx_err_c_med, frac_sigxx_err_c_std = sigma_clipped_stats(temp.flatten())
    axs[1, 0].errorbar(snr_mid, frac_sigxx_err_c_med, yerr = frac_sigxx_err_c_std, fmt = 'b.', markersize = 8)
    
    
    temp = (master_arr[loc_c_success[0], 4] - master_arr[loc_c_success[0], 16])/ master_arr[loc_c_success[0], 16]
    mean, frac_sigxy_err_c_med, frac_sigxy_err_c_std = sigma_clipped_stats(temp.flatten())
    axs[1, 1].errorbar(snr_mid, frac_sigxy_err_c_med, yerr = frac_sigxy_err_c_std, fmt = 'b.', markersize = 8)
    
    
    
    
    
    
    
    
    
    #For forced phot 
    loc_f_fail = np.where((master_arr[:,9] == -99))
    loc_f_success = np.where((master_arr[:,9] != 0) & (master_arr[:,9] != -99) )
    
    print((len(loc_f_fail[0]) + len(loc_f_success[0])))
    
    fail_rate_f = len(loc_f_fail[0])/(len(loc_f_fail[0]) + len(loc_f_success[0]))
    axs[0, 0].plot(snr_mid*1.025, fail_rate_f , 'r.', markersize = 8)
    
    temp = (master_arr[loc_f_success[0], 9] - master_arr[loc_f_success[0], 13])/ master_arr[loc_f_success[0], 13]
    mean, frac_flux_err_f_med, frac_flux_err_f_std = sigma_clipped_stats(temp.flatten())
    axs[0, 1].errorbar(snr_mid*1.025, frac_flux_err_f_med, yerr = frac_flux_err_f_std, fmt = 'r.', markersize = 8)
    print (frac_flux_err_f_med, frac_flux_err_f_std)
    #break
    
    temp = (master_arr[loc_f_success[0], 10] - master_arr[loc_f_success[0], 14])/ master_arr[loc_f_success[0], 14]
    mean, frac_sigxx_err_f_med, frac_sigxx_err_f_std = sigma_clipped_stats(temp.flatten())
    axs[1, 0].errorbar(snr_mid*1.025, frac_sigxx_err_f_med, yerr = frac_sigxx_err_f_std, fmt = 'r.', markersize = 8)
    
    
    temp = (master_arr[loc_f_success[0], 12] - master_arr[loc_f_success[0], 16])/ master_arr[loc_f_success[0], 16]
    mean, frac_sigxy_err_f_med, frac_sigxy_err_f_std = sigma_clipped_stats(temp.flatten())
    axs[1, 1].errorbar(snr_mid*1.025, frac_sigxy_err_f_med, yerr = frac_sigxy_err_f_std, fmt = 'r.', markersize = 8)
    
    
    
    
    
    
    
    
    
    #For forced measure
    loc_fail = np.where((master_arr[:,5] == -99))
    loc_success = np.where((master_arr[:,5] != 0) & (master_arr[:,5] != -99) )
    
    print((len(loc_fail[0]) + len(loc_success[0])))
    
    fail_rate = len(loc_fail[0])/(len(loc_fail[0]) + len(loc_success[0]))
    axs[0, 0].plot(snr_mid*1.04, fail_rate , 'k.', markersize = 8)
    
    temp = (master_arr[loc_success[0], 5] - master_arr[loc_success[0], 13])/ master_arr[loc_success[0], 13]
    mean, frac_flux_err_med, frac_flux_err_std = sigma_clipped_stats(temp.flatten())
    axs[0, 1].errorbar(snr_mid*1.05, frac_flux_err_med, yerr = frac_flux_err_std, fmt = 'k.', markersize = 8)
    print (frac_flux_err_med, frac_flux_err_std)
    
    temp = (master_arr[loc_success[0], 6] - master_arr[loc_success[0], 14])/ master_arr[loc_success[0], 14]
    mean, frac_sigxx_err_med, frac_sigxx_err_std = sigma_clipped_stats(temp.flatten())
    axs[1, 0].errorbar(snr_mid*1.05, frac_sigxx_err_med, yerr = frac_sigxx_err_std, fmt = 'k.', markersize = 8)
    print (frac_sigxx_err_med, frac_sigxx_err_std)
    
    temp = (master_arr[loc_success[0], 8] - master_arr[loc_success[0], 16])/ master_arr[loc_success[0], 16]
    mean, frac_sigxy_err_med, frac_sigxy_err_std = sigma_clipped_stats(temp.flatten())
    axs[1, 1].errorbar(snr_mid*1.05, frac_sigxy_err_med, yerr = frac_sigxy_err_std, fmt = 'k.', markersize = 8)
    print ('*************')
    
    
    
    
    
# =============================================================================
#     #For test
#     loc_fail_t = np.where((master_arr[:,21] == -99))
#     loc_success_t = np.where((master_arr[:,21] != 0) & (master_arr[:,21] != -99) )
#     
#     print((len(loc_fail_t[0]) ,len(loc_success_t[0])))
#     
#     fail_rate_t = len(loc_fail_t[0])/(len(loc_fail_t[0]) + len(loc_success_t[0]))
#     axs[0, 0].plot(snr_mid*1.075, fail_rate_t , 'g.', markersize = 10)
#     
#     temp = (master_arr[loc_success_t[0], 21] - master_arr[loc_success_t[0], 13])/ master_arr[loc_success_t[0], 13]
#     mean, frac_flux_err_t_med, frac_flux_err_t_std = sigma_clipped_stats(temp.flatten())
#     axs[0, 1].errorbar(snr_mid*1.075, frac_flux_err_t_med, yerr = frac_flux_err_t_std, fmt = 'g.', markersize = 8)
#     print (frac_flux_err_t_med, frac_flux_err_t_std)
#     
#     temp = (master_arr[loc_success_t[0], 22] - master_arr[loc_success_t[0], 14])/ master_arr[loc_success_t[0], 14]
#     mean, frac_sigxx_err_t_med, frac_sigxx_err_t_std = sigma_clipped_stats(temp.flatten())
#     axs[1, 0].errorbar(snr_mid*1.075, frac_sigxx_err_t_med, yerr = frac_sigxx_err_t_std, fmt = 'g.', markersize = 8)
#     print (frac_sigxx_err_t_med, frac_sigxx_err_t_std)
#     
#     temp = (master_arr[loc_success_t[0], 24] - master_arr[loc_success_t[0], 16])/ master_arr[loc_success_t[0], 16]
#     mean, frac_sigxy_err_t_med, frac_sigxy_err_t_std = sigma_clipped_stats(temp.flatten())
#     axs[1, 1].errorbar(snr_mid*1.075, frac_sigxy_err_t_med, yerr = frac_sigxy_err_t_std, fmt = 'g.', markersize = 8)
#     print ('*************')
# =============================================================================
    
    
    
axs[0, 0].set( ylabel = 'Fractional Failure')    
axs[0, 1].set( ylabel = 'Fractional error flux')
axs[1, 0].set( xlabel = r' SNR ', ylabel = r' Fractional error in $\sigma^2_{xx}$')
axs[1, 1].set( xlabel = r' SNR ', ylabel = r'Fractional error in $\sigma^2_{xy}$')

