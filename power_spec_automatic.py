#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 14:10:48 2022

@author: bruno
"""

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import pymaster as nmt
#%%
mask = hp.read_map('/home/bruno/Downloads/cross_correlations/mascara_bingo/mask_bingoGal_5yr_ns256_28h2d_apod.fits') 
nside = 256
#%%
## HI IM REAL
#cat1_HI = np.load("/home/bruno/Downloads/cross_correlations/hr4_jiajun/21cmTb_real/0.npy")
#cat2_HI = np.load("/home/bruno/Downloads/cross_correlations/hr4_jiajun/21cmTb_real/7.npy")
cat3_HI = np.load("/home/bruno/Downloads/cross_correlations/hr4_jiajun/21cmTb_real/14.npy")
# MORE 2 HI
#cat4_HI = np.load("/home/bruno/Downloads/cross_correlations/hr4_jiajun/21cmTb_real/22.npy")
#cat5_HI = np.load("/home/bruno/Downloads/cross_correlations/hr4_jiajun/21cmTb_real/29.npy")

path='/media/bruno/USB DISK/sample1/densidade/'
## Galaxy CATALOG REAL
#cat1_gal = np.load(path+"0Mhz.npy")
#cat2_gal = np.load(path+"7Mhz.npy")
cat3_gal = np.load(path+"14Mhz.npy")
# MORE 2 GALAXIES CATALOGS
#cat4_gal = np.load(path+"22Mhz.npy")
#cat5_gal = np.load(path+"29Mhz.npy")

freq_bins=[980,1045,1111,1185,1260]
#%%apagar gambiarra
path0='/media/bruno/USB DISK/sample0/densidade/'
path1='/media/bruno/USB DISK/sample1/densidade/'

cat0=np.load(path0+"14Mhz.npy")
cat1=np.load(path1+"14Mhz.npy")
cat3_gal=[]
for i in range(0,len(cat0)):
    cat3_gal.append(cat0[i]+cat1[i])


#%%
def power_spectruns(mask,nside,b_mult,cat_g,cat_h,freq):
    #tirar a média
    g = np.full(len(cat_g),np.mean(cat_g))
    cat_g = np.subtract(cat_g, g)
    
    h = np.full(len(cat_h),np.mean(cat_h))
    cat_h = np.subtract(cat_h, h)
    
    fh1=nmt.NmtField(mask,[cat_h])
    fg1=nmt.NmtField(mask,[cat_g])

    b = nmt.NmtBin.from_nside_linear(nside, b_mult)
    
    cl_auto_fh1 = nmt.compute_full_master(fh1, fh1, b)[0]
    cl_auto_fg1 = nmt.compute_full_master(fg1, fg1, b)[0]
    cl_cros_fh1_fg1 = nmt.compute_full_master(fh1, fg1, b)[0]

    ell_arr = b.get_effective_ells()
    
    u=[]
    for i in range(len(cl_auto_fh1)):
        u.append((cl_cros_fh1_fg1[i]/np.sqrt(cl_auto_fg1[i]*cl_auto_fh1[i])))

    plt.plot(ell_arr,u)
    plt.loglog()
    plt.ylabel('Cl_X/sqrt(Cl_gal * Cl_HI)')
    plt.xlabel('l')
    plt.title('Correlation Coefficient '+str(freq)+'MHz Masked map')
    plt.ylim(ymax=1,ymin=0.25)
    plt.show()

    plt.plot(ell_arr, cl_auto_fh1, 'r-', label='Auto HI_IM')
    plt.plot(ell_arr, cl_auto_fg1, 'b-', label='Auto gal_dens')
    plt.plot(ell_arr, cl_cros_fh1_fg1, 'y-', label='Cross GxH')
    plt.loglog()
    plt.xlabel('$\\ell$', fontsize=16)
    plt.ylabel('$C_\\ell$', fontsize=16)
    plt.title(str(freq)+' MHz masked')
    plt.legend(loc='upper right', ncol=2, labelspacing=0.1)
    plt.ylim(ymax=10**(-9),ymin=10**(-14))
    plt.show()
    
    ana_c1H1 = hp.sphtfunc.anafast(cat_h, nspec=None, lmax=None, mmax=None, iter=3, alm=False, pol=True, use_weights=False, datapath=None, gal_cut=0, use_pixel_weights=False)
    ana_c1gal = hp.sphtfunc.anafast(cat_g, nspec=None, lmax=None, mmax=None, iter=3, alm=False, pol=True, use_weights=False, datapath=None, gal_cut=0, use_pixel_weights=False)
    ana_H1gal1 = hp.sphtfunc.anafast(cat_h,cat_g, nspec=None, lmax=None, mmax=None, iter=3, alm=False, pol=True, use_weights=False, datapath=None, gal_cut=0, use_pixel_weights=False)

    l1h=list(range(len(ana_c1H1))) 
    l1g=list(range(len(ana_c1gal))) 
    l1gh=list(range(len(ana_H1gal1))) 
    
    plt.plot(l1h,ana_c1H1, 'r-',label='HI')
    plt.plot(l1g,ana_c1gal, 'b-',label='gal')
    plt.plot(l1gh,ana_H1gal1, 'y-',label='HI X gal')
    plt.loglog()
    plt.xlabel('$l$', fontsize=16)
    plt.ylabel('$C_l$', fontsize=16)
    plt.ylim(10**-14)
    #plt.title(str(freq)+' MHz full sky')
    plt.title(str(freq)+' MHz masked map')
    plt.legend(loc='upper right', ncol=2, labelspacing=0.1)
    plt.show()
    #plt.savefig('/home/bruno/Downloads/cross_correlations/power_spec.pdf')
    
    v=[]
    for i in range(len(ana_H1gal1)):
        v.append((ana_H1gal1[i]/np.sqrt(ana_c1H1[i]*ana_c1gal[i])))
        
    plt.plot(l1gh,v)
    plt.loglog()
    plt.ylabel('r(l)')
    plt.xlabel('l')
    #plt.title('Correlation Coefficient '+str(freq)+'MHz full sky map')
    plt.title('Correlation Coefficient '+str(freq)+'MHz masked map')
    plt.show()
    #plt.savefig('/home/bruno/Downloads/cross_correlations/coef_cor.pdf')
#%%
power_spectruns(mask,nside,10,cat3_gal,cat3_HI,freq_bins[2])
'''
power_spectruns(mask,nside,10,cat2_gal,cat2_HI,freq_bins[1])
power_spectruns(mask,nside,10,cat3_gal,cat3_HI,freq_bins[2])
power_spectruns(mask,nside,10,cat4_gal,cat4_HI,freq_bins[3])
power_spectruns(mask,nside,10,cat5_gal,cat5_HI,freq_bins[4])
'''