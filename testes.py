#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 14:15:50 2022

@author: Bruno B Bizarria

This codes do tests with the simulated catalog to verify how much it is
similar to the original DES redmagic galaxy catalog. We obtain:
    - Normalized redshift distribution without corrector factor;
    - Normalized redshift distribution with corrector factor;
    - Mean number of galaxy per sky area;
    - Results with DES binary and non-binary masks;
    
"""
#%%
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
#%%
#Define the redshift associated with the catalogs
frequency = np.linspace(980.,1260.,31)   ## Frequency range of the BINGO 
reds = 1420.40575177/frequency - 1.
x=[]
for i in range(0,30):
    x.append(reds[i]/2+reds[i+1]/2)
print(x)
#%%
#DES MASKS CONSTRUCTION 
mask_des = hp.read_map('/home/bruno/Downloads/cross_correlations/simulador_des/mascara_des/DES_256_RING.fits')
bin_des=np.zeros(len(mask_des))
for i in range(len(mask_des)):
    if mask_des[i] == -1.6375e+30:
        bin_des[i]=0
        mask_des[i]=0
    else:
        bin_des[i]=1

#APPLICATION OF THE MASKS
nbin=12*256**2
cat0m =   np.zeros(nbin)
cat0bin = np.zeros(nbin)
cat1m =   np.zeros(nbin)
cat1bin = np.zeros(nbin)
u0=[]
v0=[]
u1=[]
v1=[]


for j in range(0,30):
    cat0=np.load('/home/bruno/Downloads/cross_correlations/catalogos_22_06/sample0/contagem/'+str(j)+'Mhz.npy')
    cat1=np.load('/home/bruno/Downloads/cross_correlations/catalogos_22_06/sample1/contagem/'+str(j)+'Mhz.npy')
    for i in range(0,len(cat0)):
        cat0m[i]   = cat0[i]*mask_des[i]
        cat1m[i]   = cat1[i]*bin_des[i]
        cat0bin[i] = cat0[i]*mask_des[i]
        cat1bin[i] = cat1[i]*bin_des[i]
    u0.append(np.sum(cat0m))
    v0.append(np.sum(cat0bin))
    u1.append(np.sum(cat1m))
    v1.append(np.sum(cat1bin))
#%%
#DES DATA 
cat = fits.open('/home/bruno/Downloads/cross_correlations/real_des/y3a2_gold2.2.1_redmagic_highdens.fits')

red    = cat[1].data['zredmagic']
erro   = cat[1].data['zredmagic_e']
red_s0 = cat[1].data['zredmagic_samp_0']
red_s1 = cat[1].data['zredmagic_samp_1']
red_s2 = cat[1].data['zredmagic_samp_2']
red_s3 = cat[1].data['zredmagic_samp_3']
bin_t  = cat[1].data['bin_number']

bin_n=bin_t
red_s0b1 = red_s0[bin_n==1]
red_s1b2 = red_s1[bin_n==2]
red_s2b3 = red_s2[bin_n==3]

rr=red_s1b2
factor=0
for i in range(0,len(rr)):
    if rr[i]<0.449:
        factor=factor+1
factor=1-1/(factor/np.sum(u1))
#%%
u0=u0/((x[1]-x[2])*np.sum(u0))
v0=v0/((x[1]-x[2])*np.sum(v0))
u1=u1/((x[1]-x[2])*np.sum(u1))*factor
v1=v1/((x[1]-x[2])*np.sum(v1))*factor
#%%
n=60

d1=np.linspace(min(red_s0b1),max(red_s0b1),n)
d2=np.linspace(min(red_s1b2),max(red_s1b2),n)
d3=np.linspace(min(red_s2b3),max(red_s2b3),n)
red1,bin1=np.histogram(red_s0b1,n)
red2,bin2=np.histogram(red_s1b2,n)
red3,bin3=np.histogram(red_s2b3,n)

red1=np.multiply(red1,1/(len(red_s0b1)*(d1[2]-d1[1])))
red2=np.multiply(red2,1/(len(red_s1b2)*(d2[2]-d2[1])))
red3=np.multiply(red3,1/(len(red_s2b3)*(d3[2]-d3[1])))

plt.plot(d1,red1,color='r',label='s0b1')
plt.plot(d2,red2,color='b',label='s1b2')
plt.plot(d3,red3,color='g',label='s2b3')
plt.show()
#%%

#PLOTS

plt.plot(x,u0,'b--',label='Ponderated Mask')
plt.plot(x,u1,'b--')
plt.plot(x,v0,'g',label='Binary Mask')
plt.plot(x,v1,'g')
plt.ylabel('p(z)')
plt.xlabel('z')
plt.legend()
plt.title('Normalized Redshift Distribution simulation')
plt.show()

plt.plot(x,u0,'b--',label='Simulation DES')
plt.plot(x,u1,'b--')
plt.plot(d1,red1,'g',label='Real DES')
plt.plot(d2,red2,'g')
plt.ylabel('p(z)')
plt.xlabel('z')
plt.legend()
plt.title('DES p(z) simulation X real data')
plt.xlim(0.127,0.449)
plt.savefig('/home/bruno/f')

#%%
'''
hp.mollview(mask_des,title='DES binary mask',coord=['E','E'])
plt.savefig('/home/bruno/Imagens/disserta/des_moll_bin.png')
#hp.mollview(mask_overlap)
hp.mollview(mask_overlap_bin, title ='DES/BINGO binary overlap mask',coord=['E','E'])
plt.savefig('/home/bruno/Imagens/disserta/overdes_moll_bin.png')
'''

#%%
full=41252.96
print(full/8)
print(np.sum(red_s1[bin_n==2])/(3600*0.15*5156.62))