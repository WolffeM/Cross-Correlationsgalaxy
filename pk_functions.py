#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 10:59:20 2022

@author: Bruno Bizarria

Functions to obtain 3D field and calculate power spectruns 
"""
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.cosmology import FlatLambdaCDM, z_at_value


def Field3d(Om0,H0,freq_min,freq_max,number_maps,nside,box_size,path): 
    
    cosmo = FlatLambdaCDM(Om0=Om0,H0=H0)
    frequency = np.linspace(freq_min,freq_max,number_maps+1)
    redshift = 1420.40575177/frequency - 1.
    distance = cosmo.comoving_distance(redshift).value * cosmo.h
    
    nmesh = nside
    boxsize = box_size
    x_corr = np.linspace(0,boxsize-boxsize/nmesh,nmesh) + 0.5*boxsize/nmesh - boxsize/2.
    y_corr = np.linspace(0,boxsize-boxsize/nmesh,nmesh) + 0.5*boxsize/nmesh - boxsize/2.
    z_corr = np.linspace(0,boxsize-boxsize/nmesh,nmesh) + 0.5*boxsize/nmesh - boxsize/2.
    
    x_grid = np.zeros([nmesh,nmesh,nmesh],dtype=np.float32)
    y_grid = np.zeros([nmesh,nmesh,nmesh],dtype=np.float32)
    z_grid = np.zeros([nmesh,nmesh,nmesh],dtype=np.float32)
    d_grid = np.zeros([nmesh,nmesh,nmesh],dtype=np.float32)
    
    for i in np.arange(0,nmesh,1):
        for j in np.arange(0,nmesh,1):
            x_grid[i,j,:] = x_corr
            y_grid[i,:,j] = y_corr
            z_grid[:,i,j] = z_corr
    
    r_grid = np.sqrt(x_grid**2.+y_grid**2.+z_grid**2.)
    p_grid = hp.vec2pix(256,x_grid,y_grid,z_grid)
    
    for i in np.arange(0,number_maps,1):
        real_tb = np.load(path+str(i)+'.npy')
        delta_tb = real_tb/np.average(real_tb) - 1.
        np.float32(delta_tb)
        r_filter = (r_grid > distance[i+1]) & (r_grid <= distance[i])
        d_grid = d_grid + delta_tb[p_grid] * r_filter
    print('Full sky 3D field calculated')
    
    return d_grid

def Field3d_masked(Om0,H0,freq_min,freq_max,number_maps,nside,box_size,path,path_mask):
    
    cosmo = FlatLambdaCDM(Om0=Om0,H0=H0)
    frequency = np.linspace(freq_min,freq_max,number_maps+1)
    redshift = 1420.40575177/frequency - 1.
    distance = cosmo.comoving_distance(redshift).value * cosmo.h
    
    nmesh = nside
    boxsize = box_size
    x_corr = np.linspace(0,boxsize-boxsize/nmesh,nmesh) + 0.5*boxsize/nmesh - boxsize/2.
    y_corr = np.linspace(0,boxsize-boxsize/nmesh,nmesh) + 0.5*boxsize/nmesh - boxsize/2.
    z_corr = np.linspace(0,boxsize-boxsize/nmesh,nmesh) + 0.5*boxsize/nmesh - boxsize/2.
    
    x_grid = np.zeros([nmesh,nmesh,nmesh],dtype=np.float32)
    y_grid = np.zeros([nmesh,nmesh,nmesh],dtype=np.float32)
    z_grid = np.zeros([nmesh,nmesh,nmesh],dtype=np.float32)
    d_grid = np.zeros([nmesh,nmesh,nmesh],dtype=np.float32)
    m_grid = np.zeros([nmesh,nmesh,nmesh],dtype=np.float32)
 
    for i in np.arange(0,nmesh,1):
        for j in np.arange(0,nmesh,1):
            x_grid[i,j,:] = x_corr
            y_grid[i,:,j] = y_corr
            z_grid[:,i,j] = z_corr
    
    r_grid = np.sqrt(x_grid**2.+y_grid**2.+z_grid**2.)
    p_grid = hp.vec2pix(256,x_grid,y_grid,z_grid)
 

    for i in np.arange(0,number_maps,1):
        mask = hp.read_map(path_mask)
        real_tb = np.load(path+str(i)+'.npy')
        if np.average(real_tb) == 0:
            delta_tb=real_tb*0
            print('One of your maps was zero, but I accepted')
        delta_tb = real_tb/np.average(real_tb) - 1
        for j in range(0,len(mask)):
            if mask[j]<=0:
                delta_tb[j]=0
        np.float32(delta_tb)
        r_filter = (r_grid > distance[i+1]) & (r_grid <= distance[i])
        m_grid = m_grid + delta_tb[p_grid] * r_filter
    print('Masked 3D field calculated')
    
    return m_grid

def BAO_profile(Pk,k,nmesh,nonzero_grid,path_no_wiggle,title_plot,path_savefig,name_fig,dpi):
  pk = Pk * nmesh**3. / nonzero_grid
  
  kk = np.load(path_no_wiggle+'k.npy')
  Plin = np.load(path_no_wiggle+'pk.npy')
  
  
  p_no_bao=[]
  for j in range(0,len(k)):
    n=0
    for i in range(0,len(kk)-1):
      if kk[i]<k[j] and kk[i+1]>k[j]:
        n=i
    a = ((Plin[n+1]-Plin[n])/(kk[n+1]-kk[n]))*(k[j]-kk[n]) + Plin[n]
    p_no_bao.append(a)
  
  profile_bao = []
  for i in range(0,len(k)):
    profile_bao.append(pk[i]/p_no_bao[i])

  plt.loglog(k,profile_bao)
  plt.xlabel('k(h/Mpc)')
  plt.xlim(0.04,0.15)
  plt.axvline(0.06)
  plt.axvline(0.11)
  plt.grid()
  plt.title(title_plot)
  plt.savefig(path_savefig+name_fig,dpi=dpi)
  
  return profile_bao
  
def plot_pk2d(Pk2D,kpar,kper,title):
  intermed = []
  for i in range(0,len(kper)):
    if kpar[i]==kper[i]:
      intermed.append(i)
  kpar_max = intermed[1]-intermed[0]-1
  kper_max = round(len(kper)/kpar_max)
  print('kpar max: '+str(kpar_max))
  print('kper max: '+str(kper_max))
  
  pk2d = Pk2D.reshape(kper_max,kpar_max)
  
  plt.xlabel(r'$k_{\perp}$')
  plt.ylabel(r'$k_{||}$')
  plt.imshow(pk2d,cmap='gist_rainbow',norm=LogNorm(),origin='lower',extent=[min(kper),max(kper),min(kpar),max(kpar)])
  cbar=plt.colorbar()
  cbar.set_label(r'$P_{HI}(k_{||},k_{\perp})$')
  plt.title(title)
