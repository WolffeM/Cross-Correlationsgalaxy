#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 14:14:29 2022

@author: Bruno Benedito Bizarria
email:bruno.bizarria@inpe.br


This code takes an halo lightcone and occupes it with galaxies, following a hod obtained from a survey observation. We are using the hod obtained by DES in order to simulate its galaxy catalog. The photometric redshift error is applied, to obtain a more realistic catalog.
here we describe the code in steps:
STEP 1:
set the directories we will need to save the auxiliar files created during the execution.
STEP 2:
DES has different hod's for different redshift interval os his samples. So we need to separe the halos by redshift intervals. you can choose sample 0 or sample 1. the first takes z<0.35 and the former z>0.35. as the bingo's redshift range has its maximum at 0.45, we will not deal with more des samples. so, we need to separe the data
STEP 3:
we have to separe the halo by mass, because the hod is a mass function. as we have the mean number of galaxies ng inside a halo of mass m, and we have n halos in the interval m+dm, we need to generate a poisson distribution with mean ng for n objects, and then associate each one of the numbers generated in this poisson distribution with an halo in the sky. However, there are a lot of halos (10^9), and implement this process consumes a lot of RAM memory. So we divide the lightcone in redshift intervals and apply the algorithm for each interval. Adding, we apply the error for each galaxy generated
STEP 4:
we must reconstruct the total lightcone, because it was divided into redshift intervals
STEP 5:
with the lightcone reconstructed, we must divide it in 30 redshift channels, and create number counts catalogs and count density catalogs.
"""

#Lybraries

import os
import shutil
import numpy as np
import healpy as hp
from astropy.cosmology import FlatLambdaCDM ## Cosmological tools of astropy
import math

#STEP 1

def set_directory(path_file,sample):
    path_sample=path_file+'/sample'+str(sample)
    
    
    if os.path.exists(path_sample+'.log'):
      os.remove(path_sample+'.log')
    
    file_log = open(path_sample+'.log','a')
    
    shutil.rmtree(path_sample, ignore_errors = True)
    
    os.mkdir(path_sample)
    os.mkdir(path_sample+'/data_arrays')
    os.mkdir(path_sample+'/densidade')
    os.mkdir(path_sample+'/contagem')
    print('Directories set')
    return

#STEP 2

def separe_data(path_lightcone,div,sample):

    cosmo = FlatLambdaCDM(Om0=0.26,H0=72)
    
    #prepare for loading the data
    data = np.load(path_lightcone)
    '''
    #Run with the total lightcone
    m = data[:,0] + np.log10(cosmo.h)
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]
    red = data[:,4]
    
    '''
    #Run with the partial lightcone
    x = data[1,:]
    y = data[2,:]
    z = data[3,:]
    m = data[0,:] + np.log10(cosmo.h)
    red = data[4,:]
    
    
    x1=[]
    y1=[]
    z1=[]
    m1=[]
    red1=[]
    ft = [0 if w<=div else 1 for w in red]
    
    for j in range(0,len(red)):
        if ft[j]==sample:
            x1.append(x[j])
            y1.append(y[j])
            z1.append(z[j])
            m1.append(m[j])
            red1.append(red[j])
    print('Data separeted for sample'+str(sample))
    return m1,x1,y1,z1,red1

def hod_s0(b):      
    mmin1=11.97
    mm1=13.51
    slogm1=0.26
    alpha1=1.88
    fcen1=0.12
    n2 = (fcen1/2)*(1+math.erf((b-mmin1)/slogm1))*(1+(b/(10**mm1))**alpha1)
    return n2

def hod_s1(b):
    mmin2=12.13
    m2=13.79
    slogm2=0.34
    alpha2=2.61
    fcen2=0.13
    n2 = (fcen2/2)*(1+math.erf((b-mmin2)/slogm2))*(1+(b/(10**m2))**alpha2)
    return n2

#STEP 3
#this step can be modified to do a parallel computation.

def gera_gal(red_divisions,mass_bins,m1,x1,y1,z1,red1,sample,path_file):
    path=path_file+'/sample'+str(sample)+'/data_arrays'
    
    if sample==0:
        hod=hod_s0
    elif sample==1:
        hod=hod_s1
    
    u = np.linspace(0.12,0.45,red_divisions)
    
    for i in range(0,red_divisions-1):
        x=[]
        y=[]
        z=[]
        m=[]
        red=[]
        ft = [1 if u[i]<=w<=u[i+1] else 0 for w in red1]
        
        for j in range(0,len(red1)):
            if ft[j]==1:
                x.append(x1[j])
                y.append(y1[j])
                z.append(z1[j])
                m.append(m1[j])
                red.append(red1[j])
                                        
        hist_m,bins1=np.histogram(m,mass_bins)  #hist_m is the number of halos in a bins1 interval 
        
        aux=[]                                  #aux is an array of poisson arrays
        for k in range(0,len(hist_m)):
            aux.append(np.random.poisson(hod(bins1[k]),hist_m[k]))
            
        inter=bins1[2]-bins1[1]
         
        #this part of the algorithm puts an label in each halo. If the halo is in the bin 1, the it receives a correspondent label.
        #It will help to redistribute the galaxies generated, that are in the aux array.    
        #uu = list(map(lambda n: math.floor(n/inter), m))
        uu = list(map(lambda n: n/inter, m))
        hist_label,b=np.histogram(uu,mass_bins)
        numero_gal=[]
        
        for j in range(0,len(m)):
            for k in range(0,len(b)-1):
                if ((b[k]<=uu[j]<=b[k+1]) and(len(aux[k])!=0)): #ERRO ESTA AQUI
                    numero_gal.append(aux[k][0])
                    aux[k]=np.delete(aux[k],[0])
        print('galaxies generated for division '+str(i))
        
        #We must account for halos that has more than one galaxy. The way that we do is to duplicate the coordinates of the host halo. Example: A halo hosts 2 galaxies, and its position are (xh,yh,zh,redh). Then we will have two galaxies with the positions (xh,yh,zh,redh). It helps us to apply the errors, because even though the two galaxies are in the same postion, the erros will be different, because they are different objects
        mf=[]
        xf=[]
        yf=[]
        zf=[]
        redf=[]
        ngf=[]
                  
        for j in range(0,len(red)):
            if numero_gal[j]!=0 and numero_gal[j]!=1:
                for l in range(0,numero_gal[j]):
                    mf.append(m[j])
                    xf.append(x[j])
                    yf.append(y[j])
                    zf.append(z[j])
                    redf.append(red[j])
            elif numero_gal[j]==1:
                mf.append(m[j])
                xf.append(x[j])
                yf.append(y[j])
                zf.append(z[j])
                redf.append(red[j])
        
        print('galaxies expansion for division '+str(i))
        print('massa'+str(len(mf)))
        for j in range(0,len(redf)):
            redf[j]=redf[j]+0.0126*(1+redf[j])*np.random.normal(0,1)
        print('error applied for division '+str(i)+'\n')    
        light_lite=np.transpose([mf,xf,yf,zf,redf])
        print(light_lite.shape)
        np.save(path+'/intermed_cat'+str(i)+'.npy',light_lite)
        
    return


#STEP 4

def reconstruct_data(path_file,sample,red_divisions):
    print('\nReconstructing data\n')
    m=[]
    x=[]
    y=[]
    z=[]
    red=[]
    for i in range(0,red_divisions-1):
        data=np.load(path_file+'/sample'+str(sample)+'/data_arrays/intermed_cat'+str(i)+'.npy')
        if len(data[:,0]) != 0:
            m1=data[:,0]
            x1=data[:,1]
            y1=data[:,2]
            z1=data[:,3]
            red1=data[:,4]
            print('data loaded '+str(i))
        
            for j in range(0,len(m1)):
                m.append(m1[j])
                x.append(x1[j])
                y.append(y1[j])
                z.append(z1[j])
                red.append(red1[j])

    catalogo=np.transpose([m,x,y,z,red])
    np.save(path_file+'/sample'+str(sample)+'/galaxy_catalog.npy',catalogo)
    

    return

#STEP 5

def create_catalogues(path_file,sample):
    cosmo = FlatLambdaCDM(Om0=0.26,H0=72)
    data=np.load(path_file+'/sample'+str(sample)+'/galaxy_catalog.npy',allow_pickle=True)
    
    m=data[:,0]
    x=data[:,1]
    y=data[:,2]
    z=data[:,3]
    red=data[:,4]
    
    print('data recovered')

        
    frequency = np.linspace(980.,1260.,31)   ## Frequency range of the BINGO 
    redshift = 1420.40575177/frequency - 1.  ## Redshift compare with the HI frequency
    #distance = cosmo.comoving_distance(redshift).value * cosmo.h ## The distance at a given redshift, without unitis,
                                                                 ## times the Hubble constant of FlatLambdaCDM
    #os.mkdir(path_catalogues+'simu'+str(i))
  
    k=0
    while(k<30):
        rx=[]
        ry=[]
        rz=[]
        for j in range(0,len(red)):
            if red[j]<redshift[k] and red[j]>redshift[k+1]:
                rx.append(x[j])
                ry.append(y[j])
                rz.append(z[j])

    	# vec2pix function translate the vector x,y,z into the location in healpix pixel
        rpix = hp.vec2pix(256,rx,ry,rz)
        npix = hp.pixelfunc.nside2npix(256)
        # the pixel location histogram calculate the galaxy number in each pixel, divided by the volume in each pixel, which is the pixel area times the depth  
        rmap = np.histogram(rpix,bins=np.arange(0,npix+1,1))[0]
        volume = cosmo.comoving_volume(redshift[k]) - cosmo.comoving_volume(redshift[k+1]) ## Deep of pixel in redshift
        volume = volume.value * cosmo.h**3./hp.nside2npix(256)                             ## Divided by area?
        rrho = rmap/volume 
        np.save(path_file+'/sample'+str(sample)+'/densidade/'+str(k)+'Mhz.npy',rrho)
        np.save(path_file+'/sample'+str(sample)+'/contagem/'+str(k)+'Mhz.npy',rmap)
        k=k+1
        print('catalogos '+str(k)+' prontos')
                

    return







