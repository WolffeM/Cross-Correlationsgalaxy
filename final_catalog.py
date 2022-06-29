#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 22:01:47 2022

@author: bruno
"""

import numpy as np

k=0
while k<30:
    u_dens_1=np.load("/home/bruno/Downloads/cross_correlations/simulador_des/catalogos/z_menor_035_des/densidade_light/"+str(k)+"Mhz.npy")
    u_dens_2=np.load("/home/bruno/Downloads/cross_correlations/simulador_des/catalogos/z_maior_035_des/densidade_light/"+str(k)+"Mhz.npy")
    
    u_dens=[ (a + b) for a, b in zip(u_dens_1, u_dens_2) ]

    np.save("/home/bruno/Downloads/cross_correlations/simulador_des/catalogos/catalogos_finais/densidade/"+str(k)+"Mhz.npy",u_dens)
    print('catalogo densidade '+str(k)+' pronto')
    
    u_cont_1=np.load("/home/bruno/Downloads/cross_correlations/simulador_des/catalogos/z_menor_035_des/contagem_light/"+str(k)+"Mhz.npy")
    u_cont_2=np.load("/home/bruno/Downloads/cross_correlations/simulador_des/catalogos/z_maior_035_des/contagem_light/"+str(k)+"Mhz.npy")
    
    u_cont=[ (a + b) for a, b in zip(u_cont_1, u_cont_2) ]

    np.save("/home/bruno/Downloads/cross_correlations/simulador_des/catalogos/catalogos_finais/contagem/"+str(k)+"Mhz.npy",u_cont)
    print('catalogo contagem '+str(k)+' pronto')
    
    k=k+1