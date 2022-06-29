from def_funcs import *

#Define parameters
path='/home/bruno/Downloads/cross_correlations/simulador_des/des_17_06'
path_lightcone = '/home/bruno/Downloads/cross_correlations/lightcone/HR4_fof_z01-05/fof_lightcone_lite.npy'
edge_bin=0.35
redshift_divisions=100
mass_divisions=250

#sample 0

sample=0

set_directory(path,sample)

m0,x0,y0,z0,red0=separe_data(path_lightcone,edge_bin,sample)

gera_gal(redshift_divisions,mass_divisions,m0,x0,y0,z0,red0,sample,path)

reconstruct_data(path,sample,redshift_divisions)

create_catalogues(path,sample)

#sample 1

sample=1

set_directory(path,sample)

m1,x1,y1,z1,red1=separe_data(path_lightcone,edge_bin,sample)

gera_gal(redshift_divisions,mass_divisions,m1,x1,y1,z1,red1,sample,path)

reconstruct_data(path,sample,redshift_divisions)

create_catalogues(path,sample)
