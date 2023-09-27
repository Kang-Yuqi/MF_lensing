import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import sys
import MFLensingTools as mf

nside = 512
lmax = 400
lmin = 0

As = 2.1 # 1e-9
wcdm = 0.12
Mnu = 0.06

fnl = 3
gnl = -1

ls_camb,cl_camb = mf.camb_wl.Cl_camb(As,wcdm,Mnu,lmax)

map_G = hp.sphtfunc.synfast(cl_camb, nside, lmax=lmax)
alm_G = hp.sphtfunc.map2alm(map_G,lmax=lmax)
map_nG = mf.generate_nG_map.generate_nG_map(map_G,fnl=fnl,gnl=gnl)
alm_nG = hp.sphtfunc.map2alm(map_nG,lmax=lmax)

[sigma_G,sigma1_G,S_G,S1_G,S2_G,K_G,K1_G,K12_G,K22_G,sigma_sq_map_G,sigma_1_sq_map_G,S_map_G,S1_map_G,S2_map_G,K_map_G,K1_map_G,K12_map_G,K22_map_G] = mf.SK_parameters.S_K_parameter(alm_G,nside,lmax,local=True)
[sigma_nG,sigma1_nG,S_nG,S1_nG,S2_nG,K_nG,K1_nG,K12_nG,K22_nG,sigma_sq_map_nG,sigma_1_sq_map_nG,S_map_nG,S1_map_nG,S2_map_nG,K_map_nG,K1_map_nG,K12_map_nG,K22_map_nG] = mf.SK_parameters.S_K_parameter(alm_nG,nside,lmax,local=True)

max_map = max(map_G)
min_map = min(map_G)
hp.projview(map_G, title='Gaussian map', norm="symlog2", sub=211,max=max_map,min=min_map)
hp.projview(map_nG, title='non-Gaussian map', norm="symlog2", sub=212,max=max_map,min=min_map)
plt.savefig('G_nG_map')
plt.close()
hp.projview(map_nG-map_G, title='non-Gaussian map - Gaussian map', norm="symlog2")
plt.savefig('G_nG_map_residual')
plt.close()

hp.projview(S_map_G, title=r'$S$ local map', norm="symlog2", sub=421)
hp.projview(S1_map_G, title=r'$S_1$ local map', norm="symlog2", sub=422)
hp.projview(S2_map_G, title=r'$S_2$ local map', norm="symlog2", sub=423)
hp.projview(K_map_G, title=r'$K$ local map', norm="symlog2", sub=424)
hp.projview(K1_map_G, title=r'$K_1$ local map', norm="symlog2", sub=425)
hp.projview(K12_map_G, title=r'$K_{12}$ local map', norm="symlog2", sub=426)
hp.projview(K22_map_G, title=r'$K_{22}$ local map', norm="symlog2", sub=427)
plt.savefig('G_map_sk_local_map')
plt.close()

hp.projview(S_map_nG, title=r'$S$ local map', norm="symlog2", sub=421)
hp.projview(S1_map_nG, title=r'$S_1$ local map', norm="symlog2", sub=422)
hp.projview(S2_map_nG, title=r'$S_2$ local map', norm="symlog2", sub=423)
hp.projview(K_map_nG, title=r'$K$ local map', norm="symlog2", sub=424)
hp.projview(K1_map_nG, title=r'$K_1$ local map', norm="symlog2", sub=425)
hp.projview(K12_map_nG, title=r'$K_{12}$ local map', norm="symlog2", sub=426)
hp.projview(K22_map_nG, title=r'$K_{22}$ local map', norm="symlog2", sub=427)
plt.savefig('nG_map_sk_local_map')
plt.close()


## Planck data
module_dir = '../Planck_lensing_data/'
mask_file = '../Planck_lensing_data/baseline_MV/mask.fits'
smooth_angle = 50

[ls_bias,all_bias_noise] = mf.Planck_data.get_noise_bias(lmax=lmax,lmin=lmin,dirc=module_dir)
alm_planck = mf.Planck_data.get_wl(lmax=lmax,lmin=lmin,nside=nside,dirc=module_dir,mask=False,mask_file=None,return_type='alm')
alm_planck_smoothed = hp.sphtfunc.smoothalm(alm_planck, fwhm=np.radians(smooth_angle/60))

[sigma_plk,sigma1_plk,S_plk,S1_plk,S2_plk,K_plk,K1_plk,K12_plk,K22_plk,sigma_sq_map_plk,sigma_1_sq_map_plk,S_map_plk,S1_map_plk,S2_map_plk,K_map_plk,K1_map_plk,K12_map_plk,K22_map_plk] = mf.SK_parameters.S_K_parameter(alm_planck_smoothed,nside,lmax,local=True,mask_file=mask_file)

print('sigma_0:','Gaussian=',sigma_G,'non-Gaussian=',sigma_nG,'Planck=',sigma_plk)
print('sigma_1:','Gaussian=',sigma1_G,'non-Gaussian=',sigma1_nG,'Planck=',sigma1_plk)
print('S:','Gaussian=',S_G,'non-Gaussian=',S_nG,'Planck=',S_plk)
print('S1:','Gaussian=',S1_G,'non-Gaussian=',S1_nG,'Planck=',S1_plk)
print('S2:','Gaussian=',S2_G,'non-Gaussian=',S2_nG,'Planck=',S2_plk)
print('K:','Gaussian=',K_G,'non-Gaussian=',K_nG,'Planck=',K_plk)
print('K1:','Gaussian=',K1_G,'non-Gaussian=',K1_nG,'Planck=',K1_plk)
print('K12:','Gaussian=',K12_G,'non-Gaussian=',K12_nG,'Planck=',K12_plk)
print('K22:','Gaussian=',K22_G,'non-Gaussian=',K22_nG,'Planck=',K22_plk)

hp.projview(S_map_plk, title=r'$S$ local map', norm="symlog2", sub=421)
hp.projview(S1_map_plk, title=r'$S_1$ local map', norm="symlog2", sub=422)
hp.projview(S2_map_plk, title=r'$S_2$ local map', norm="symlog2", sub=423)
hp.projview(K_map_plk, title=r'$K$ local map', norm="symlog2", sub=424)
hp.projview(K1_map_plk, title=r'$K_1$ local map', norm="symlog2", sub=425)
hp.projview(K12_map_plk, title=r'$K_{12}$ local map', norm="symlog2", sub=426)
hp.projview(K22_map_plk, title=r'$K_{22}$ local map', norm="symlog2", sub=427)
plt.savefig('Planck_map_sk_local_map_smooth_50deg')
plt.close()
