import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import sys
import MFLensingTools as mf

B = 1.95  # The parameter B, a positive number controlling the tiling of needlets
j_min = 4  # Minimum scale index
j_max = 10  # Maximum scale index
nside = 512
lmax = 400
lmin = 8
smooth_angle = 50 # arcmin

As = 2.119 # 1e-9
wcdm = 0.1203
Mnu = 0.06

ls = np.arange(lmax + 1)
ls_camb,cl_camb = mf.camb_wl.Cl_camb(As,wcdm,Mnu,lmax,lmin=lmin)
map_G = hp.sphtfunc.synfast(cl_camb, nside, lmax=lmax)
map_G = hp.sphtfunc.smoothing(map_G, fwhm=np.radians(smooth_angle / 60))
alm_G = hp.sphtfunc.map2alm(map_G, lmax=lmax)

hp.projview(map_G, title= 'sample map' , 
                norm="symlog2", cmap='planck')
plt.savefig('orgin_map')
plt.close()


needlet_window_list = []
i = 0
for j in range(j_min, j_max):
    needlet_window,filtered_alm = mf.needlets.needlet_filter(nside, j, B=B,lmax =lmax,alm=alm_G)
    needlet_map = hp.sphtfunc.alm2map(filtered_alm, nside)
    needlet_window_list.append(needlet_window)
    hp.projview(needlet_map, title= r'$B = %s$, $j=%s$' %(B,j), 
                norm="symlog2", sub=321+i, cmap='planck')
    i +=1
plt.savefig('camb_needlets_filter_map')
plt.close()


plt.figure(figsize=(10, 6))
i = 0
for j in range(j_min, j_max):
    plt.plot(ls, needlet_window_list[i], label=r"$j=%s$" %j)
    i +=1
plt.xlabel(r"multipole $\ell$")
plt.ylabel(r"needlet window function $b$")
plt.legend()
plt.grid()
plt.savefig('needlets_B_195')
# plt.savefig('needlets_B_195.pdf', format='pdf')
plt.close()
