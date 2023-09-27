import matplotlib
import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import MFLensingTools as mf

nside =512
lmax = 3*nside-1
lmin = 8
threshold_num = 31
smooth_angle = 50

As = 2.1 # 1e-9
wcdm = 0.12
Mnu = 0.06

fnl = 3
gnl = 1


ls_camb,cl_camb = mf.camb_wl.Cl_camb(As,wcdm,Mnu,lmax,lmin=lmin)
beam_l = hp.sphtfunc.gauss_beam(np.radians(smooth_angle/60), lmax=lmax)
cl_camb_smoothed = cl_camb*beam_l**2

map_G = hp.sphtfunc.synfast(cl_camb, nside,lmax=lmax)
map_G = hp.pixelfunc.remove_dipole(map_G)
map_G = hp.sphtfunc.smoothing(map_G, fwhm=np.radians(smooth_angle/60)) #100 arcmin
alm_G = hp.sphtfunc.map2alm(map_G,lmax=lmax)

map_nG = mf.generate_nG_map.generate_nG_map(map_G,fnl=fnl,gnl=gnl)
alm_nG = hp.sphtfunc.map2alm(map_nG,lmax=lmax)

sigma_sk_nG = mf.SK_parameters.S_K_parameter(alm_nG,nside,lmax)
sigma_sk_G = mf.SK_parameters.S_K_parameter(alm_G,nside,lmax)

[sigma_nG,sigma1_nG,S_nG,S1_nG,S2_nG,K_nG,K1_nG,K12_nG,K22_nG] = sigma_sk_nG
[sigma_G,sigma1_G,S_G,S1_G,S2_G,K_G,K1_G,K12_G,K22_G] = sigma_sk_G

print('sigma_0:','Gaussian=',sigma_G,'non-Gaussian=',sigma_nG)
print('sigma_1:','Gaussian=',sigma1_G,'non-Gaussian=',sigma1_nG)
print('S:','Gaussian=',S_G,'non-Gaussian=',S_nG)
print('S1:','Gaussian=',S1_G,'non-Gaussian=',S1_nG)
print('S2:','Gaussian=',S2_G,'non-Gaussian=',S2_nG)
print('K:','Gaussian=',K_G,'non-Gaussian=',K_nG)
print('K1:','Gaussian=',K1_G,'non-Gaussian=',K1_nG)
print('K12:','Gaussian=',K12_G,'non-Gaussian=',K12_nG)
print('K22:','Gaussian=',K22_G,'non-Gaussian=',K22_nG)


[v0_cl,v1_cl,v2_cl] = mf.theory_MFs.analytical_MFs(cl_camb_smoothed,lmin,lmax,3*sigma_nG,threshold_num)
[v0_map_G,v1_map_G,v2_map_G] = mf.get_MFs.MFs_map(nside,alm_G,3*sigma_nG,threshold_num,lmax,width_divide=10)
[v0_map_nG,v1_map_nG,v2_map_nG] = mf.get_MFs.MFs_map(nside,alm_nG,3*sigma_nG,threshold_num,lmax,width_divide=10)

[v0_corr,v1_corr,v2_corr] = mf.nG_correction.MF_correction(3*sigma_nG,threshold_num,sigma_sk=sigma_sk_nG,order=2)

v0_nG_corr = v0_cl + v0_corr
v1_nG_corr = v1_cl + v1_corr
v2_nG_corr = v2_cl + v2_corr

threshold = np.linspace(-3, 3, num=threshold_num)
v_names = ['V0', 'V1', 'V2']
v_latex = [r'$V_0$', r'$V_1$', r'$V_2$']
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12,6), sharex=True)
for i, v in enumerate(v_names):
    axs[0,i].plot(threshold, eval(f'{v.lower()}_map_nG'), label=f'{v} map (nG)')
    axs[0,i].plot(threshold, eval(f'{v.lower()}_nG_corr'), label=f'{v} analytical\n(corrected)')
    axs[0,i].plot(threshold, eval(f'{v.lower()}_cl'), label=f'{v} analytical')
    axs[0,i].set_ylabel(f'{v_latex[i]}')

axs[0,2].legend(loc='upper left', bbox_to_anchor=(1,1))
for i, v in enumerate(v_names):
    axs[1,i].plot(threshold, eval(f'{v.lower()}_map_nG-{v.lower()}_nG_corr'), label=f'{v} analytical\n(corrected)')
    axs[1,i].plot(threshold, eval(f'{v.lower()}_map_nG-{v.lower()}_cl'), label=f'{v} analytical')
    axs[1,i].set_ylabel(f'{v_latex[i]} nG - {v_latex[i]} analytical')
    axs[1,i].set_xlabel(r'threshold (in $\sigma_0$)')
axs[1,2].legend(loc='upper left', bbox_to_anchor=(1,1))

fig.tight_layout()
fig.subplots_adjust(right=0.85)
plt.savefig('MF_sk_corr_nG_map')
plt.close()
