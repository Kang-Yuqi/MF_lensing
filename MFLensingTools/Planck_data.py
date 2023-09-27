import numpy as np

def get_wl(lmax=None,lmin=0,nside=None,dirc=None,mask=False,mask_file=None,return_type='alm'):
	import healpy as hp
	klm_file_root = dirc + 'baseline_MV/'
	alm_planck_input = hp.fitsfunc.read_alm(klm_file_root + 'dat_klm.fits')

	ls,_ = hp.Alm.getlm(lmax=lmax)
	alm_planck_input = hp.almxfl(alm_planck_input, (ls > lmin).astype(int))

	lmax_planck = hp.Alm.getlmax(len(alm_planck_input))
	ls_new, ms_new = hp.Alm.getlm(lmax=lmax)
	idxs_old = hp.Alm.getidx(lmax=lmax_planck, l=ls_new, m=ms_new)
	idxs_new = hp.Alm.getidx(lmax=lmax, l=ls_new, m=ms_new)
	alm_planck = np.zeros(hp.Alm.getsize(lmax), dtype=complex)
	alm_planck[idxs_new] = alm_planck_input[idxs_old]

	if return_type == 'alm' and mask==False:
		return alm_planck

	elif return_type == 'alm' and mask==True:
		raise ValueError("Better to not mask a alm")

	elif return_type == 'map' and mask==True:
		if mask_file == None:
			raise ValueError("Please provide the mask file")
		else:
			map_planck = hp.sphtfunc.alm2map(alm_planck,nside)
			mask_oringnal = hp.read_map(mask_file)
			mask = hp.ud_grade(mask_oringnal, nside)
			n_unmasked = np.count_nonzero(mask)
			map_planck = hp.ma(map_planck)
			map_planck.mask = np.logical_not(mask)

			fsky = n_unmasked/len(map_planck)

			return map_planck,fsky
	elif return_type == 'map' and mask==False:
		map_planck = hp.sphtfunc.alm2map(alm_planck,nside)
		return map_planck


def get_noise_bias(lmax=None,lmin=0,dirc=None):
	klm_file_root = dirc + 'baseline_MV/'
	cl_planck_bias = np.loadtxt(klm_file_root + 'nlkk_bias.dat')

	ls_bias = cl_planck_bias[:,0]
	RD_N0 = cl_planck_bias[:,1] # same as the "noise" var
	MC_N0 = cl_planck_bias[:,2]
	N1 = cl_planck_bias[:,3]
	PS = cl_planck_bias[:,4]
	diff = cl_planck_bias[:,5]

	RD_N0 = RD_N0[:lmax+1]
	diff = diff[:lmax+1]
	N1 = N1[:lmax+1]
	PS = PS[:lmax+1]

	RD_N0[:lmin] = 0
	diff[:lmin] = 0
	N1[:lmin] = 0
	PS[:lmin] = 0

	all_bias_noise = RD_N0 + diff + N1 + PS

	return [ls_bias,all_bias_noise]

def get_Planck_fidutial(lmax=None,lmin=0,dirc=None):
	klm_file_root = dirc + 'baseline_MV/'
	cl_planck_signal_noise = np.loadtxt(klm_file_root + 'nlkk.dat')
	cl_planck_bias = np.loadtxt(klm_file_root + 'nlkk_bias.dat')

	ls_Planck = cl_planck_signal_noise[:,0]
	noise = cl_planck_signal_noise[:,1]
	noise_signal =cl_planck_signal_noise[:,2]

	noise_signal = noise_signal[:lmax+1]
	noise = noise[:lmax+1]
	noise_signal[:lmin] = 0
	noise[:lmin] = 0
	signal = noise_signal - noise
	ls_Planck = ls_Planck[:lmax+1]

	return [ls_Planck,signal,noise_signal]