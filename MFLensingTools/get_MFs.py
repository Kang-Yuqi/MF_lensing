import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import sys

def first_dp(alm,lmax):
    
    alm_result = np.ones(len(alm))
    alm_result_dp = alm_result*(0+0j)
    ls = np.arange(0, lmax)
    for l in ls:
        ms = np.arange(0, l+1)
        idx = hp.Alm.getidx(lmax, l, ms)   
        alm_result_dp[idx] = 1j*ms*alm[idx]
    return alm_result_dp

def second_dp(alm,lmax):
    
    alm_result = np.ones(len(alm))
    alm_result = alm_result*(0+0j)
    for l in range(0,lmax+1): 
        ms = np.arange(0, l+1)
        idx = hp.Alm.getidx(lmax,l,ms)
        alm_result[idx] = -ms**2*alm[idx]
    return alm_result

def MFs_map(nside,input_map,range_std,threshold_num,lmax,width_divide=20,mask_file=None):

	threshold = np.linspace(-range_std, range_std, num=threshold_num)

	if len(input_map)==hp.pixelfunc.nside2npix(nside):
		alm = hp.sphtfunc.map2alm(input_map, lmax=lmax, iter=5)
	elif len(input_map)==hp.Alm.getsize(lmax, mmax=None):
		alm = input_map
	else:
		raise ValueError("length of map/alm is not consistance with nside/lmax")

	field,map_dt,_ = hp.sphtfunc.alm2map_der1(alm,nside)
	alm_dt = hp.sphtfunc.map2alm(map_dt, lmax=lmax, iter=5)
	_,map_dtt,_ = hp.sphtfunc.alm2map_der1(alm_dt,nside)

	alm_dp = first_dp(alm,lmax)
	alm_dpp = second_dp(alm,lmax)
	alm_dtp = first_dp(alm_dt,lmax)

	ipix = np.arange(12*(nside**2)) # number of pixels for a nside number
	(theta, phi)=hp.pixelfunc.pix2ang(nside, ipix, nest=False, lonlat=False) # define the polar coordinated

	map_dp = hp.sphtfunc.alm2map(alm_dp,nside)/np.sin(theta)
	map_dpp = hp.sphtfunc.alm2map(alm_dpp,nside)/(np.sin(theta)**2)
	map_dtp = hp.sphtfunc.alm2map(alm_dtp,nside)/np.sin(theta)

	if mask_file is not None:
		mask_oringnal = hp.read_map(mask_file).astype(np.bool_)
		mask = hp.ud_grade(mask_oringnal, nside)
		n_unmasked = np.count_nonzero(mask)
		field = hp.ma(field)
		field.mask = np.logical_not(mask)
		map_dt = hp.ma(map_dt)
		map_dt.mask = np.logical_not(mask)
		map_dp = hp.ma(map_dp)
		map_dp.mask = np.logical_not(mask)
		map_dtt = hp.ma(map_dtt)
		map_dtt.mask = np.logical_not(mask)
		map_dpp = hp.ma(map_dpp)
		map_dpp.mask = np.logical_not(mask)
		map_dtp = hp.ma(map_dtp)
		map_dtp.mask = np.logical_not(mask)
	else:
		n_unmasked = hp.nside2npix(nside)

	# print(n_unmasked/hp.nside2npix(nside)) # f_sky

	deltanu = np.std(field) / width_divide
	index_thresholds = [(field > (t - (deltanu/2.))) & (field < t + (deltanu/2.)) for t in threshold]

	V0_list = []
	for n in range(len(threshold)):
		if mask_file is not None: 
			index_threshold = np.where(field.filled(-1e30) > threshold[n])
		else:
			index_threshold = np.where(field > threshold[n])
			
		V0 = len(index_threshold[0])/n_unmasked # normalisation
		V0_list.append(V0)

	V1_list = []
	V2_list = []

	for n, thresh in enumerate(index_thresholds):
		index_threshold = np.where(thresh)[0]

		if mask_file is not None: 
			dt = map_dt[index_threshold].filled(0)
			dp = map_dp[index_threshold].filled(0)
		else:
			dt = map_dt[index_threshold]
			dp = map_dp[index_threshold]
		dtp = map_dtp[index_threshold]
		dtt = map_dtt[index_threshold]
		dpp = map_dpp[index_threshold]

		counter1 = np.sum(np.sqrt(dt**2 + dp**2)) / (4*deltanu)
		counter2 = np.sum((2*dt*dp*dtp - dt**2*dpp - dp**2*dtt) / (2*np.pi*(dt**2 + dp**2))) / deltanu
		V1_list.append(counter1 /n_unmasked)
		V2_list.append(counter2 /n_unmasked)

	return [V0_list, V1_list, V2_list]