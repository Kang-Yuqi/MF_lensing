import numpy as np
import camb
import matplotlib.pyplot as plt

def Cl_camb(As,w,nu,lmax,lmin=0,H0=67,ombh2=0.02216,tau=0.06,nnu=3.044,omk=0,ns=0.96,lens_potential_accuracy=8):
	pars = camb.CAMBparams()
	pars.set_cosmology(H0=H0,
						ombh2=ombh2,
						omch2=w,
						tau=tau,
						mnu=nu,
						nnu=nnu,
						omk=omk)
	pars.InitPower.set_params(As=As*1e-9,ns=ns)

	pars.set_for_lmax(lmax,lens_potential_accuracy=lens_potential_accuracy)
	pars.NonLinearModel.set_params(halofit_version='mead')
	results = camb.get_results(pars)
	cl_camb=results.get_lens_potential_cls(lmax)
	cl_phi = cl_camb[:,0]
	cl_phi[:lmin]=0
	ls = np.linspace(0,lmax,lmax+1)

	return ls,cl_phi
