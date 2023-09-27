import healpy as hp
import numpy as np
import mtneedlet as nd

def needlet_filter(nside, j, B=2,lmax =None,alm=None):
    bl = nd.standardneedlet(B, j, lmax)
    if alm is not None:
        filtered_alm = hp.almxfl(alm, bl)
        return bl,filtered_alm
    else:
        return bl