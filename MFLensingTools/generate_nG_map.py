import numpy as np

def generate_nG_map(map_G,fnl=0,gnl=0):
	## add fnl (first order)
	map_nG_fnl = fnl*(map_G**2-np.mean(map_G**2))
	## add gnl (second order)
	map_nG_gnl = gnl*(map_G**3)

	map_nG = map_G + map_nG_fnl + map_nG_gnl # result nG map

	return map_nG