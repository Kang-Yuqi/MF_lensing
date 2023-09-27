import numpy as np
from scipy.special import erfc

def analytical_MFs(cl,lmin,lmax,range_std,threshold_num):

    threshold = np.linspace(-range_std, range_std, num=threshold_num)

    ls = np.arange(lmin, lmax+1)
    sigma = np.sqrt(np.sum(((2*ls+1)*cl[ls])/(4*np.pi)))
    sigma1 = np.sqrt(np.sum(ls*(ls+1)*(2*ls+1)*cl[ls]/((4*np.pi))))

    # Gaussian MFs
    x = threshold / sigma
    v0 = (1./2.)*erfc(x/np.sqrt(2))
    v1 = (1./8.)*sigma1/(np.sqrt(2)*sigma)*np.exp(-x**2/2)
    v2 = 1/((2*np.pi)**(3./2.))*sigma1**2/(2*sigma**2)*np.exp(-x**2/2)*x

    # print('analytical sigma sigma1',sigma,sigma1)
    return [v0,v1,v2]


