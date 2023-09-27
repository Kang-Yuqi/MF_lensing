import sympy as sp
import numpy as np

def Hpol(order):
    x = sp.Symbol('x')
    der = sp.exp(-x**2/2)
    for i in range(order):
        der = sp.diff(der, x)
    der = der*(-1)**order*sp.exp(x**2/2)
    f = sp.lambdify(x, sp.simplify(der), 'numpy')
    return f

def MF_correction(range_std,threshold_num,sigma_sk=None,order=2):

    threshold = np.linspace(-range_std, range_std, num=threshold_num)

    [sigma,sigma1,S,S1,S2,K,K1,K12,K22] = sigma_sk

    x = threshold / sigma

    Pre0 = np.sqrt(2*np.pi)**(-1)*np.exp(-x**2/2)
    Pre1 = (1/8)*sigma1/(np.sqrt(2)*sigma)*np.exp(-x**2/2)
    Pre2 = 1/(2*np.pi)**(3/2)*sigma1**2/(2*sigma**2)*np.exp(-x**2/2)

    # first order 
    v0_first = (1/6)*S*Hpol(2)(x)
    v1_first = (1/6)*S*Hpol(3)(x) + (1/3)*S1*Hpol(1)(x)
    v2_first = (1/6)*S*Hpol(4)(x) + (2/3)*S1*Hpol(2)(x) + (1/3)*S2*Hpol(0)(x)

    if order == 1:
        MF_corr = [Pre0*v0_first*sigma, Pre1*v1_first*sigma, Pre2*v2_first*sigma]
    elif order == 2:
        # second order
        v0_second = (1/72)*S**2*Hpol(5)(x) + (1/24)*K*Hpol(3)(x)

        v1_second = (1/72)*S**2*Hpol(6)(x) + ((1/24)*K + (1/18)*S*S1)*Hpol(4)(x) + ((1/8)*K1-(1/18)*S1**2)*Hpol(2)(x) + ((-1/16)*K12+(1/16)*K22)*Hpol(0)(x)

        v2_second = (1/72)*S**2*Hpol(7)(x) + ((1/24)*K + (2/18)*S*S1)*Hpol(5)(x) + 2*((1/8)*K1+(1/36)*S*S2)*Hpol(3)(x) + 2*((2/16)*K22-(1/9)*S1*S2)*Hpol(1)(x)

        MF_corr = [Pre0*(v0_first*sigma + v0_second*sigma**2), 
                    Pre1*(v1_first*sigma + v1_second*sigma**2), 
                    Pre2*(v2_first*sigma + v2_second*sigma**2)]
    else:
        raise ValueError("order can only be 1 or 2")

    return MF_corr


