import numpy as np
from scipy.special import j0,j1

'''
here are the calculations of some basic geometry shape and peak functions
the calculation corresponds to the 1D intensity profile
'''

def spherical_form_factor(q,r,density=1):
    V = 4./3*np.pi*r**3
    return density*((3*V*(np.sin(q*r)-q*r*np.cos(q*r))/((q*r)**3)))**2/V

def box_form_factor(q,width,density=1):
    return density*np.sinc(q*width)

def cylinder_form_factor(qr,r,qz=0,height=1,density=1):
    V = 2*np.pi*r**2*height
    F = 2*j0(qz*height/2)*j1(qr*r)/qr/r
    F *= density*V
    return F

def center_center_interference(q,d):
    return 2*(1+j0(q*d))

def norm_func(x):
    #guassian function with mean = 0 and sigma =1
    return np.exp(-x**2/2)/np.sqrt(2*np.pi)
    
def poly_dispersity(func,q,r,dr,dsigma=1,nsigma=3):
    # for poly_diversity, pdf using the norm_func with sigma = 1
    # sample points depends on both nsigma (largest sigma reached) and 
    # the sampling span within each sigma region
    for i in range(int(nsigma/dsigma)*2+1):
        p = i - int(nsigma/dsigma)
        if i == 0:
            I = func(q,r+dr*dsigma*p)*norm_func(p*dsigma)
        else:
            I += func(q,r+dr*dsigma*p)*norm_func(p*dsigma)
    return I
