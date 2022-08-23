#!/usr/bin/env python
# coding: utf-8
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import glob

sys.path.append('/data/id13/inhouse12/jiliang/code_v3/xs_proc')
from visual_func import *
from xs_data_proc import *
from proc_data_ana import *


def load_proc_qphi_size(fn,proc_type="integrate2d"):
    # this return the scan shape saved in jiliang's customized "*_proc.h5"
    # fn is the file name of '*_proc.h5' data
    # kw is keyword of h5 dataset or attribute
    with h5py.File(fn,'r') as f:
        fc = f[proc_type]
        qphi_shape = fc['map_qphi'].shape
        return qphi_shape
    
def load_proc_single_Iq(fn,r,c,proc_type="integrate1d"):
    # this return the qphi data saved in jiliang's customized "*_proc.h5"
    # fn is the file name of '*_proc.h5' data
    # kw is keyword of h5 dataset or attribute
    with h5py.File(fn,'r') as f:
        fc = f[proc_type]
        Iq = fc['map_Iq'][r,c]
        return Iq

def load_proc_single_qphi(fn,r,c,proc_type="integrate2d"):
    # this return the qphi data saved in jiliang's customized "*_proc.h5"
    # fn is the file name of '*_proc.h5' data
    # kw is keyword of h5 dataset or attribute
    with h5py.File(fn,'r') as f:
        fc = f[proc_type]
        qphi = fc['map_qphi'][r,c]
        return qphi
    
def single_qphi_sum_roi(r,c,fn=None,a=None,q=None,qmin=None,
                        qmax=None,amin=None,amax=None,mask=None,
                        bkgd=None,vs_axis='q'):
    '''
    this returns the average 1D intensity profile of roi for desired 
    q and azimuth range.
    
    Parameters
    ___________
    r : row position of scan
    c : column position of scan
    fn : the file name of correlated processed h5 file, which is
         usualy with subfix "*_proc.h5" and contains qphi map of scan
    a : angle coordinate of qphi map
    q : q coordinate of qphi map
    qmin : minimum of q for roi
    qmax : maximum of q for roi
    amin : minimum of a for roi
    amax : maximum of a for roi
    mask : mask is boolean array, which let the masked qphi region having value NaN
    bkgd : background array, which is subtracted from origin qphi
    vs_axis : determine the averaging direction of roi. vs_axis = 'q' will average along
              azimuthal axis, vs_axis = 'a' will average along the q axis.
    Returns
    ___________
    coordinate,intensity_profile
    coordinate: if vs_axis = 'q', return the q coordinate, else return a coordinate
    intensity_profile: average intensity profile
    '''
    qphi = np.copy(load_proc_single_qphi(fn,r,c,proc_type="integrate2d"))
    if not isinstance(mask,type(None)):
        qphi[mask] = np.nan
    if not isinstance(bkgd,type(None)):
        qphi -= bkgd
    if isinstance(amin,type(None)):
        aid1 = 0
    else:
        aid1 = np.argmin(np.abs(a-amin))
    if isinstance(amax,type(None)):
        aid2 = len(a)
    else:
        # this is because the numpy array index start from 0, thus, end index should add 1 to have correct length
        aid2 = np.argmin(np.abs(a-amax))+1
    if isinstance(qmin,type(None)):
        qid1 = 0
    else:
        qid1 = np.argmin(np.abs(q-qmin))
    if isinstance(qmax,type(None)):
        qid2 = len(q)
    else:
        # this is because the numpy array index start from 0, thus, end index should add 1 to have correct length
        qid2 = np.argmin(np.abs(q-qmax))+1
    if aid1 > aid2:
        print('amax must not smaller amin')
        return
    elif aid1 == aid2:
        aid2 += 1
    if qid1 > qid2:
        print('qmax must not smaller qmin')
        return
    elif qid1 == qid2:
        qid2 += 1
    if vs_axis == 'q':
        return q[qid1:qid2],np.nanmean(qphi[aid1:aid2,qid1:qid2],axis=0)
    elif vs_axis == 'a':
        return a[aid1:aid2],np.nanmean(qphi[aid1:aid2,qid1:qid2],axis=1)
    
def gaussian(x,am,mu,sigma):
    I = am*np.exp(-(x-mu)**2/(2*sigma**2))#/(sigma*(2*np.pi)**.5)
    return I

def lorentz(x,am,mu,sigma):
    I = am/(np.pi*sigma*(1+(x-mu)**2/sigma**2))
    return I

def azi_distr(a,am,mu,sigma,azi_diff=180,bkgd=0,peak_shape='lorentz'):#'gaussian'):
    # a is angle coord for azimuth distribution, unit for a is degree
    # am,mu,sigma are required keyword inputs
    if peak_shape == 'gaussian':
        I = gaussian(a,am,mu,sigma)+gaussian(a,am,mu+azi_diff,sigma)+\
            gaussian(a,am,mu+360,sigma)+gaussian(a,am,mu+azi_diff+360,sigma)+\
            gaussian(a,am,mu-360,sigma)+gaussian(a,am,mu+azi_diff-360,sigma)+bkgd
    elif peak_shape == 'lorentz':
        I = lorentz(a,am,mu,sigma)+lorentz(a,am,mu+azi_diff,sigma)+\
            lorentz(a,am,mu+360,sigma)+lorentz(a,am,mu+azi_diff+360,sigma)+\
            lorentz(a,am,mu-360,sigma)+lorentz(a,am,mu+azi_diff-360,sigma)+bkgd
    return np.abs(I)

def amax_intensity_profile(qphi,angle,qw,azi_max,da=2,sym_corr=True):
    q,Ia = sum_roi(qphi,angle,qw,amin=azi_max-da,amax=azi_max+da,vs_axis='q')
    if sym_corr:
        if azi_max < 0:
            q,Ia1 = sum_roi(qphi,angle,qw,amin=azi_max+180-da,amax=azi_max+180+da,vs_axis='q') 
        else:
            q,Ia1 = sum_roi(qphi,angle,qw,amin=azi_max-180-da,amax=azi_max-180+da,vs_axis='q')
        return q,np.nanmean(np.vstack((Ia,Ia1)),axis=0)
    else:
        return q,Ia

from lmfit import Model
from scipy.signal import medfilt

def azi_fit_diff(angle,Ia,mu_initial=None):
    '''
    this returns the fitting of azimuth intensity profie    
    Parameters
    ___________
    angle : angle coordinate of intensity distribution
    Ia    : azimuth intensity profile for fitting
    mu_initial : initial guess of peak position for fitting. If initial not 
                 properly chosen, the fitting may fail due to trap in the local minimum
    
    Returns
    ___________
    mu : maximum of gaussian of fitting
    sigma : standarad devation of gaussian of fittig
    am : amplitude of gaussian of fitting
    bkgd : constant background for fitting
    azi_diff : angle difference between two gaussain of fitting, should be around 180 degrees
    '''
    amod = Model(azi_distr)
    
    if isinstance(mu_initial,type(None)):
        Ia_medfilt = medfilt(Ia[np.isnan(Ia)==0])
        angle_medfilt = angle[np.isnan(Ia)==0]
        mu_initial = angle_medfilt[np.nanargmax(Ia_medfilt)]
    else:
        mu_initial = mu_initial
    amod.set_param_hint('am',value=np.nanmax(medfilt(Ia[np.isnan(Ia)==0])),min=0)
    amod.set_param_hint('mu',value=mu_initial,min=-180,max=180)
    amod.set_param_hint('sigma',value=40,min=1,max=90)
    amod.set_param_hint('azi_diff',value=180,min=160,max=200)
    amod.set_param_hint('bkgd',value=0.1,min=0)
    params = amod.make_params()
    res   = amod.fit(Ia[np.isnan(Ia)==0],params,a=angle[np.isnan(Ia)==0])
    mu    = res.params['mu'].value
    sigma = res.params['sigma'].value
    am    = res.params['am'].value
    bkgd  = res.params['bkgd'].value
    azi_diff = res.params['azi_diff'].value
    return mu,sigma,am,bkgd,azi_diff

def equator_intensity_profile(qphi,angle,qw,mu):
    qphi = np.copy(qphi)   
    q,I = amax_intensity_profile(qphi,angle,qw,mu)
    return I

def single_qphi_fit_paral(num,fn,scan_size,qmin,qmax,int_threshold=1,
                          std_threshold=2,mu_initial=90,mask=None,
                          bkgd_pos=None):
    _  = int(num/scan_size[1])
    __ = int(num%scan_size[1])
    if isinstance(bkgd_pos,type(None)):
        bkgd_pos = int(0)
    # here is a problem only consider the certain column of each row as correlated background of that row
    #
    bkgd = load_proc_single_qphi(fn,_,bkgd_pos)
    a1,I1 = single_qphi_sum_roi(_,__,fn,a=azi,q=q,qmin=qmin,qmax=qmax,vs_axis='a',bkgd=bkgd,mask=mask)
        
    if np.max(medfilt(I1[np.isnan(I1)==0])) < int_threshold:
        #print(1,_,__)
        mu = np.nan
        sigma = np.nan
        azi_diff = np.nan
    else:
        mu,sigma,am,bkgd,azi_diff = azi_fit_diff(a1,I1,mu_initial=mu_initial)
        if am > np.nanstd(I1)*std_threshold:
            mu = mu
            sigma = sigma
            azi_diff = azi_diff
        else:
            #print(2,_,__)
            #print(am,np.nanstd(I1)*std_threshold)
            mu = np.nan
            sigma = np.nan
            azi_diff = np.nan
    return mu,sigma,azi_diff 

from multiprocessing import Pool

def single_qphi_ori_determ(num,fn,q,a,qmin,qmax,scan_shape,
                          ll_thrhd=None,hl_thrhd=None,bkgd=None):
    total_num = scan_shape[0]*scan_shape[1]
    r = int(num/scan_shape[1])
    c = int(num%scan_shape[1])
    with h5py.File(fn,'r') as f:
        fc = f['integrate2d']
        qphi = fc['map_qphi'][r,c,:]
    if not isinstance(ll_thrhd,type(None)):
        qphi[qphi<=ll_thrhd] = np.nan
    if not isinstance(hl_thrhd,type(None)):
        qphi[qphi>=hl_thrhd] = np.nan
    if not isinstance(bkgd,type(None)):
        qphi = qphi - bkgd
    a,I = sum_roi(qphi,a,q,qmin,qmax,vs_axis='a')
    amax,awid = ori_determ_max(a,I)
    return amax,awid

def ori_determ2d_para(
                 fn,
                 azi,
                 q,
                 #qid,
                 qmin,
                 qmax,
                 mask= None,
                 ll_thrhd = None,
                 hl_thrhd = None,
                 num_cores=16,
                 bkgd=None
                ):
    scan_shape = load_proc_qphi_size(fn)
    #qphi = np.copy(qphi)
    qid1 = np.argmin(np.abs(q-qmin))
    qid2 = np.argmin(np.abs(q-qmax))+1

    ori_mat = np.zeros((scan_shape[:2]))*np.nan
    wid_mat = np.zeros((scan_shape[:2]))*np.nan
    res = parallel_func(single_qphi_ori_determ,num_cores,
                    np.arange(scan_shape[0]*scan_shape[1]),fn=fn,
                    a=azi,q=q,qmin=qmin,qmax=qmax,scan_shape=scan_shape,
                    ll_thrhd=ll_thrhd,hl_thrhd=hl_thrhd,bkgd=bkgd)
    for _ in range(scan_shape[0]):
        for __ in range(scan_shape[1]):
            if isinstance(mask,type(None)):
                (ori_mat[_,__],wid_mat[_,__]) = res[int(_*scan_shape[1]+__)]
            else:
                if mask[_,__]:
                    (ori_mat[_,__],wid_mat[_,__]) = res[int(_*scan_shape[1]+__)]
    return ori_mat,wid_mat


def ori_determin(fn,qmin,qmax,q,azi,
                 int_threshold=1,std_threshold=2,
                 mu_initial=90,num_cores=16,mask=None,bkgd_pos=125):
    scan_size = load_proc_qphi_size(fn)
    ori_mat  = np.zeros(scan_size[:2])
    wid_mat  = np.zeros(scan_size[:2])
    diff_mat = np.zeros(scan_size[:2])
    
    res = parallel_func(single_qphi_fit_paral,num_cores,np.arange(scan_size[0]*scan_size[1]),
                        fn=fn,scan_size=scan_size,qmin=qmin,qmax=qmax,
                        int_threshold=int_threshold,std_threshold=std_threshold,
                        mu_initial=mu_initial,mask=mask,bkgd_pos=bkgd_pos)
    
    for _ in range(scan_size[0]):
        for __ in range(scan_size[1]):
            ori_mat[_,__]  = res[int(_*scan_size[1]+__)][0]
            wid_mat[_,__]  = res[int(_*scan_size[1]+__)][1]
            diff_mat[_,__] = res[int(_*scan_size[1]+__)][2]
    return ori_mat,wid_mat,diff_mat


def plot_quiver(ori_mat,wid_mat,
                scale=0.2,clim=None,
               ):
    # ori_mat is anlge with unit of degrees
    # wid_mat is similar to fwhm, inverse to sharpness
    # length_adjust, try to correlate the length to wid_mat
    ori = np.radians(ori_mat)
    if isinstance(clim,type(None)):
        clim = (np.nanmin(1/wid_mat),np.nanmax(1/wid_mat))
    else:
        clim = clim
    wid_mat = 1/wid_mat
    #wid_mat /= np.nanmean(wid_mat)
    #wid_std  = np.nanstd(wid_mat)
    
    #wid_mean = np.mean(wid_mat)
    #wid_mat[wid_mat<=(wid_mean-wid_std*3)] = wid_mean-wid_std*3
    #wid_mat[wid_mat>=(wid_mean+wid_std*3)] = wid_mean+wid_std*3
    (x,y) = np.meshgrid(np.arange(ori.shape[1]),
                        np.arange(ori.shape[0]))
    (u,v) = (np.cos(ori),np.sin(ori))#*wid_mat
    plt.figure()#100,figsize=(12,9))
    ax = plt.gca()
    ax.quiver(x,y,u,v,wid_mat,
              units='xy',angles='uv',clim=clim,
              scale=scale,headwidth=0,pivot='mid',cmap='hsv')
    ax.set_ylim(ax.get_ylim()[::-1])
    
    
# here is degree for vals is -90 to 90, orientation had been adjusted within this range
def circle_sine_ramp(r_max=20, r_min=0, amp=np.pi/5, cycles=50,
                     power=2, nr=50, ntheta=1025):
    r, t = np.mgrid[r_min:r_max:nr*1j, 0:2*np.pi:ntheta*1j]
    r_norm = (r - r_min)/(r_max - r_min)
    vals = amp * r_norm**power * np.sin(cycles*t) + t
    vals = np.mod(vals, 2*np.pi) - np.pi
    vals[vals<0] += np.pi
    vals -= np.pi/2
    return t, r, vals

def cylinder_form_factor(qr,r,qz=0,height=1,density=1):
    V = 2*np.pi*r**2*height
    F = 2*j0(qz*height/2)*j1(qr*r)/qr/r + 1j*0
    F *= density*V
    return np.abs(F)

def center_center_interference(q,d):
    return np.abs(2*(1+j0(q*d)))

def norm_func(x):
    #guassian function with mean = 0 and sigma =1
    return np.exp(-x**2/2)/np.sqrt(2*np.pi)

def poly_diversity(func,q,r,dr,dsigma=1,nsigma=3):
    '''
    for poly_diversity, pdf using the norm_func with sigma = 1
    sample points depends on both nsigma (largest sigma reached) and 
    the sampling span within each sigma region
    Parameter
    ______________
    func : function for polydiversity calculation, function could be common
           form factor, such as sphere, cyinder or cubic function or structure factor
    q : q coordinates
    r : radius of sphere or cylinder, or width of cubic
    dr : variation of radius for polydiversity calculation
    dsigma : sampling rate of polydiversity, 1 means that calculation samples at every standard
             deviation postion, 0.5 means that every standrad devation range is equally sampled twice.
    nsigma : how much standard deviations will be covered for calculation.
    Return
    _______________
    I : form factor or structure factor including polydiversity
    '''
    for i in range(int(nsigma/dsigma)*2+1):
        p = i - int(nsigma/dsigma)
        if i == 0:
            I = func(q,r+dr*dsigma*p)*norm_func(p*dsigma)
        else:
            I += func(q,r+dr*dsigma*p)*norm_func(p*dsigma)
    return I

# define function for cylinderical volume normalization

def volume_normalization(r,x,dx=0.1,xc=None,hist_interval=2.5):
    # the choice of hist_interval and x should be compatible with each, othewise hist
    # will be distorted due to uneven distribution of interval
    if isinstance(xc,type(None)):
        xc = x/2
    else:
        xc = xc
    xx,yy = np.meshgrid(np.arange(0,x,dx),np.arange(0,x,dx))
    xx_new = xx-xc
    # here should notice that center of fiber could displace to the horiztontal center 
    # but not affect the beam direction, which is the projection direction
    yy_new = yy-(x/2)
    xr = np.sqrt(xx_new**2+yy_new**2)
    w  = (xr<=r).astype(np.int32)
    I,dx = np.histogram(xx,bins=np.arange(0,x+hist_interval,hist_interval),weights=w)
    return I/np.min(I),dx,w

# function need to improve

def plot_spot(x,y,vmin=0,vmax=30,log=False,
              qphi_plot=True,pttn_plot=False,line_plot=True,
              qphi=None,path_list=None,
              qmin=1.83,qmax=1.89,amin=None,amax=None,mask=None):
    if pttn_plot:
        if not isinstance(path_list,type(None)):
            plt.subplots()
            plot_pttn(path_list[path_idx[y,x]],
                      pttn_idx[y,x],log=log,
                      vmin=vmin,vmax=vmax,mask=pttn_mask,cmap='jet')
        else:
            print("path list is missing")
            return
    if qphi_plot:
        if not isinstance(qphi,type(None)):
            plt.subplots()
            plot_qphi(qphi[y,x],q=q,a=azi,vmin=vmin,
                      vmax=vmax,log=log,cmap='jet')
            if isinstance(amin,type(None)): amin = np.min(azi)
            if isinstance(amax,type(None)): amax = np.max(azi)
            plt.vlines(qmin,amin,amax,linewidth=1.5)
            plt.vlines(qmax,amin,amax,linewidth=1.5)
        else:
            print("qphi is missing")
            return
    if line_plot:
        a,Ia = sum_roi(qphi[y,x],azi,q,qmin=qmin,qmax=qmax,mask=mask,vs_axis='a')
        plt.subplots()
        plt.plot(a,Ia)
        plt.title(r'$reflection\,\,(0\,\,1\,\,3)')


def asymetry_compare(x,y,qmin=1.83,qmax=1.89):
    # compare the intensity difference of maximum at (-180,0) and (0,180)
    a1,Ia1 = sum_roi(qphi[y,x],azi,q,qmin=1.8,qmax=1.9,mask=mask,
                 amin=-180,amax=0,vs_axis='a')
    i1 = np.nanmax(Ia1)
    a2,Ia2 = sum_roi(qphi[y,x],azi,q,qmin=1.8,qmax=1.9,mask=mask,
                 amin=0,amax=180,vs_axis='a')
    i2 = np.nanmax(Ia2)
    if i1 >= i2:
        sign = 1.
    else:
        sign = -1.
    
    return sign*(np.max([i1,i2])/np.min([i1,i2])-1)

def asym_mat(qphi,qmin,qmax):
    mat = np.zeros(qphi.shape[:2])
    for y in range(qphi.shape[0]):
        for x in range(qphi.shape[1]):
            mat[y,x] = asymetry_compare(x,y,qmin,qmax)
    return mat


#----------------------
#parallel roi sum of qphi avoid loading the whole qphi data


def single_qphi_roi_sum(r,c,fn=None,
                        a=None,q=None,
                        qmin=None,qmax=None,
                        amin=None,amax=None,mask=None,
                        low_thrd=0,high_thrd=np.inf,
                        bkgd=None):
    '''
    this returns the mean intensity of q and angle roi of correlated qphi 
    q and azimuth range.
    
    Parameters
    ___________
    r : row position of scan
    c : column position of scan
    fn : the file name of correlated processed h5 file, which is
         usualy with subfix "*_proc.h5" and contains qphi map of scan
    a : angle coordinate of qphi map
    q : q coordinate of qphi map
    qmin : minimum of q for roi
    qmax : maximum of q for roi
    amin : minimum of a for roi
    amax : maximum of a for roi
    mask : mask is boolean array, which let the masked qphi region having value NaN
    bkgd : background array, which is subtracted from origin qphi
    Returns
    ___________
    mean intensity of roi of qphi
    '''
    qphi = np.copy(load_proc_single_qphi(fn,r,c,proc_type="integrate2d"))
    qphi[qphi<=low_thrd] = np.nan
    qphi[qphi>=high_thrd] = np.nan
    if not isinstance(mask,type(None)):
        qphi[mask] = np.nan
    if not isinstance(bkgd,type(None)):
        qphi -= bkgd
    if isinstance(amin,type(None)):
        aid1 = 0
    else:
        aid1 = np.argmin(np.abs(a-amin))
    if isinstance(amax,type(None)):
        aid2 = len(a)
    else:
        # this is because the numpy array index start from 0, thus, end index should add 1 to have correct length
        aid2 = np.argmin(np.abs(a-amax))+1
    if isinstance(qmin,type(None)):
        qid1 = 0
    else:
        qid1 = np.argmin(np.abs(q-qmin))
    if isinstance(qmax,type(None)):
        qid2 = len(q)
    else:
        # this is because the numpy array index start from 0, thus, end index should add 1 to have correct length
        qid2 = np.argmin(np.abs(q-qmax))+1
    if aid1 > aid2:
        print('amax must not smaller amin')
        return
    elif aid1 == aid2:
        aid2 += 1
    if qid1 > qid2:
        print('qmax must not smaller qmin')
        return
    elif qid1 == qid2:
        qid2 += 1
    return np.nanmean(np.nanmean(qphi[aid1:aid2,qid1:qid2],axis=0),axis=0)


def single_qphi_roi_sum_paral(num,fn,scan_size,
                              q,azi,
                              qmin,qmax,amin,amax,
                              low_thrd=0,high_thrd=np.inf,
                              mask=None,bkgd_pos=None,
                              bkgd_direction='c'):
    _  = int(num/scan_size[1])
    __ = int(num%scan_size[1])
    if not isinstance(bkgd_pos,type(None)):
        bkgd_pos = int(bkgd_pos)
        if bkgd_direction == 'c':
            # bkgd is chosen at certain column
            bkgd = load_proc_single_qphi(fn,_,bkgd_pos)
        elif bkgd_direction == 'r':
            bkgd = load_proc_single_qphi(fn,bkgd_pos,__)
    else:
        bkgd = None
    I = single_qphi_roi_sum(_,__,fn,a=azi,q=q,qmin=qmin,qmax=qmax,bkgd=bkgd,mask=mask)
    return I        

from multiprocessing import Pool

def qphi_roi_sum_2dmap(fn,qmin,qmax,amin,amax,
                       q,azi,num_cores=36,
                       low_thrd=0,high_thrd=np.inf,
                       mask=None,bkgd_pos=None):
    scan_size = load_proc_qphi_size(fn)
    qphi_roi_map  = np.zeros(scan_size[:2])
    res = parallel_func(single_qphi_roi_sum_paral,
                        num_cores,
                        np.arange(scan_size[0]*scan_size[1]),
                        fn=fn,scan_size=scan_size,
                        q=q,azi=azi,
                        qmin=qmin,qmax=qmax,
                        amin=amin,amax=amax,
                        low_thrd=low_thrd,high_thrd=high_thrd,
                        mask=mask,bkgd_pos=None)
    for _ in range(scan_size[0]):
        for __ in range(scan_size[1]):
            qphi_roi_map[_,__]  = res[int(_*scan_size[1]+__)]#[0]
    return qphi_roi_map
