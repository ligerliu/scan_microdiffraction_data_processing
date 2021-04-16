import numpy as np
from skimage.draw import polygon
import matplotlib.pyplot as plt
from functools import partial
import h5py
import fabio


def load_proc(fn,proc_type="integrate2d"):
    # load the data from proc hdf, specific formate for JL analysis
    # maybe change the input to allow load partially later
    with h5py.File(fn,"r") as f:
        fc = f[proc_type]
        azi       = np.array(fc['angle'])#[:]
        ct        = np.array(fc['beam_intensity'])#[:]
        qphi      = np.array(fc['map_qphi'])#[:]
        q         = np.array(fc['q'])#[:]
        # str could not be saved as dataset, thus save an attrs
        path_list = np.array(fc.attrs['origin_h5_path'])#[:]
        path_idx  = np.array(fc['path_idx'])#[:]
        pttn_idx  = np.array(fc['pttn_idx'])#[:]
    return azi,ct,qphi,q,path_list,path_idx,pttn_idx
    
def load_proc_dataset(fn,kw,proc_type="integrate2d"):
    # fn is the file name of '*_proc.h5' data
    # kw is keyword of h5 dataset or attribute
    with h5py.File(fn,'r') as f:
        fc = f[proc_type]
        if kw == 'origin_h5_path':
            return fc.attrs['origin_h5_path'][:]
        else:
            return np.array(fc[kw])
    
def sum_roi_2dmap(qphi,a,q,
            amin=None,
            amax=None,
            qmin=None,
            qmax=None,
            mask=None,
            ):
    # qphi is qphi_map from proc.h5 data has shape (y,x,phi,q)
    # y,x is scan position of microdiffraction. phi is azimuth, q has unit A-1
    # mask should have same shape to (qphi)
    if not isinstance(mask,type(None)):
        qphi[:,:,mask] = np.nan
    if isinstance(amin,type(None)):
        aid1 = 0
    else:
        aid1 = np.argmin(np.abs(a-amin))
    if isinstance(amax,type(None)):
        aid2 = len(a)
    else:
        aid2 = np.argmin(np.abs(a-amax))
    if isinstance(qmin,type(None)):
        qid1 = 0
    else:
        qid1 = np.argmin(np.abs(q-qmin))
    if isinstance(qmax,type(None)):
        qid2 = len(q)
    else:
        qid2 = np.argmin(np.abs(q-qmax))
    return np.nanmean(np.nanmean(qphi[:,:,aid1:aid2,qid1:qid2],axis=-1),axis=-1)


def sum_roi(qphi,a,q,
            amin=None,
            amax=None,
            qmin=None,
            qmax=None,
            mask=None,
            vs_axis='q',
            ):
    # qphi is 2d remesh of original diffraction pattern to phi,q
    # phi also called chi for azimuth meaning
    # q has unit A-1
    # mask should have same shape to (qphi)
    if not isinstance(mask,type(None)):
        qphi[mask] = np.nan
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


def ct_normalization(data,ct,scale=None):
    # here data is 2D, 3D or 4D array, for example the processed qphi data
    # ct is beam intensity, which had been saved in original or processed data
    data = np.copy(data)
    shape = data.shape
    dim   = np.ones(len(shape)).astype(int)
    if isinstance(scale,type(None)):
        scale = np.nanmean(ct)
    else:
        scale = scale
    if np.size(ct) == shape[0]:
        dim[0] = np.size(ct)
        ct = ct.reshape(dim)
        for _ in range(1,len(dim)):
            dim[_] = shape[_]
        dim[0] = 1
        ct = np.tile(ct,dim)
        data /= ct
        return data*scale
    elif np.size(ct) >= shape[0]*shape[1]:
        dim[0] = shape[0]
        dim[1] = shape[1]
        ct = ct[:(shape[0]*shape[1])]
        ct = ct.reshape(dim)
        for _ in range(2,len(dim)):
            dim[_] = shape[_]
        dim[0] = 1
        dim[1] = 1
        ct = np.tile(ct,dim)
        data /= ct
        return data*scale
    else:
        print("dimension of beam intensity is not consist with scan dimension of data")
        return

from scipy.signal import medfilt2d
def medfilt_pttn(data,kernel_size=3):
    # medfilt2d could directly ignore np.nan value
    # larger kernel size will leads to more smoother diffuse data
    smooth_data    = medfilt2d(data,kernel_size)
    high_frequence = data-smooth_data
    return smooth_data,high_frequence

def qphi_high_freq_separation(qphi,kernel_size):
    # this is slow should be parallelized to accelerate
    # also very ram consuming, this has no good solution yet
    diffuse = np.copy(qphi)*0
    mineral = np.copy(qphi)*0
    for i in range(qphi.shape[0]):
        for j in range(qphi.shape[1]):
            (diffuse[i,j],mineral[i,j]) = medfilt_pttn(qphi[i,j],kernel_size)
    return diffuse,mineral

def qphi_bkgd_sub(bkgd,qphi):
    # could be accelerated by parallelization
    data = np.copy(qphi)
    for i in range(qphi.shape[0]):
        for j in range(qphi.shape[1]):
            data[i,j] = data[i,j] - bkgd
    return data

#######################
def plot_qphi(qphi,
              vmin=None,
              vmax=None,
              log=True,
              q=None,
              a=None,
              mask=None,
              bkgd=None,
              **kwargs):
    if isinstance(q,type(None)):
        q = np.arange(qphi.shape[1])
    if isinstance(a,type(None)):
        a = np.arange(qphi.shape[0])
    if not isinstance(mask,type(None)):
        qphi[mask] = np.nan
    if not isinstance(bkgd,type(None)):
        qphi -= bkgd
    if log:
        qphi = np.log(qphi)
    if isinstance(vmin,type(None)):
        vmin = np.nanmean(qphi)-np.nanstd(qphi)*3
    if isinstance(vmax,type(None)):
        vmax = np.nanmean(qphi)+np.nanstd(qphi)*3
    plt.subplots()
    plt.imshow(qphi,vmin=vmin,vmax=vmax,
               extent=(q[0],q[-1],a[-1],a[0]),
               aspect=np.abs(q[-1]-q[0])/np.abs(a[-1]-a[0]),
               **kwargs)
    plt.ylabel(r'$\phi\,\,(^{o})$')
    plt.xlabel(r'$Q\,\,(\AA^{-1})$')

def plot_pttn(path,
              pttn,
              data_path='entry_0000/measurement/data',
              vmin = None,
              vmax = None,
              log  = True,
              mask = None,
              bkgd = None,
              **kwargs):
    with h5py.File(path,'r') as f:
        d = np.copy(f[data_path][pttn]).astype(float)
    if log:
        d = np.log(d)
    if isinstance(vmin,type(None)):
        vmin = np.nanmean(d)-np.nanstd(d)*3
    if isinstance(vmax,type(None)):
        vmax = np.nanmean(d)+np.nanstd(d)*3
    if not isinstance(mask,type(None)):
        d[mask] = np.nan
    if not isinstance(bkgd,type(None)):
        d -= bkgd
    plt.subplots()
    plt.imshow(d,vmin=vmin,vmax=vmax,
               **kwargs)

def pttn_of_int_map(int_map,
                    path_list,
                    path_idx,
                    pttn_idx,
                    data_path='entry_0000/measurement/data',
                    vmin1 = None,
                    vmax1 = None,
                    vmin2 = None,
                    vmax2 = None,
                    log  = True,
                    mask = None,
                    zoom_xlim=(0,-1),
                    zoom_ylim=(0,-1),
                    bkgd = None,
                    **kwargs):
    zoom_xlim = np.array(zoom_xlim)
    zoom_ylim = np.array(zoom_ylim)
    fig1,ax1 = plt.subplots()
    if isinstance(vmin1,type(None)):
        vmin1 = np.nanmean(int_map)-3*np.nanstd(int_map)
    if isinstance(vmax1,type(None)):
        vmax1 = np.nanmean(int_map)+3*np.nanstd(int_map)
    ax1.imshow(int_map,vmin=vmin1,vmax=vmax1)
    def onclick(event,
                data_path='entry_0000/measurement/data',
                vmin = vmin2,
                vmax = vmax2,
                log  = log,
                mask = mask,
                bkgd = bkgd,
                **kwargs):
        fig2,ax2 = plt.subplots()
        col = np.round(event.xdata).astype(int)
        row = np.round(event.ydata).astype(int)
        path = path_list[path_idx[row,col]]
        pttn = pttn_idx[row,col]
        #print(path,pttn)
        with h5py.File(path,'r') as f:
            d = np.copy(f[data_path][pttn]).astype(float)
        if not isinstance(bkgd,type(None)):
            d -= bkgd
        if log:
            d = np.log(d)
        if isinstance(vmin,type(None)):
            vmin = np.nanmean(d)-np.nanstd(d)*3
        if isinstance(vmax,type(None)):
            vmax = np.nanmean(d)+np.nanstd(d)*3
        if not isinstance(mask,type(None)):
            d[mask] = np.nan
        ax2.imshow(d,vmin=vmin,vmax=vmax,
                   **kwargs)
        if zoom_xlim[1] == -1:
            zoom_xlim[1] = d.shape[1]
        if zoom_ylim[1] == -1:
            zoom_ylim[1] = d.shape[0]
        ax2.axis([zoom_xlim[0],zoom_xlim[1],zoom_ylim[0],zoom_ylim[1]])
    cid = fig1.canvas.mpl_connect('button_press_event', onclick)

def gaussian(x,am,mu,sigma):
    I = am/(sigma*(2*np.pi)**.5)*np.exp(-(x-mu)**2/(2*sigma**2))
    return I

def lorentz(x,am,mu,sigma):
    I = am/(np.pi*sigma*(1+(x-mu)**2/sigma**2))
    return I

def azi_distr(a,am,mu,sigma,azi_diff=180,bkgd=0):
    # a is angle coord for azimuth distribution, unit for a is degree
    # am,mu,sigma are required keyword inputs
    # include plus minus one periodicity of azimuth to enusre correctness of curve
    I = gaussian(a,am,mu,sigma)+gaussian(a,am,mu+azi_diff,sigma)+\
        gaussian(a,am,mu+360,sigma)+gaussian(a,am,mu+azi_diff+360,sigma)+\
        gaussian(a,am,mu-360,sigma)+gaussian(a,am,mu+azi_diff-360,sigma)+bkgd
    return I

def ori_determ_max(a,d):
    # this is function to determine the orientation by angle of maximum 
    # intensity position, 
    # the sharpness is determined by the sum(d>ave(d))/max(d)
    # a correlate to azimuth coordinate
    # d is azimuth distrion of intenisty
    # for return, aid had been limited in [0,180]
    if len(a) != len(d):
        print("length of intensity should be same to azimuth coord.")
    idx = np.argmax(d[np.isnan(d)==0])
    aid = a[np.isnan(d)==0][idx]
    if (aid>=-180) and (aid<0):
        aid += 180
    width = np.nansum(d[d>=np.nanmean(d)])/(np.nanmax(d)-np.nanmean(d))
    return aid,width

def ori_determ2d(qphi,
                 azi,
                 qid,
                 mask=None,
                ):
    ori_mat = np.zeros((qphi.shape[:2]))*np.nan
    wid_mat = np.zeros((qphi.shape[:2]))*np.nan
    for _ in range(qphi.shape[0]):
        for __ in range(qphi.shape[1]):
            if isinstance(mask,type(None)):
                ori_mat[_,__],wid_mat[_,__] = ori_determ_max(azi,
                                                     qphi[_,__,:,qid])
            else:
                if mask[_,__]:
                    ori_mat[_,__],wid_mat[_,__] = ori_determ_max(azi,
                                                         qphi[_,__,:,qid])
    return ori_mat,wid_mat

def plot_quiver(ori_mat,wid_mat,
                length_adjust=True,
                scale=0.2,
               ):
    # ori_mat is anlge with unit of degrees
    # wid_mat is similar to fwhm, inverse to sharpness
    # length_adjust, try to correlate the length to wid_mat
    ori = np.radians(ori_mat)
    wid_mat = 1/wid_mat
    wid_mat /= np.nanmean(wid_mat)
    wid_std  = np.std(wid_mat)
    wid_mean = np.mean(wid_mat)
    wid_mat[wid_mat<=(wid_mean-wid_std*3)] = wid_mean-wid_std*3
    wid_mat[wid_mat>=(wid_mean+wid_std*3)] = wid_mean+wid_std*3
    (x,y) = np.meshgrid(np.arange(ori_mat.shape[1]),
                        np.arange(ori_mat.shape[0]))
    (u,v) = (np.cos(ori),np.sin(ori))
    plt.figure()#100,figsize=(12,9))
    ax = plt.gca()
    if length_adjust:
        u *= wid_mat**2
        v *= wid_mat**2
    ax.quiver(x,y,u,v,wid_mat,
               units='xy',
               scale=scale,headwidth=0,pivot='mid',cmap='hsv')
    ax.set_ylim(ax.get_ylim()[::-1])


def plot_spot(x,y,vmin=0,vmax=30,log=False,
              qphi_plot=True,pttn_plot=False,line_plot=True,
              qmin=1.83,qmax=1.89,amin=None,amax=None,mask=None,
              bkgd=None,):
    #pttn_mask = np.load('msk_file/ndet_m640.npz')['mask']
    if pttn_plot:
        plt.subplots()
        plot_pttn(path_list[path_idx[y,x]],
                  pttn_idx[y,x],log=log,bkgd=bkgd,
                  vmin=vmin,vmax=vmax,mask=mask,cmap='jet')
    if qphi_plot:
        plt.subplots()
        plot_qphi(qphi[y,x],q=q,a=azi,vmin=vmin,
                  bkgd=bkgd,mask=mask,
                  vmax=vmax,log=log,cmap='jet')
        if isinstance(amin,type(None)): amin = np.min(azi)
        if isinstance(amax,type(None)): amax = np.max(azi)
        plt.vlines(qmin,amin,amax,linewidth=1.5)
        plt.vlines(qmax,amin,amax,linewidth=1.5)
    if line_plot:
        Ia = sum_roi(qphi[y,x],a=azi,q=q,qmin=qmin,qmax=qmax,mask=mask,vs_axis='a')
        plt.subplots()
        plt.plot(azi,Ia)
        plt.title(r'$reflection\,\,(0\,\,1\,\,3)')     

def asymetry_compare(x,y,qmin=1.83,qmax=1.89):
    # compare the intensity difference of maximum at (-180,0) and (0,180)
    Ia1 = sum_roi(qphi[y,x],azi,q,qmin=1.8,qmax=1.9,mask=mask,
                 amin=-180,amax=0,vs_axis='a')
    i1 = np.nanmax(Ia1)
    Ia2 = sum_roi(qphi[y,x],azi,q,qmin=1.8,qmax=1.9,mask=mask,
                 amin=0,amax=180,vs_axis='a')
    i2 = np.nanmax(Ia2)
    if i1 >= i2:
        sign = 1.
    else:
        sign = -1.
    
    return sign*np.max([i1,i2])/np.min([i1,i2])

def asym_mat(qphi,qmin,qmax):
    mat = np.zeros(qphi.shape[:2])
    for y in range(qphi.shape[0]):
        for x in range(qphi.shape[1]):
            mat[y,x] = asymetry_compare(x,y,qmin,qmax)
    return mat         


def save_as_edf(fn,data,force_type='float',fit2dMode=True):
    edf_obj = fabio.edfimage.EdfImage(data=data)
    edf_obj.write(fn,force_type=force_type,fit2dMode=fit2dMode)   

from PIL import Image
def save_as_tif(fn,data):
    img = Image.fromarray(data)
    img.save(fn,format='TIFF')
