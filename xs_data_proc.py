import numpy as np
import h5py
import pyFAI
from multiprocessing import Pool
from functools import partial
import warnings
warnings.filterwarnings("ignore")#ignore warnings

#-------------------------------------------------------------------
# try to calcualte the 2D azimuthal integral for further processing, 
# azimuthal integral was conducted by the pyFAI function
# calculation is cpu paralllelized by multiprocessing.Pool
#-------------------------------------------------------------------


def chunk_idx(total_num,slice_num):
    total_num = int(total_num)
    slice_num = int(slice_num)
    idx  =  np.arange(0,total_num)[slice(0,total_num,slice_num)]
    idx1 =  idx
    idx2 =  np.append(idx[1:],total_num)
    #idx2.append(total_nem)
    return list(zip(idx1,idx2))

def partial_sum(chunk,
                h5_path   = None,
                data_path = None):
    start,end = chunk
    with h5py.File(h5_path,"r") as f:
        data = f[data_path][start:end]
        data = np.sum(data,axis=0)
    return data

def partial_mean(chunk,
                h5_path   = None,
                data_path = None):
    start,end = chunk
    with h5py.File(h5_path,"r") as f:
        data = f[data_path][start:end]
        data = np.mean(data,axis=0)
    return data

def calculate_Iqphi(num, 
                   h5_path = None,
                   data_path = None,
                   pyfai_obj = None,
                   method  = None,
                   q_npts=300,
                   a_npts=180,
                   q_unit="q_A^-1",
                   **kwargs):
    if isinstance(method,type(None)):
        method = 'csr'
    else:
        method = method
    with h5py.File(h5_path,"r") as f:
        data = f[data_path][num]
        qphi,q,a = pyfai_obj.integrate2d(
                                    data,
                                    npt_rad  = q_npts,
                                    npt_azim = a_npts,
                                    unit     = q_unit,
                                    method   = method,
                                    **kwargs)
    return qphi,q,a

def calculate_Iq(num,
                 h5_path   = None,
                 data_path = None,
                 pyfai_obj = None,
                 q_npts    = 300 ,
                 **kwargs):
    with h5py.File(h5_path,"r") as f:
        data = f[data_path][num]
        q,I  = pyfai_obj.integrate1d(data,
                                    q_npts,
                                    **kwargs)
    return q,I

def pttn_roi_sum(num,
                 h5_path   = None,
                 data_path = None,
                 left_top     = [0,0],
                 right_bottom = [-1,-1] ,
                 **kwargs):
    with h5py.File(h5_path,"r") as f:
        data = f[data_path][num][left_top[1]:right_bottom[1],
                                 left_top[0]:right_bottom[0]].astype(float)
        data[data>1e6] = np.nan
        I = np.nanmean(data)
    return I
   
def parallel_func(func,num_cores,args,**kwargs):
    func = partial(func,**kwargs)
    with Pool(num_cores) as pool:
        res = pool.map(func,args) 
        pool.close()
    return res
