import sys
sys.path.append('../../xs_proc')
from xs_data_proc import *
from proc_data_ana import *
import matplotlib.pyplot as plt
from h5_data_search import id13_h5_search
import pyFAI
import numpy as np
import h5py
import os
import fabio
import time
import re
import warnings
warnings.filterwarnings('ignore')

def jl_pyFAI_setup(obj,
                   poni      = None,
                   mask      = None,
                   save_path = None,
                   ):
    # obj is the "id13_h5_search object includes id13 h5 data info
    if isinstance(save_path,type(None)):
        save_path = os.getcwd()
    if np.sign(obj.ndetx) <= 0:
        det_n = 'm{}'.format(int(np.abs(obj.ndetx)))
    else:
        det_n = 'p{}'.format(int(np.abs(obj.ndetx)))

    if isinstance(mask,type(None)):
        if os.path.isfile('{}/msk_file/ndet_{}.npz'.format(save_path,det_n)):
            mask = np.load('{}/msk_file/ndet_{}.npz'.format(save_path,det_n))['mask']
    else:
        mask = mask
    if isinstance(poni,type(None)):
        if os.path.isfile('{}/poni_file/ndet_{}.poni'.format(save_path,det_n)):
            poni = "{}/poni_file/ndet_{}.poni".format(save_path,det_n)
        else:
            print('poni file is needed for pyFAI processing')
            sys.exit()
    else:
        poni = poni
    return poni,mask

def save_qphi_as_h5(obj,save_path,name,
                    q,
                    azi,
                    qphi,
                    path_idx,
                    pttn_idx,
                    h5_path_list):
    os.chdir(save_path)
    if not os.path.exists("hdf_file"):
        os.mkdir("hdf_file")
    os.chdir("hdf_file")
    #print("\nhdf file writinge start time:\n%5.2f sec" %(time.time()-t))
    h5_name = "{}_proc.h5".format(name)
    fc = h5py.File(h5_name,"w")
    fc.create_dataset("beam_intensity",data=obj.ct34)
    fc.create_dataset("q",data=q)
    fc.create_dataset("angle",data=azi)
    fc.create_dataset("map_qphi",
                      data=qphi,compression="gzip",
                      #compression_opts=9
                     )
    fc.attrs["origin_h5_path"]=h5_path_list
    fc.create_dataset("path_idx",data=path_idx)
    fc.create_dataset("pttn_idx",data=pttn_idx)
    fc.create_dataset("detector_distance",data=obj.ndetx)
    fc.close()

def auto_proc(obj,
              samples,
              save_path,
              num_core= 8,
              mask    = None,
              poni    = None,
              q_npts  = 200,
              a_npts  = 120,
              save    = True,
              data_path = 'entry_0000/measurement/data',
              **kwargs
              ):
    failed_pattern = []
    for ii in samples.keys():
        obj.keyword_search(*samples[ii])
        poni,mask = jl_pyFAI_setup(obj,poni=poni,mask=mask,save_path=save_path)
        total_pttns,scan_shape,idx_list = scan_info(obj)
        h5_list,path_idx,pttn_idx = scan_h5_data_info(obj,scan_shape,idx_list)
        ai = pyFAI.load(poni)    
        #tm = time.time()
        try:
            res = parallel_func(scan_calculate_Iqphi,
                           num_core,
                           np.arange(len(path_idx.flatten())),
                           h5_list   = h5_list,
                           path_idx  = path_idx.flatten(),
                           pttn_idx  = pttn_idx.flatten(),
                           data_path = data_path,
                           pyfai_obj = ai,
                           mask      = mask,
                           q_npts    = q_npts,
                           a_npts    = a_npts,
                           **kwargs
                           )
            #print(time.time()-tm)
            q    = res[0][1]
            azi  = res[0][2]
            qphi = np.zeros((scan_shape[0],
                             scan_shape[1],
                             res[0][0].shape[0],
                             res[0][0].shape[1])) 
            for _ in range(len(res)):
                if len(idx_list) > 1:
                    i1 = int((_+i*t1.single_h5_shape[0])/scan_shape[1])
                    i2 = int((_+i*t1.single_h5_shape[0])%scan_shape[1])
                else:
                    i1 = int(_/scan_shape[1])
                    i2 = int(_%scan_shape[1])
                qphi[i1,i2,:]   = res[_][0]
            name = obj._data_name[0].split('.')[0]
            if save:
                save_qphi_as_h5(obj,save_path,name,
                    q    = q,
                    azi  = azi,
                    qphi = qphi,
                    path_idx = path_idx,
                    pttn_idx = pttn_idx,
                    h5_path_list = h5_list)
            break
        except:
            failed_pattern.append(obj._data_name)
            pass
        print(failed_pattern) 
