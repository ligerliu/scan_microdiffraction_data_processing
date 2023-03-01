import sys
#sys.path.append('../../xs_proc')
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
    try:
        if np.sign(obj.ndetx) <= 0:
            det_n = 'm{}'.format(int(np.abs(obj.ndetx)))
        else:
            det_n = 'p{}'.format(int(np.abs(obj.ndetx)))
    except:
        pass

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

#def save_proc_h5(num,h5_path_list,name,proc_h5_folder,
#                idx_list,res):
#       proc_h5_name = "{}_{:05d}_proc.h5".format(name,num)
#       proc_h5_name = os.path.join(proc_h5_folder,proc_h5_name)
#       with h5py.File(proc_h5_name,'a') as k:
#           if "integrate2d" in list(k):
#               del k["integrate2d"]
#           k.create_group("integrate2d")
#           qphi_data = []
#           for __ in range(idx_list[num][0],idx_list[num][1]):
#               qphi_data.append(res[__][0])
#           k["integrate2d"].create_dataset("map_qphi",
#                                   data = np.array(qphi_data),
#                                   compression="gzip",
#                                   compression_opts=9)
#       return proc_h5_name

def save_qphi_as_h5(obj,save_path,name,
                    q,
                    azi,
                    res,
                    path_idx,
                    pttn_idx,
                    h5_path_list,
                    total_pttn_num,
                    single_h5_pttn_num):
    hdf_folder = os.path.join(save_path,"hdf_file")
    if not os.path.exists(hdf_folder):
        os.mkdir(hdf_folder)
    #print("\nhdf file writinge start time:\n%5.2f sec" %(time.time()-t))
    h5_name = os.path.join(hdf_folder,"{}_proc.h5".format(name))
    with h5py.File(h5_name,'a') as f:
        if "integrate2d" in list(f):
            del f["integrate2d"]
        f.create_group("integrate2d")
        fc = f['integrate2d']
        fc.create_dataset("beam_intensity",data=obj.ct34)
        fc.create_dataset("q",data=q)
        fc.create_dataset("angle",data=azi)
        fc.attrs["origin_h5_path"]=h5_path_list
        fc.create_dataset("path_idx",data=path_idx)
        fc.create_dataset("pttn_idx",data=pttn_idx)
        #fc.create_dataset("detector_distance",data=obj.ndetx)
        proc_h5_folder = os.path.join(hdf_folder,name)
        if not os.path.exists(proc_h5_folder):
            os.mkdir(proc_h5_folder)
        idx_list = chunk_idx(total_pttn_num,single_h5_pttn_num)
        proc_h5_name_list = []
        # here it works, but should be parallelized to accelerate
        #name_res = parallel_func(save_proc_h5,12,
        #            np.arange(len(h5_path_list)),
        #            h5_path_list=h5_path_list,name=name,
        #            proc_h5_folder=proc_h5_folder,
        #            idx_list=idx_list,res=res)
        #for _ in name_res:
        #    proc_h5_name_list.append(_)
        for _ in range(len(h5_path_list)):
            proc_h5_name = "{}_{:05d}_proc.h5".format(name,_)
            proc_h5_name = os.path.join(proc_h5_folder,proc_h5_name)
            proc_h5_name_list.append(proc_h5_name)
            with h5py.File(proc_h5_name,'a') as k:
                if "integrate2d" in list(k):
                    del k["integrate2d"]
                k.create_group("integrate2d")
                qphi_data = []
                for __ in range(idx_list[_][0],idx_list[_][1]):
                    qphi_data.append(res[__][0])
                k["integrate2d"].create_dataset("map_qphi",
                                        data = np.array(qphi_data),
                                        compression="gzip",
                                        compression_opts=9)
        fc.attrs["proc_h5_list"] = proc_h5_name_list
    print('file save finished')

def auto_proc_qphi(obj,
              samples,
              save_path,
              num_core= 8,
              mask    = None,
              poni    = None,
              q_npts  = 200,
              a_npts  = 120,
              save    = True,
              data_path = 'entry_0000/measurement/data',
              radial_range = None,
              **kwargs
              ):
    failed_pattern = []
    for ii in samples.keys():
        try:
            obj.keyword_search(*samples[ii])
            poni_file,mask_file = jl_pyFAI_setup(obj,poni=poni,
                                  mask=mask,save_path=save_path)
            total_pttns,scan_shape,idx_list = scan_info(obj)
            h5_list,path_idx,pttn_idx = scan_h5_data_info(obj,scan_shape,idx_list)
            ai = pyFAI.load(poni_file)    
            tm = time.time()
            
            if not isinstance(obj.ct34, type(None)):
                ct34 = np.copy(obj.ct34)
                cur_pos = np.argwhere(ct34<(np.nanmean(ct34)*0.6)) 
                ct34[cur_pos] = ct34[cur_pos-1]
                if np.nanmean(ct34) <= 0:
                    ct34 = np.ones((len(pttn_idx.flatten()),))
                else:
                    ct34 /= np.nanmean(ct34)
            else:
                ct34 = np.ones((len(pttn_idx.flatten()),))
            
            res = parallel_func(scan_calculate_Iqphi,
                           num_core,
                           np.arange(len(path_idx.flatten())),
                           h5_list   = h5_list,
                           path_idx  = path_idx.flatten(),
                           pttn_idx  = pttn_idx.flatten(),
                           data_path = data_path,
                           pyfai_obj = ai,
                           mask      = mask_file,
                           q_npts    = q_npts,
                           a_npts    = a_npts,
                           ct        = ct34,
                           radial_range = radial_range,
                           **kwargs
                           )
            print(time.time()-tm)
            q    = res[0][1]
            azi  = res[0][2]
            name = obj._data_name[0].split('.')[0]
            if save:
                save_qphi_as_h5(obj,save_path,name,
                    q    = q,
                    azi  = azi,
                    res = res,
                    path_idx = path_idx,
                    pttn_idx = pttn_idx,
                    h5_path_list = h5_list,
                    total_pttn_num = total_pttns,
                    single_h5_pttn_num = obj.single_h5_shape[0])
            print(time.time()-tm)
            #break
        except Exception as e:
            print(e)
            failed_pattern.append(obj._data_name)
            pass
    print(failed_pattern) 
