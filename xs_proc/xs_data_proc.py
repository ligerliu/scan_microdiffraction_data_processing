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

############### additional code 

def samples_for_process(obj,kw_list):
    # obj is the "id13_h5_search" object include dataset structure information, defined in the "h5_data_search"
    # kw_list is key word list, should have form like: [kw] for one sample, 
    # or [[kw1],[kw2]...] for multiple samples
    sample_num = 1
    samples    = {}
    for kw in kw_list:
        if isinstance(kw,list):
            obj.keyword_search(*kw)
        elif isinstance(kw,str):
            obj.keyword_search(kw)
        
        for _ in range(len(obj._data_name)):
            samples['sample{}'.format(sample_num)] = [obj._data_name[_]]
        
        sample_num += 1
        
    return samples

def scan_info(obj):
    # obj is the "id13_h5_search object includes id13 h5 data info
    h5_path_list = []
    h5_shape     = []
    h5_path_list = obj.data_h5
    h5_shape.append(obj.shape)
    
    if np.size(h5_shape) == 1:
        if np.size(obj.ct34) == h5_shape[0][0]:
            total_pttns = h5_shape[0][0]
            scan_shape  = [total_pttns,1]
        elif np.size(obj.ct34) == (h5_shape[0][0]+1):
            total_pttns = h5_shape[0][0]+1
            scan_shape  = [total_pttns,1]
        else:
            print("scan shape is inconsisted with planned")
            #sys.exit()
            total_pttns = np.size(obj.ct34)
            if np.size(obj.ct34) == h5_shape[0][0]:
                scan_shape = [h5_shape[0][0],1]
            elif np.size(obj.ct34) == (h5_shape[0][0]+1):
                scan_shape = [h5_shape[0][0]+1,1]
    elif np.size(h5_shape) == 2:
        if np.size(obj.ct34) == h5_shape[0][0]*h5_shape[0][1]:
            total_pttns = h5_shape[0][0]*h5_shape[0][1]
            scan_shape  = h5_shape[0]
        elif  np.size(obj.ct34) == (h5_shape[0][0]+1)*(h5_shape[0][1]+1):
            total_pttns = (h5_shape[0][0]+1)*(h5_shape[0][1]+1)
            scan_shape  = (h5_shape[0]+1)
        else:
            print("scan shape is inconsisted with planned")
            #sys.exit()
            total_pttns = np.size(obj.ct34)
            if np.size(obj.ct34) == h5_shape[0][0]*h5_shape[0][1]:
                scan_shape = [h5_shape[0][0],1]
            elif np.size(obj.ct34) == (h5_shape[0][0]+1)*(h5_shape[0][1]+1):
                scan_shape = [h5_shape[0][0]+1,1]
    
    # a soft problem here, we assum every h5 file is not including too much images. 
    # If the h5 is too large to load at once, this chunk size will cause memory problem.
    # Here JL assume the frame per file will be less than the default for h5 compression, 2000 per file.
    # a soft solution is try to further chunk the data for each h5 file data
    if np.size(h5_shape) == 1:
        if obj.single_h5_shape[0] < (total_pttns):
            idx_list = chunk_idx((total_pttns),
                                  obj.single_h5_shape[0])
        else:
            idx_list = [[0,(total_pttns)]]
    if np.size(h5_shape) == 2:
        if obj.single_h5_shape[0] < (total_pttns):
            idx_list    = chunk_idx((total_pttns),
                                    obj.single_h5_shape[0])
        else:
            idx_list    = [[0,(total_pttns)]]

    if len(idx_list) > len(h5_path_list):
        # since the total pttn had been determined by the ion chamber counting, this should never happen
        # these give a 
        print('\n',"the scan plan to collect {} patterns".format(total_pttns))
        print('\n',"but only {} h5 files, each contains {}".format(len(h5_path_list),
                                                                  obj.single_h5_shape[0]))
        print('\n','collected patterns less than required patterns, the scan is incomplete')
        sys.exit()

    return total_pttns,scan_shape,idx_list

def scan_h5_data_info(obj,scan_shape,idx_list):
    path_idx = np.zeros((scan_shape[0]*scan_shape[1]),dtype=np.int16)
    pttn_idx = np.zeros((scan_shape[0]*scan_shape[1]),dtype=np.int16)
    for i in range(len(idx_list)):
        if len(idx_list) > 1:
            idx1 = idx_list[i][0] - i*obj.single_h5_shape[0]
            idx2 = idx_list[i][1] - i*obj.single_h5_shape[0]
        else:
            idx1 = idx_list[i][0]
            idx2 = idx_list[i][1]
        path_idx[idx_list[i][0]:idx_list[i][1]] = i
        pttn_idx[idx_list[i][0]:idx_list[i][1]] = np.arange(idx1,idx2).astype(np.int16)
    path_idx = path_idx.reshape((scan_shape[0],scan_shape[1]))
    pttn_idx = pttn_idx.reshape((scan_shape[0],scan_shape[1]))
    return  obj.data_h5,path_idx,pttn_idx

def scan_calculate_Iqphi(num, 
                         h5_list  = None,
                         path_idx = None,
                         pttn_idx = None,
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
    h5_path  = h5_list[path_idx[num]]
    pttn_num = pttn_idx[num]
    
    with h5py.File(h5_path,"r") as f:
        data = f[data_path][pttn_num]
        qphi,q,a = pyfai_obj.integrate2d(
                                    data,
                                    npt_rad  = q_npts,
                                    npt_azim = a_npts,
                                    unit     = q_unit,
                                    method   = method,
                                    **kwargs)
    return qphi,q,a

def scan_calculate_Iq(num,
                 h5_list   = None,
                 path_idx  = None,
                 pttn_idx  = None,
                 data_path = None,
                 pyfai_obj = None,
                 q_npts    = 300 ,
                 **kwargs):
    h5_path  = h5_list[path_idx[num]]
    pttn_num = pttn_idx[num]
    with h5py.File(h5_path,"r") as f:
        data = f[data_path][pttn_num]
        q,I  = pyfai_obj.integrate1d(data,
                                    q_npts,
                                    **kwargs)
    return q,I

def scan_pttn_roi_sum(num,
                 h5_list   = None,
                 path_idx  = None,
                 pttn_idx  = None,
                 data_path = None,
                 left_top     = [0,0],
                 right_bottom = [-1,-1] ,
                 **kwargs):
    h5_path  = h5_list[path_idx[num]]
    pttn_num = pttn_idx[num]
    with h5py.File(h5_path,"r") as f:
        data = f[data_path][pttn_num][left_top[1]:right_bottom[1],
                                 left_top[0]:right_bottom[0]].astype(float)
        data[data>1e6] = np.nan
        I = np.nanmean(data)
    return I
