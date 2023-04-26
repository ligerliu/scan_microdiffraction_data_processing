#import sys
#sys.path.append('./')
from xs_data_proc import calculate_Iqphi,chunk_idx,parallel_func,pttn_roi_sum
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

def auto_load_xrf(obj,
              #samples,
              #save_path,
              #save    = False,
              output  = True,
              element    = None,
              energy_roi = (0,40),
              energy_calibration = 5e-3,
              **kwargs
             ):
    '''
    func is process funciton, for 2D integral, we use "calculate_Iqphi" from "xs_scan_proc.py"
    obj is the "id13_h5_search" object
    '''
    t = time.time()
    nn = 0
    h5_path_list = []
    h5_shape     = []
    # here my personal habit for saving the poni and mask file
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
            print("scan shape is inconsisted with planned, scan may not be complete.")
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
            print("scan shape is inconsisted with planned, scan may not be complete.")
            #sys.exit()
            total_pttns = np.size(obj.ct34)
            if np.size(obj.ct34) == h5_shape[0][0]*h5_shape[0][1]: 
                scan_shape = [h5_shape[0][0],1]
            elif np.size(obj.ct34) == (h5_shape[0][0]+1)*(h5_shape[0][1]+1):
                scan_shape = [h5_shape[0][0]+1,1]
    
    
    with h5py.File(obj.fn,'r') as f:
        if ('xmap3_det0' in f[obj._data_name[nn]]['measurement']):
            XRF_m = np.array(f[obj._data_name[nn]]['measurement']['xmap3_det0'])
            energy = np.arange(XRF_m.shape[1])*energy_calibration
            XRF_m = np.nanmean(XRF_m[:,((energy>=energy_roi[0])&(energy<=energy_roi[1]))],axis=-1)
            if len(XRF_m) < scan_shape[0]*scan_shape[1]:
                XRF_roi = np.zeros((scan_shape[0]*scan_shape[1],))
                XRF_roi[:len(XRF_m)] = XRF_m
            else:
                XRF_roi = np.copy(XRF_m) 
            XRF_roi = XRF_roi.reshape(scan_shape[0],scan_shape[1])
            #XRF_roi = np.nanmean(XRF_m[:,:,((energy>=energy_roi[0])&(energy<=energy_roi[1]))],axis=-1)
        elif ('xmap2_det0' in f[obj._data_name[nn]]['measurement']):
            XRF_m =  np.array(f[obj._data_name[nn]]['measurement']['xmap2_det0'])
            energy = np.arange(XRF_m.shape[1])*energy_calibration
            XRF_m = np.nanmean(XRF_m[:,((energy>=energy_roi[0])&(energy<=energy_roi[1]))],axis=-1)
            if len(XRF_m) < scan_shape[0]*scan_shape[1]:
                XRF_roi = np.zeros((scan_shape[0]*scan_shape[1],))
                XRF_roi[:len(XRF_m)] = XRF_m 
            else:
                XRF_roi = np.copy(XRF_m)
            XRF_roi = XRF_roi.reshape(scan_shape[0],scan_shape[1])
            #XRF_m = XRF_m.reshape(scan_shape[0],scan_shape[1],XRF_m.shape[1])
            #XRF_roi = np.nanmean(XRF_m[:,:,((energy>=energy_roi[0])&(energy<=energy_roi[1]))],axis=-1)
        else:
            print('xmap3 and xmap2 are not define in the measurement,/n please ensure xmap is included')
            return 
        #if not isinstance(element,type(None)):
        #    if (element in f[obj._data_name[nn]]['measurement']):
        #        XRF_e = np.array(f[obj._data_name[nn]]['measurement'][element])
        #        XRF_e = XRF_e.reshape(scan_shape[0],scan_shape[1])
        #    else:
        #        print('input element name is not correct, please find correct name below')
        #        print('/n',np.array(f[obj._data_name[nn]]['measurement']))
        #        return 
        #else:
        #    XRF_e = np.zeros(scan_shape)
    #nn += 1
 
    #if save:
    #    os.chdir(save_path)
    #    if not os.path.exists("xrf_file"):
    #        os.mkdir("xrf_file")
    #    os.chdir("xf_file")
    #    print("\nhdf file writinge start time:\n%5.2f sec" %(time.time()-t))
        
    if output == True:
        #print(XRF_m.shape,XRF_m[:,:,1305],XRF_m[:,:,100],energy_roi)
        return XRF_roi,energy#XRF_m,XRF_roi,XRF_e,h5_path_list,path_idx,pttn_idx
    
            
            
