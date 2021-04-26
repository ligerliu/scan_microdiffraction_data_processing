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

proposal_path = '/data/visitor/sc5121/id13/'
proposal_num  = 'sc5121'
data_path     = 'entry_0000/measurement/data'
# the directory could be replaced by user with specific path
save_path     = os.getcwd()


t1 = id13_h5_search(proposal_path,proposal_num)


sample_keywords = [
                  #['cell1_exsitu_in55_26.1'],
                  #['cell1_test1_3.1'],
                  #['cell1_exsitu_in55_10.1'],
                  ['cell1_exsitu_horgeo_in55_wet_4.1'],
                  ]
samples = samples_for_process(t1,sample_keywords)

total_pttns,scan_shape,idx_list = scan_info(t1)
h5_list,path_idx,pttn_idx = scan_h5_data_info(t1,scan_shape,idx_list)

def bkgd_cal(left_top,bottom_right,h5_list,path_idx,pttn_idx,
             data_path= 'entry_0000/measurement/data'):
    #left_top and bottom_right are the coordinate with (x,y) or (col,row)
    data = []
    for i in range(left_top[1],(bottom_right[1]+1)):
        for j in range(left_top[0],(bottom_right[0]+1)):
            with h5py.File(h5_list[path_idx[i,j]],'r') as f:
                data.append(f[data_path][pttn_idx[i,j]])
    data = np.array(data)
    data[data>1e6] = 0
    return np.nanmean(data,axis=0)

bkgd = bkgd_cal((0,250),(0,300),h5_list,path_idx,pttn_idx)

import time
tm = time.time()
res = parallel_func(scan_pttn_roi_sum,12,np.arange(len(path_idx.flatten())),
                   h5_list=h5_list,
                   path_idx=path_idx.flatten(),
                   pttn_idx=pttn_idx.flatten(),
                   data_path=data_path,
                   left_top=(1150,1150),
                   right_bottom=(1450,1450),
                   bkgd=bkgd,
                   #mask = mask[1150:1450,1150:1450],
                   )
print(time.time()-tm)
roi = np.array(res).reshape(scan_shape)
#pttn_of_int_map(roi,h5_list,path_idx,pttn_idx,log=False)
