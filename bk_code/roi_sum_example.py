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

proposal_path = "/data/visitor/ls2995/id13"
proposal_num  = "ls2995"
data_path     = 'entry_0000/measurement/data'
# the directory could be replaced by user with specific path
save_path     = os.getcwd()


t1 = id13_h5_search(proposal_path,proposal_num)


sample_keywords = [
                  #['MoBS_A_C1_9.1'],
                  #['WTCS1','8'],
                  #['WTCS1_b_mac1_strip1lr1_1.1'],
                  #['WTCS1_b_mac1_strip1lr1_2.1'],
                  ['WTCS1_c_mac1_strip1lr1_1.1'],
                  #['mac1','sup']
                  ]
samples = samples_for_process(t1,sample_keywords)

total_pttns,scan_shape,idx_list = scan_info(t1)
h5_list,path_idx,pttn_idx = scan_h5_data_info(t1,scan_shape,idx_list)

#import time
#tm = time.time()
#res = parallel_func(scan_pttn_roi_sum,12,np.arange(len(path_idx.flatten())),
#                   h5_list=h5_list,
#                   path_idx=path_idx.flatten(),
#                   pttn_idx=pttn_idx.flatten(),
#                   data_path=data_path,
#                   left_top=(1068,1140),
#                   right_bottom=(1222,1322),)
#print(time.time()-tm)
#roi = np.array(res).reshape(scan_shape)
