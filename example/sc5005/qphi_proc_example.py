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
from multi_scan_qphi_proc import *

if __name__ == '__main__':
    proposal_path = "/data/visitor/sc5005/id13"
    proposal_num  = "sc5005"
    data_path     = 'entry_0000/measurement/data'
    # the directory could be replaced by user with specific path
    save_path     = os.getcwd()
    
    t1 = id13_h5_search(proposal_path,proposal_num)
    
    sample_keywords = [
                  ['BM1N'],
                  ['BM1D'],
                  #['mount01_DS_a_mount01_x50_M1_BM1N_p01_0013_2.1'],
                  ]
    samples = samples_for_process(t1,sample_keywords)
    tm = time.time()
    auto_proc(t1,samples,save_path,num_core=16)
    print(time.time()-tm)