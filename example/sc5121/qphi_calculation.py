import sys
sys.path.append('../../xs_proc/')
from visual_func import *
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
from multi_scan_Iq_proc import *

proposal_path = "/data/visitor/sc5121/id13/"
proposal_num  = "sc5121"
data_path     = 'entry_0000/measurement/data'

save_path     = os.getcwd()

t1 = id13_h5_search(proposal_path,proposal_num)

sample_keywords = [
                   #['cell1_exsitu_in55_26.1'],
                   #['cell1_test1_3.1'],
                   #['cell1_exsitu_in55_10.1'],
                   ['cell1','exsitu','horgeo','wet','4.1']
                   ]

# get hpyer link for all the sample keywords, more keywords could be included for locating specific pattern
samples = samples_for_process(t1,sample_keywords)
tm = time.time()
#print(samples.keys())
#sys.exit()
auto_proc_qphi(t1,samples,save_path,num_core=16,q_npts=1200,a_npts=360)
