import sys
sys.path.append('../xs_proc/')
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

proposal_path = "/gpfs/jazzy/data/id13/inhouse3/DATAPOLICY_I3_1/eh3/inhouse/ihsc1698/id13"
proposal_num  = "ihsc1698"
data_path     = 'entry_0000/measurement/data'

save_path = '/data/id13/inhouse12/jiliang/code_v6/example/ihsc1698'


poni = '/data/id13/inhouse12/jiliang/code_v6/example/ihsc1698/poni_file/ndetx_m300.poni'
mask = np.load('/data/id13/inhouse12/jiliang/code_v6/example/ihsc1698/msk_file/ndet_m300_m1.npz')['mask']

t1 = id13_h5_search(proposal_path,proposal_num)
tim = time.time()
sample_keywords = [
['calibrant_Al2O3_ndet_m300_2_1'],    
]
samples = samples_for_process(t1,sample_keywords)
auto_proc_qphi(t1,samples,save_path,num_core=25,q_npts=720,poni=poni,mask=mask)

print(f'\n\n{time.time()-tim}')
