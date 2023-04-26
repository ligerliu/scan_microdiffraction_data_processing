import sys
sys.path.append('/data/id13/inhouse12/jiliang/code_v7/xs_proc/')
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

proposal_path = "/data/visitor/ihls3485/id13/20230307/"
proposal_num  = "ihls3485"
data_path     = 'entry_0000/measurement/data'



# for saving integration results, save_path is needed. save_path is defined by user

save_path = '/data/visitor/ihls3485/id13/XXPROCESS/'

poni = '/data/visitor/ihls3485/id13/XXPROCESS/poni_file/udetx_m820.poni'
msk_file = '/data/visitor/ihls3485/id13/XXPROCESS/edf_sum_file/udetx_m820.edf'
mask    = fabio.open(msk_file).data
#mask = np.savez(msk_file)

#mask    = np.load(msk_file)['mask']
processed_list = os.listdir(os.path.join(save_path,'hdf_file'))
t1 = id13_h5_search(proposal_path,proposal_num)

tim = time.time()

sample_keywords = [
    ['DG'],
    ['EC'],
]

samples = samples_for_process(t1,sample_keywords)

samples_key_list = list(samples.keys())

for _ in range(len(samples)):

    if samples[samples_key_list[_]][0][:-2] in processed_list:

        del samples[samples_key_list[_]]


tim = time.time()

auto_proc_qphi(t1,samples,save_path,num_core=24,q_npts=720, a_npts=180,poni=poni,mask=mask,)#radial_range=(0,0.2))
print('qphi calculation finished')
auto_proc_Iq(t1,samples,save_path,num_core=24,q_npts=720,poni=poni,mask=mask,)#radial_range=(0,2))
print('Iq calculation finished')
print(f'\n\n{time.time()-tim}')
