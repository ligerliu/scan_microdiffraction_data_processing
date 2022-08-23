#get_ipython().run_line_magic('run', '-i ../../../../xs_proc/proc_data_ana.py')

fn = '../hdf_file/hum_sample_hum_sample_fiber1_c_p01_1_proc.h5'
azi,ct,qphi,q,path_list,path_idx,pttn_idx = load_proc(fn)
qphi_norm = ct_normalization(qphi,ct)
qphi_norm_new = np.copy(qphi_norm)
for _ in range(qphi_norm.shape[0]):
    for __ in range(qphi_norm.shape[1]):
        qphi_norm_new[_,__] -= (qphi_norm[_,6]+qphi_norm[_,22])/2
    
qphi_norm_cellulose = np.copy(qphi_norm_new)
del qphi_norm_new

fn = '../hdf_file/hum_sample_hum_sample_fiber2_l5_p01_1_proc.h5'
azi,ct,qphi,q,path_list,path_idx,pttn_idx = load_proc(fn)
qphi_norm = ct_normalization(qphi,ct)
qphi_norm_new = np.copy(qphi_norm)
for _ in range(qphi_norm.shape[0]):
    for __ in range(qphi_norm.shape[1]):
        qphi_norm_new[_,__] -= (qphi_norm[_,6]+qphi_norm[_,22])/2
    
qphi_norm_lignin50 = np.copy(qphi_norm_new)
del qphi_norm_new

fn = '../hdf_file/hum_sample_hum_sample_fiber3_l3_p01_1_proc.h5'
azi,ct,qphi,q,path_list,path_idx,pttn_idx = load_proc(fn)
qphi_norm = ct_normalization(qphi,ct)
qphi_norm_new = np.copy(qphi_norm)
for _ in range(qphi_norm.shape[0]):
    for __ in range(qphi_norm.shape[1]):
        qphi_norm_new[_,__] -= (qphi_norm[_,6]+qphi_norm[_,22])/2
    
qphi_norm_lignin30 = np.copy(qphi_norm_new)
del qphi_norm_new