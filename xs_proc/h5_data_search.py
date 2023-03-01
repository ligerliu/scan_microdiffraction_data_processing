import numpy as np
import h5py
import os
import sys
import re
import glob
import warnings
warnings.filterwarnings("ignore")#ignore warnings

# ------------------------------------------------------------------
# try to locate the 
# azimuthal integral was conducted by the pyFAI function
# ------------------------------------------------------------------


class id13_h5_search:
    def __init__(self,
                proposal_path = None,
                proposal      = None,
                data_h5path   = "entry_0000/measurement/data",
                file_sys      = "linux",
                detector      ='eiger',
                ):
        if proposal_path == None:
            self.proposal_path = os.getcwd()
        else:
            if os.path.exists(proposal_path):
                self.proposal_path = proposal_path
            else:
                raise ValueError("input proposal path is not correct, pls check")
            
        if proposal == None:
            raise ValueError("Input proposal is not correct")
        else:
            self.fn = "{}/{}_id13.h5".format(self.proposal_path,proposal)
            self.check_path(self.fn,'proposal {}'.format(proposal))
        self.file_sys = file_sys
        self.data_h5path = data_h5path
        self.detector = detector
    
    def check_path(self,path,name):
        if not os.path.exists(path):
            raise ValueError("input path :{}, is not correct path for {}".format(path,name))
    
    def keyword_search(self,*args,output_whole_list=False):
        keywords = []
        for _ in args:
            keywords.append(str(_))
        num = 0
        self._data_name = []
        with h5py.File(self.fn,"r") as f:
            print('found dataset:')
            for _ in list(f):
                if all(str(x) in _ for x in keywords):
                    print('{}'.format(_))
                    self._data_name.append(_)
                    num += 1
            if num == 0:
                print("No data found, One or more keywods are incorrect, please check input of proposal, sample, dataset, scan")
                return
            elif num == 1:
                p  = re.split('_',self._data_name[0])
                for path in os.listdir(self.proposal_path):
                    if any(x in path for x in p):
                        subdir = os.path.join(self.proposal_path,path)
                        if os.path.isdir(subdir):
                            for subpath in os.listdir(subdir):
                                if str.join("_",p[:-1]) == subpath:
                                    self.sample  = path
                                    self.dataset = subpath
                                    self.scannum = "scan%05d" %int(
                                                  p[-1].split('.')[0])
                                                  
                                    self.data_path = os.path.join(
                                                self.proposal_path,
                                                self.sample,
                                                self.dataset,
                                                self.scannum,
                                                self.detector)
                                    if os.path.isdir(self.data_path):
                                        pass
                                    else:
                                        self.scannum = "scan%04d" %int(
                                                      p[-1].split('.')[0])
                                        
                                        self.data_path = os.path.join(
                                                    self.proposal_path,
                                                    self.sample,
                                                    self.dataset,
                                                    self.scannum,
                                                    self.detector)
                                    print('\nsample:\n{}'.format(path))
                                    print('\ndataset:\n{}'.format(subpath))
                                    print('\nscan:\n{}'.format(self.scannum))
                self.data_info = {}
                a  = np.array(f["{}/title".format(
                                    self._data_name[0])])
                if type(a.item()) == bytes:
                    a = a.item().decode("utf-8")
                elif isinstance(a.item(),u''.__class__):
                    a = a.item()
                elif type(a.item()) == str:
                    a = a.item()
                u  = re.split(" ",a)
                if u[0][-1] != '(':
                    self.data_info = {"scan_type":u[0]}
                elif u[0][-1] == '(':
                    self.data_info = {"scan_type":u[0][:-1]}
                    
                #should revise ct34 to beam intensity and ndetx to detx
                try:
                    try:
                        self.ct34 = np.array(
                            f["{}/measurement/ct34".format(self._data_name[0])])
                    except:
                        self.ct34 = np.array(
                           f["{}/measurement/ct24".format(self._data_name[0])])
                except:
                    self.ct34  = None
                
                try:                
                    try:
                        self.ndetx = np.array(
                            f["{}/instrument/positioners/ndetx".format(self._data_name[0])])
                    except:
                        self.ndetx = np.array(
                            f["{}/instrument/positioners/udetx".format(self._data_name[0])])
                except:
                    self.ndetx = None
                
                if self.data_info["scan_type"] == u'akmap':
                    self.data_info["fast_axis"]=u[1][:-1]
                    self.data_info["fast_axis_start_pos"]=float(u[2][:-1])
                    self.data_info["fast_axis_end_pos"]=float(u[3][:-1])
                    self.data_info["fast_axis_num_step"]=int(u[4][:-1])
                    self.data_info["slow_axis"]=u[5][:-1]
                    self.data_info["slow_axis_start_pos"]=float(u[6][:-1])
                    self.data_info["slow_axis_end_pos"]=float(u[7][:-1])
                    self.data_info["slow_axis_num_step"]=int(u[8][:-1])
                    self.data_info["exposure_time"]=float(u[9][:])
                    self.shape = np.array((int(u[8][:-1]),
                                           int(u[4][:-1])))
                elif self.data_info["scan_type"] == u'ct':
                    self.data_info["exposure_time"]=u[1]
                    self.shape = np.array((1,))
                elif self.data_info["scan_type"] == u'loopscan':
                    #self.data_info["exposure_time"]=u[1]
                    #self.shape = np.array((1,))
                    try:
                        self.data_info["exposure_time"] = f["{}/measurement/acq_time_3".format(self._data_name[0])][0]
                        self.shape = np.array((len(f["{}/measurement/acq_time_3".format(self._data_name[0])]),))
                    except:    
                        self.data_info["exposure_time"] = f["{}/measurement/acq_time_2".format(self._data_name[0])][0]
                        self.shape = np.array((len(f["{}/measurement/acq_time_2".format(self._data_name[0])]),))
                elif (self.data_info["scan_type"] == u'dscan') or \
                     (self.data_info["scan_type"] == u'ascan'):
                    self.data_info["scan_axis"]= u[1]
                    self.data_info["start_pos"]= u[2]  
                    self.data_info["end_pos"  ]= u[3]
                    self.data_info["num_step" ]= u[4]
                    self.data_info["exposure_time"]= u[5] 
                    self.shape = np.array((int(u[4]),))
                elif (self.data_info["scan_type"]  == u'dmesh') or \
                     (self.data_info["scan_type"]  == u'amesh'):
                    self.data_info["axis1"]=u[1]
                    self.data_info["axis1_start_pos"]=float(u[2])
                    self.data_info["axis1_end_pos"]=float(u[3])
                    self.data_info["axis1_num_step"]=int(u[4])
                    self.data_info["axis2"]=u[5]
                    self.data_info["axis2_start_pos"]=float(u[6])
                    self.data_info["axis2_end_pos"]=float(u[7])
                    self.data_info["axis2_num_step"]=int(u[8])
                    self.data_info["exposure_time"]=float(u[9][:])
                    self.shape = np.array((int(u[8]),
                                           int(u[4])))
                elif self.data_info["scan_type"]  == u'akmap_lut':
                    self.data_info["fast_axis"]=u[1][:-1]
                    self.data_info["fast_axis_start_pos"]=float(u[2][:-1])
                    self.data_info["fast_axis_end_pos"]=float(u[3][:-1])
                    self.data_info["fast_axis_num_step"]=int(u[4][:-1])
                    self.data_info["slow_axis"]=u[5][:-1]
                    self.data_info["slow_axis_start_pos"]=float(u[6][:-1])
                    self.data_info["slow_axis_end_pos"]=float(u[7][:-1])
                    self.data_info["slow_axis_num_step"]=int(u[8][:-1])
                    self.data_info["exposure_time"]=float(u[9][:])
                    self.shape = np.array((int(u[8][:-1]),
                                           int(u[4][:-1])))
                else:
                    raise ValueError("the type of scan format have not been included, please check")
                #dp = os.path.join(self.data_path,
                #                 "*%06d.h5" %(int(p[-1].split('.')[1])-1))
                dp = os.path.join(self.data_path,
                                 "*.h5")
                data_list = glob.glob(dp)
                if len(data_list) == 0:
                    print("No hdf data found")
                    return
                elif len(data_list) == 1:
                    if self.file_sys == 'linux':
                        self.data_h5 = data_list
                    elif self.file_sys == 'windows':
                        data_list[0] = data_list[0].replace("\\","/")
                        self.data_h5 = data_list
                elif len(data_list) > 1:
                    print("\nmore than one hdf data found")
                    #data_list = sorted(data_list,
                    #       key=lambda x:int(re.split("_",x)[-1][:-3]))
                    data_list = sorted(data_list)
                    if self.file_sys == 'linux':
                        self.data_h5 = data_list
                    elif self.file_sys == 'windows':
                        for _ in range(len(data_list)):
                            data_list[_] = data_list[_].replace("\\",'/')
                        self.data_h5 = data_list
                    #return
                print('\ndata information:\n{}'.format(self.data_info))
                print('\ndetector position:\n{}'.format(self.ndetx))
                #print('\nh5 data path:\n{}'.format(self.data_h5))
                if output_whole_list:
                    print('\nh5 data path:')
                    for lll in self.data_h5:
                        print("{}\n".format(lll))
                else:
                    print('\n%d h5 files in this scan' %len(self.data_h5))
                with h5py.File(self.data_h5[0],'r') as f:
                    self.single_h5_shape = f[self.data_h5path].shape
                print('\n\n')
            else:
                print("\n\nMore than one scan matched, please give more specific input")



if __name__ == '__main__':
    paras = []
    for _ in sys.argv[1:]:paras.append(_)
    t1 = id13_h5_search(paras[0],paras[1])
    if len(paras) > 2:
        t1.keyword_search(*paras[2:])
    else:
        t1.keyword_search()    
