import numpy as np
import math
import h5py
import os
import sys
import re
import glob
import warnings
warnings.filterwarnings("ignore")#ignore warnings

class scan_pos_info:
    def __init__(self,
                h5_fn,
                frame_num,
                ):
        self.fn = h5_fn
        self.frame_num = int(frame_num)
    
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
 
                ###############

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
                elif self.data_info["scan_type"]  == u'dmesh':
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
                #should revise ct34 to beam intensity and ndetx to detx
                try:
                    #print(self.data_info['scan_type'])
                    if 'kmap' in self.data_info['scan_type']:
                        self.mot1 = self.data_info['fast_axis']
                        self.mot2 = self.data_info['slow_axis']
                        self.ncols = self.data_info['fast_axis_num_step']
                        self.nrows = self.data_info['slow_axis_num_step']
                        self.pos1 = np.array(
                           f["{}/instrument/{}_position/value".format(self._data_name[0],self.mot1)])
                        self.pos2 = np.array(
                           f["{}/instrument/{}_position/value".format(self._data_name[0],self.mot2)])
                    elif 'dmesh' in self.data_info['scan_type']:
                        self.mot1 = self.data_info['axis1']
                        self.mot2 = self.data_info['axis2']
                        self.ncols = self.data_info['axis1_num_step']
                        self.nrows = self.data_info['axis2_num_step']
                        self.pos1 = np.array(
                           f["{}/instrument/{}_position/value".format(self._data_name[0],self.mot1)])
                        self.pos2 = np.array(
                           f["{}/instrument/{}_position/value".format(self._data_name[0],self.mot2)])
                    print('\n####################\n',
                          f'{self.mot1}, {self.mot2}, {self.pos1[self.frame_num]}, {self.pos2[self.frame_num]}',
                          '\n####################\n')
                           
                except Exception as e:
                    print(e,' or this is not a mesh scan')
            else:
                print("\n\nMore than one scan matched, please give more specific input")



if __name__ == '__main__':
    paras = []
    for _ in sys.argv[1:]:paras.append(_)
    t1 = scan_pos_info(paras[0],paras[-1])
    if len(paras) > 2:
        t1.keyword_search(*paras[1:-1])
    else:
        t1.keyword_search()    
