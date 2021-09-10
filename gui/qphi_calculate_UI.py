import sys
sys.path.append('../xs_proc')
from h5_data_search import *
from xs_data_proc import *
from proc_data_ana import *     
from visual_func import * 
from multi_scan_qphi_proc import * 
from multi_scan_Iq_proc import *
 
import os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import h5py
import numpy as np
import fabio

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import file_tree as ft
from time import time

from silx.gui.plot import *
from silx.gui.plot.PlotWidget import PlotWidget
from silx.gui import qt
from silx.gui.plot.PlotWindow import PlotWindow 
from silx.gui.plot.PlotTools import PositionInfo
from silx.gui.plot.ImageView import ImageView
from silx.gui.colors import Colormap

from colormap_setup import setup_colormap                
    
class qphi_calculate(QWidget):
    def __init__(self,obj):#h5_list,path_idx,pttn_idx,scan_shape,obj):
        super().__init__()
        #self.h5_list  = h5_list
        #self.path_idx = path_idx
        #self.pttn_idx = pttn_idx
        #self.scan_shape = scan_shape
        self.poni  = None
        self.mask_file = None
        self.obj = obj
        self.qphi_cal_UI()
        self.show()
    
    def qphi_cal_UI(self):
        self.setGeometry(500,40,900,650)
        self.setWindowTitle('Automatic 2D integrate')
        
        poni_load = QWidget(self)
        poni_load.setGeometry(20,20,240,80)
        poni_file_load = QPushButton('load poni',poni_load)
        poni_file_load.clicked.connect(self.load_poni_window)
        
        mask_load = QWidget(self)
        mask_load.setGeometry(20,60,240,80)
        mask_file_load = QPushButton('load mask',mask_load)
        mask_file_load.clicked.connect(self.load_mask_window)
        #if isinstance(self.mask_file,type(None)):
        #    print("mask is needed, at current develop version")
        #    return
        #else:
        #    self.mask = np.load(self.mask_file)['mask']
        
        qphi_shape = QWidget(self)
        qphi_shape.setGeometry(20,120,240,160)
        qphi_q_label = QLabel("q_npts:",qphi_shape)
        qphi_q_label.setGeometry(10,10,60,20)
        self.qphi_q_shape = QLineEdit("720",qphi_shape)
        self.qphi_q_shape.setGeometry(80,10,120,20)
        qphi_a_label = QLabel("a_npts:",qphi_shape)
        qphi_a_label.setGeometry(10,40,60,20)
        self.qphi_a_shape = QLineEdit("360",qphi_shape)
        self.qphi_a_shape.setGeometry(80,40,120,20)
        core_num_label = QLabel("num_core:",qphi_shape)
        core_num_label.setGeometry(10,70,60,20)
        self.core_num = QLineEdit("12",qphi_shape)
        self.core_num.setGeometry(80,70,120,20)
        #data_path_label = QLabel("data_path:",qphi_shape)
        #data_path_label.setGeometry(10,100,60,20)
        #self.data_path_input = QLineEdit("entry_0000/measurement/data",qphi_shape)
        #self.data_path_input.setGeometry(80,100,120,20)
        
        save_path_label = QLabel("save_path:",qphi_shape)
        save_path_label.setGeometry(10,130,60,20)
        self.save_path_input = QLineEdit("/data/id13/inhouse12/jiliang/code_v3/example/sc5005/",qphi_shape)
        self.save_path_input.setGeometry(80,130,120,20)
        
        process_button = QPushButton("process",self)
        process_button.setGeometry(20,450,60,20)
        process_button.clicked.connect(self.qphi_cal)
        
    
    def load_mask_window(self):
        try:
            self.mask_file,_ = QFileDialog.getOpenFileName(self,"load mask","","")
            self.mask = np.load(self.mask_file)['mask']
        except Exception as e:
            print(e)
            pass
    
    def load_poni_window(self):
        try:
            self.poni,_ = QFileDialog.getOpenFileName(self,"load poni","","")
            self.ai = pyFAI.load(self.poni)
        except Exception as e:
            print(e)
            pass
    
    def qphi_cal(self,**kwargs):
        # currently still in form of ID13 customization, need to be generized
        try:
            t = time()
            self.total_pttns,self.scan_shape,self.idx_list = \
                    scan_info(self.obj)
            self.h5_list,self.path_idx,self.pttn_idx = \
                    scan_h5_data_info(self.obj,
                               self.scan_shape,self.idx_list)
            print(self.h5_list[0],self.obj._data_name)
            self.q_npts = int(self.qphi_q_shape.text())
            self.a_npts = int(self.qphi_a_shape.text())
            self.data_path = self.obj.data_h5path#str(self.data_path_input.text())
            num_core = int(self.core_num.text())
            
            res = parallel_func(scan_calculate_Iqphi,
                           num_core,
                           np.arange(len(self.path_idx.flatten())),
                           h5_list   = self.h5_list,
                           path_idx  = self.path_idx.flatten(),
                           pttn_idx  = self.pttn_idx.flatten(),
                           data_path = self.data_path,
                           pyfai_obj = self.ai,
                           mask      = self.mask,
                           q_npts    = self.q_npts,
                           a_npts    = self.a_npts,
                           **kwargs
                           )
            q    = res[0][1]
            azi  = res[0][2]
            scan_shape = self.scan_shape
            qphi = np.zeros((scan_shape[0],
                             scan_shape[1],
                             res[0][0].shape[0],
                             res[0][0].shape[1]))
            for _ in range(len(res)):
                #if len(idx_list) > 1:
                #    i1 = int((_+i*t1.single_h5_shape[0])/scan_shape[1])
                #    i2 = int((_+i*t1.single_h5_shape[0])%scan_shape[1])
                #else:
                i1 = int(_/scan_shape[1])
                i2 = int(_%scan_shape[1])
                qphi[i1,i2,:]   = res[_][0]
            name = self.obj._data_name[0].split('.')[0]
            save_path =self.save_path_input.text()
            print(time()-t)
            save_qphi_as_h5(self.obj,save_path,name,
                    q    = q,
                    azi  = azi,
                    qphi = qphi,
                    path_idx = self.path_idx,
                    pttn_idx = self.pttn_idx,
                    h5_path_list = self.h5_list)
            print(time()-t)#,'\n\nfuck complete')
        except Exception as e:
            print(e)
            pass
