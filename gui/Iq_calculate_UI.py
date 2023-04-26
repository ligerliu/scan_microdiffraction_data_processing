import sys
sys.path.append('../xs_proc')
from h5_data_search import *
from xs_data_proc import *
from proc_data_ana import *     
from visual_func import * 
from gui_scan_xrd_analysis import *
#from multi_scan_qphi_proc import * 
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
    
class Iq_calculate(QWidget):
    def __init__(self,obj):#h5_list,path_idx,pttn_idx,scan_shape,obj):
        super().__init__()
        self.obj = obj
        self.poni  = None
        self.mask_file = None
        self.Iq_cal_UI()
        self.show()
    
    def Iq_cal_UI(self):
        self.setGeometry(500,40,900,650)
        self.setWindowTitle('Automatic 1D integrate')
        
        poni_load = QWidget(self)
        poni_load.setGeometry(20,20,240,80)
        poni_file_load = QPushButton('load poni',poni_load)
        poni_file_load.clicked.connect(self.load_poni_window)
        
        mask_load = QWidget(self)
        mask_load.setGeometry(20,120,240,80)
        mask_file_load = QPushButton('load mask',mask_load)
        mask_file_load.clicked.connect(self.load_mask_window)
        
        Iq_shape = QWidget(self)
        Iq_shape.setGeometry(20,240,240,120)
        Iq_q_label = QLabel("q_npts:",Iq_shape)
        Iq_q_label.setGeometry(10,10,60,20)
        self.Iq_q_shape = QLineEdit("720",Iq_shape)
        self.Iq_q_shape.setGeometry(80,10,120,20)
        core_num_label = QLabel("num_core:",Iq_shape)
        core_num_label.setGeometry(10,40,60,20)
        self.core_num = QLineEdit("12",Iq_shape)
        self.core_num.setGeometry(80,40,120,20)
        #data_path_label = QLabel("data_path:",Iq_shape)
        #data_path_label.setGeometry(10,70,60,20)
        #self.data_path_input = QLineEdit("entry_0000/measurement/data",Iq_shape)
        #self.data_path_input.setGeometry(80,70,120,20)

        save_path_load = QPushButton('save path load',Iq_shape)
        save_path_load.setGeometry(10,70,120,20)
        save_path_load.clicked.connect(self.load_path)
                
        save_path_label = QLabel("save_path:",Iq_shape)
        save_path_label.setGeometry(10,100,60,20)
        self.save_path_input = QLineEdit("",Iq_shape)
        self.save_path_input.setGeometry(80,100,120,20)
        
        process_button = QPushButton("process",self)
        process_button.setGeometry(20,420,60,20)
        process_button.clicked.connect(self.Iq_cal)
        
        load_Iq_button = QPushButton("load_Iq",self)
        load_Iq_button.setGeometry(20,450,60,20)
        load_Iq_button.clicked.connect(self.load_Iq)
        
        Iq_2d_map_button = QPushButton("Iq map",self)
        Iq_2d_map_button.setGeometry(20,480,60,20)
        Iq_2d_map_button.clicked.connect(self.Iq_map_show)
    
        self.add_map_box = QCheckBox('add Iq map',self)
        self.add_map_box.move(20,510)
        self.map_id = 0
        
    
    def load_path(self):
        try:
            path = QFileDialog.getExistingDirectory()
            self.save_path_input.setText(path)
        except:
            pass
        
    def load_Iq(self):
        try:
            self.Iq_h5,_ = QFileDialog.getOpenFileName(self,"load Iq","","")
            self.q  = load_proc_dataset(self.Iq_h5,'q',proc_type="integrate1d")
            num_core = int(self.core_num.text())
            self.Iq = Iq_serries_2dmap(self.Iq_h5,num_cores=num_core)
        except Exception as e:
            print(e)
            pass
        
    def Iq_map_show(self):
        Iq_show = Plot2D(self)#PlotWindow(self)#,position=True)
        Iq_show.setGeometry(240,20,640,480)
        data = np.copy(self.Iq)
        #print(data.shape)
        #data = data.reshape((self.Iq.shape[0]*self.Iq.shape[1],
        #                     self.Iq.shape[2]))
        if self.add_map_box.isChecked():
            self.img = np.vstack((data,self.img))
        else:
            self.img = np.copy(data) 
        colormap = setup_colormap(self.img)
        Iq_show.addImage(self.img,
                        colormap=colormap,
                        scale = ((self.q[-1]-self.q[0])/len(self.q),1),
                        ) 
        Iq_show.setGraphXLabel(label=r'$Q\,\,(\AA^{-1})$')
        toolBar = qt.QToolBar()
        Iq_show.addToolBar(
        qt.Qt.BottomToolBarArea,
        toolBar)
        position = PositionInfo(plot= \
        Iq_show,
        converters=[('X', 
        lambda x,y: x),
        ('Y',
        lambda x,y: y),
        ('value',
        lambda x,y: self.img[int(np.round(y)),
                        np.argmin(np.abs(self.q-x))]),
        ])
        toolBar.addWidget(position)
        Iq_show.show()
        
    def load_mask_window(self):
        try:
            self.mask_file,_ = QFileDialog.getOpenFileName(self,"load mask","","")
            if '.npz' in self.mask_file:
                self.mask = np.load(self.mask_file)['mask']
            elif '.edf' in self.mask_file:
                self.mask = (fabio.open(self.mask_file).data).astype(bool)
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
    
    def Iq_cal(self,**kwargs):
        # currently still in form of ID13 customization, need to be generized
        try:
            t = time()
            self.total_pttns,self.scan_shape,self.idx_list = \
                    scan_info(self.obj)
            self.h5_list,self.path_idx,self.pttn_idx = \
                    scan_h5_data_info(self.obj,
                               self.scan_shape,self.idx_list)
            print(self.h5_list[0],self.obj._data_name)
            self.q_npts = int(self.Iq_q_shape.text())
            self.data_path = self.obj.data_h5path#str(self.data_path_input.text())
            num_core = int(self.core_num.text())
            
            if not isinstance(self.obj.ct34,type(None)):
                self.ct34 = np.copy(self.obj.ct34)
                cur_pos   = np.argwhere(self.ct34<(np.nanmean(self.ct34)*0.6))
                self.ct34[cur_pos] = self.ct34[cur_pos-1]
                if np.nanmean(self.ct34) <= 0:
                    self.ct34 = np.ones((len(self.pttn_idx.flatten()),))
                else:
                    self.ct34 /= np.nanmean(self.ct34)
            else:
                    self.ct34 = np.ones((len(self.pttn_idx.flatten()),))
            #print(self.ct34,len(self.ct34))
            
            res = parallel_func(scan_calculate_Iq,
                           num_core,
                           np.arange(self.total_pttns),
                           h5_list   = self.h5_list,
                           path_idx  = self.path_idx.flatten(),
                           pttn_idx  = self.pttn_idx.flatten(),
                           data_path = self.data_path,
                           pyfai_obj = self.ai,
                           mask      = self.mask,
                           q_npts    = self.q_npts,
                           ct        = self.ct34,
                           **kwargs
                           )
            q    = res[0][0]
            q   /= 10
            scan_shape = self.scan_shape
            #Iq   = np.zeros((scan_shape[0],
            #                 scan_shape[1],
            #                 len(res[0][1]),
            #                ))
            #for _ in range(len(res)):
            #    #if len(idx_list) > 1:
            #    #    i1 = int((_+i*t1.single_h5_shape[0])/scan_shape[1])
            #    #    i2 = int((_+i*t1.single_h5_shape[0])%scan_shape[1])
            #    #else:
            #    i1 = int(_/scan_shape[1])
            #    i2 = int(_%scan_shape[1])
            #    Iq[i1,i2,:]   = res[_][1]
            name = self.obj._data_name[0].split('.')[0]
            save_path =self.save_path_input.text()
            print(time()-t)
            save_Iq_as_h5(self.obj,save_path,name,
                    q    = q,
                    res   = res,
                    path_idx = self.path_idx,
                    pttn_idx = self.pttn_idx,
                    total_pttn_num=self.total_pttns,
                    single_h5_pttn_num = self.obj.single_h5_shape[0],
                    h5_path_list = self.h5_list)

            print(time()-t)#,'\n\nfuck complete')
            self.q  = q
            #self.Iq = Iq 
        except Exception as e:
            print(e)
            pass
