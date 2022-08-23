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

def sum_single_h5_eiger(num,h5_list,idx_list,
                        single_h5_data_size,data_path):
    #print(num,h5_list[num],idx_list)
    start = int(idx_list[num][0] - num*single_h5_data_size)
    end = int(idx_list[num][1] - num*single_h5_data_size)
    with h5py.File(h5_list[num],'r') as f:
        data = f[data_path][start:end,:]
    ave_data = np.nanmean(data,axis=0)
    max_data = np.nanmax(data,axis=0)
    return ave_data,max_data
            
class total_sum(QWidget):
    # average and maximum project all the pattern
    def __init__(self,obj):
        super().__init__()
        self.obj = obj
        self.total_sum_UI()
        self.show()
        
    def total_sum_UI(self):
        self.setGeometry(500,40,900,650)
        self.setWindowTitle('total_sum')
        layout = QGridLayout()
        
        panel_layout = QWidget() 
        num_core_ui = QWidget(panel_layout)
        num_core_ui.setGeometry(20,100,200,20)
        num_core_label = QLabel('num cores',num_core_ui)
        num_core_label.setGeometry(0,0,80,20)
        self.num_core_input = QLineEdit('12',num_core_ui)
        self.num_core_input.setGeometry(100,0,80,20)

        proc_button = QPushButton('process',panel_layout)
        proc_button.setGeometry(20,130,80,20)
        proc_button.clicked.connect(self.total_sum_proc)
        
        pttn_vmin_input_ui = QWidget(panel_layout)
        pttn_vmin_input_ui.setGeometry(20,300,200,20)
        pttn_vmin_label = QLabel('pttn_vmin',pttn_vmin_input_ui)
        pttn_vmin_label.setGeometry(0,0,80,20)
        self.pttn_vmin_input = QLineEdit('0',pttn_vmin_input_ui)
        self.pttn_vmin_input.setGeometry(100,0,80,20)
        
        pttn_vmax_input_ui = QWidget(panel_layout)
        pttn_vmax_input_ui.setGeometry(20,330,200,20)
        pttn_vmax_label = QLabel('pttn_vmax',pttn_vmax_input_ui)
        pttn_vmax_label.setGeometry(0,0,80,20)
        self.pttn_vmax_input = QLineEdit('1000',pttn_vmax_input_ui)
        self.pttn_vmax_input.setGeometry(100,0,80,20) 
        
        self.ave_img_widget = QWidget(self)
        self.max_img_widget = QWidget(self)
        
        layout.addWidget(panel_layout,0,0,10,4)
        self.setLayout(layout)
        
        self.img_type = QComboBox()
        self.img_type.addItems(["ave img","max img"])
        self.img_type.activated.connect(self.switchPage)
        layout.addWidget(self.img_type,0,5,1,10)
        self.stack_img =QStackedWidget()
        self.stack_img.addWidget(self.ave_img_widget)
        self.stack_img.addWidget(self.max_img_widget)
        layout.addWidget(self.stack_img,1,5,9,10)
        #img_layout.addWidget(self.ave_img_widget)
        #layout.addLayout(img_layout)
        
        #self.setLayout(layout)
        self.show()
    
    def switchPage(self):
        self.stack_img.setCurrentIndex(self.img_type.currentIndex())
         
    def total_sum_proc(self):
        try:
            if hasattr(self.obj,'data_h5'):
                t = time()
            
                self.total_pttns,self.scan_shape,self.idx_list = \
                        scan_info(self.obj)
                self.h5_list,self.path_idx,self.pttn_idx = \
                        scan_h5_data_info(self.obj,
                                   self.scan_shape,self.idx_list)
                data_path = self.obj.data_h5path
                num_cores = int(self.num_core_input.text())
                    
                with h5py.File(self.h5_list[0],'r') as f:
                    single_h5_num_pttn =  f[data_path].shape[0]
                    data_shape         =  f[data_path][0].shape
                idx_list = chunk_idx(self.total_pttns,single_h5_num_pttn) 
                res = parallel_func(sum_single_h5_eiger,num_cores,
                        np.arange(len(self.h5_list)),h5_list=self.h5_list,
                        idx_list=idx_list,
                        single_h5_data_size=single_h5_num_pttn,
                        data_path=data_path)
                ave_data = []
                max_data = []
                for _ in res:
                    ave_data.append(_[0])
                    max_data.append(_[1])
                ave_data = np.array(ave_data)
                ave_data = np.nanmean(ave_data,axis=0).astype(float)
                self.ave_img = np.copy(ave_data)
                max_data = np.array(max_data)
                max_data = np.nanmax(max_data,axis=0).astype(np.uint16)
                self.max_img = np.copy(max_data)
                img1 = fabio.edfimage.EdfImage(ave_data)
                img2 = fabio.edfimage.EdfImage(max_data)
                save_path = QFileDialog.getExistingDirectory()
                save_folder = os.path.join(save_path,'edf_sum_file')
                #os.chdir(save_path)
                if not os.path.exists(save_folder):
                    os.mkdir(save_folder)
                save_name_ave = os.path.join(save_folder,
                              "{}_ave.edf".format(self.obj._data_name[0][:-4]))
                save_name_max = os.path.join(save_folder,
                              "{}_max.edf".format(self.obj._data_name[0][:-4]))
                #os.chdir("edf_sum_file")
                img1.write(save_name_ave)
                img2.write(save_name_max)
                self.img_show(self.ave_img,self.ave_img_widget)
                self.img_show(self.max_img,self.max_img_widget)
                super().__init__()
                print("\nsum processing time:\n%5.2f sec" %(time()-t))
        except Exception as e:
            print(e)
            pass
 
    def img_show(self,img,widget):
        try:
            self.roi_show = PlotWindow(widget)#,position=True)
            self.pttn_vmin = float(self.pttn_vmin_input.text())
            self.pttn_vmax = float(self.pttn_vmax_input.text())
            colormap = setup_colormap(img,
                                      vmin=self.pttn_vmin,
                                      vmax=self.pttn_vmax)
            self.roi_show.addImage(img,colormap=colormap)
            self.roi_show.setKeepDataAspectRatio(flag=True)
            self.roi_show.show()
            
        except Exception as e:
            print(str(e))
            pass
    
        
        
