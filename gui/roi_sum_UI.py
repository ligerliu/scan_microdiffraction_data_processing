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

            
class roi_sum(QWidget):
    def __init__(self,obj):
        super().__init__()
        self.obj = obj
        self.roi_sum_UI()
        self.show()
        
    def roi_sum_UI(self):
        self.setGeometry(500,40,900,650)
        self.setWindowTitle('roi_sum')
        
        up_corner = QWidget(self)
        up_corner.setGeometry(20,20,200,20)
        up_corner_label = QLabel('up corner',up_corner)
        up_corner_label.setGeometry(0,0,80,20)
        up_corner_input = QLineEdit('1100,1100',up_corner)
        self.ul_corner_set(up_corner_input.text())
        up_corner_input.setGeometry(100,0,80,20)
        up_corner_input.textChanged[str].connect(self.ul_corner_set)
        
        dr_corner = QWidget(self)
        dr_corner.setGeometry(20,60,200,20)
        dr_corner_label = QLabel('down corner',dr_corner)
        dr_corner_label.setGeometry(0,0,80,20)
        dr_corner_input = QLineEdit('1300,1300',dr_corner)
        self.dr_corner_set(dr_corner_input.text())
        dr_corner_input.setGeometry(100,0,80,20)
        dr_corner_input.textChanged[str].connect(self.dr_corner_set)
        
        num_core_ui = QWidget(self)
        num_core_ui.setGeometry(20,100,200,20)
        num_core_label = QLabel('num cores',num_core_ui)
        num_core_label.setGeometry(0,0,80,20)
        self.num_core_input = QLineEdit('24',num_core_ui)
        self.num_core_input.setGeometry(100,0,80,20)

        proc_button = QPushButton('process',self)
        proc_button.setGeometry(20,130,80,20)
        proc_button.clicked.connect(self.roi_sum_process)
        
        pttn_roi_sum_button = QPushButton('roi sum pttn',self)
        pttn_roi_sum_button.setGeometry(20,160,80,20)
        pttn_roi_sum_button.clicked.connect(self.pttn_roi_sum_process)
        self.pttn_roi_sum_window = None
         
        self.bkgd_load_bttn = QPushButton('bkgd load',self)
        self.bkgd_load_bttn.move(20,220)
        self.bkgd_load_bttn.clicked.connect(self.load_bkgd)
        
        self.bkgd_sub_box = QCheckBox('bkgd subtract',self)
        self.bkgd_sub_box.move(20,250)
        
        pttn_vmin_input_ui = QWidget(self)
        pttn_vmin_input_ui.setGeometry(20,300,200,20)
        pttn_vmin_label = QLabel('pttn_vmin',pttn_vmin_input_ui)
        pttn_vmin_label.setGeometry(0,0,80,20)
        self.pttn_vmin_input = QLineEdit('0',pttn_vmin_input_ui)
        self.pttn_vmin_input.setGeometry(100,0,80,20)
        
        pttn_vmax_input_ui = QWidget(self)
        pttn_vmax_input_ui.setGeometry(20,330,200,20)
        pttn_vmax_label = QLabel('pttn_vmax',pttn_vmax_input_ui)
        pttn_vmax_label.setGeometry(0,0,80,20)
        self.pttn_vmax_input = QLineEdit('0.3',pttn_vmax_input_ui)
        self.pttn_vmax_input.setGeometry(100,0,80,20) 
        
        self.show()
         
    def load_bkgd(self):
        try:
            self.bkgd_fn,_ = QFileDialog.getOpenFileName(self,"open file","","")
            self.bkgd = fabio.open(self.bkgd_fn).data
        except Exception as e:
            print(e)
            pass
    
    def ul_corner_set(self,text):
        try:
            t = text.split(',')
            self.ul =  [int(t[0]),int(t[1])]
        except:
            pass
    
    def dr_corner_set(self,text):
        try:
            t = text.split(',')
            self.dr =  [int(t[0]),int(t[1])]
        except:
            pass
    
    def roi_sum_process(self):
        try:
            t = time()
            self.total_pttns,self.scan_shape,self.idx_list = \
                    scan_info(self.obj)
            self.h5_list,self.path_idx,self.pttn_idx = \
                    scan_h5_data_info(self.obj,
                               self.scan_shape,self.idx_list)
            data_path = self.obj.data_h5path
            num_core = int(self.num_core_input.text())
            res = parallel_func(scan_pttn_roi_sum,num_core,
                                np.arange(len(self.path_idx.flatten())),
                                h5_list   = self.h5_list,
                                path_idx  = self.path_idx.flatten(),
                                pttn_idx  = self.pttn_idx.flatten(),
                                data_path = data_path,
                                left_top  = self.ul,
                                right_bottom = self.dr,
                                )
            self.roi_map = np.zeros(self.scan_shape)
            for _ in range(self.scan_shape[0]):
                for __ in range(self.scan_shape[1]):
                    self.roi_map[_,__] = res[int(_*self.scan_shape[1]+__)]
            #self.roi_map = np.array(res).reshape(self.scan_shape)    
            self.roi_show = PlotWindow(self)#,position=True)
            colormap = setup_colormap(self.roi_map)
            self.roi_show.addImage(self.roi_map,colormap=colormap)

            toolBar = qt.QToolBar()
            self.roi_show.addToolBar(
            qt.Qt.BottomToolBarArea,
            toolBar)
            position = PositionInfo(plot= \
            self.roi_show,
            converters=[('X', 
            lambda x,y: int(np.round(x))),
            ('Y',
            lambda x,y: int(np.round(y)))])
            toolBar.addWidget(position)

            self.position = PositionInfo(plot=self.roi_show)
            self.roi_show.sigPlotSignal.connect(self.roi_map_clicked)
            self.roi_show.sigPlotSignal.connect(self.roi_map_polygon)
            self.roi_show.move(200,20)
            self.subwindow1 = None
            self.roi_show.show()
            print(time() - t)
        except Exception as e:
            print(str(e))
            pass
    
    def roi_map_clicked(self,event):
        widget = self.roi_show.getWidgetHandle()
        if 'Clicked' in event['event']:
            try:
                if isinstance(self.subwindow1,type(None)):
                    self.subwindow1 = Plot2D()#PlotWindow(position=True)
                else:
                    self.subwindow1.clear()
                    
                self.subwindow1.show()
                position = widget.mapFromGlobal(qt.QCursor.pos())
                xPixel,yPixel = position.x(),position.y()
                dataPos = self.roi_show.pixelToData(xPixel,yPixel,check=True)
                col,row = (int(dataPos[0]),int(dataPos[1]))
                with h5py.File(self.h5_list[self.path_idx[row,col]],'r') as f:
                    data = np.copy(f['entry_0000/measurement/data'][self.pttn_idx[row,col]])
                    data = data.astype(np.int)
                    if self.bkgd_sub_box.isChecked():
                        if not isinstance(self.bkgd,type(None)):
                            bkgd = np.copy(self.bkgd).astype(np.int)
                            data = data - bkgd
                    data[data<1] = 0 #data = data.astype(np.uint32)
                    #data = data.astype(np.uint32)
                    vmin = float(self.pttn_vmin_input.text())
                    vmax = float(self.pttn_vmax_input.text())
                    colormap = setup_colormap(data,vmin=vmin,vmax=vmax)
                    self.subwindow1.addImage(data,colormap=colormap)
                    self.subwindow1.setYAxisInverted(flag=True)
                    self.subwindow1.setKeepDataAspectRatio(True)                      
            except Exception as e:
                print(e)
        else:
            pass
    
    def roi_map_polygon(self,event):
        try:
            #self.mask_roi 
            mask = (self.roi_map*0+1)
            self.vertex = []
            if event['event'] == 'drawingFinished':
                coord = np.array(event['points'])
                coord = coord.astype(int)
                self.vertex.append(coord)
                self.mask_roi = mask_making(mask,self.vertex)
        except Exception as e:
            print(e)
            pass

    def pttn_roi_sum_process(self,high_thrshd=1e6):
        try:
            t = time()
            pttn_idx = self.pttn_idx[self.mask_roi].flatten()
            path_idx = self.path_idx[self.mask_roi].flatten()
            #mask_check = Plot2D()
            #mask_check.addImage(self.mask_roi)
            #mask_check.show()
            if self.bkgd_sub_box.isChecked():
                if not isinstance(self.bkgd,type(None)):
                    bkgd = np.copy(self.bkgd).astype(np.float)
                else:
                    bkgd = 0
            else:
                bkgd = 0
            for i,(p1,p2) in enumerate(zip(pttn_idx,path_idx)):
                with h5py.File(self.h5_list[p2],'r') as f:
                    data = np.copy(f['entry_0000/measurement/data'][p1]).astype(np.float)
                    data[data>1e6] = 0
                #print(bkgd)
                if i == 0:
                    sum_pttn = data - bkgd
                else:
                    sum_pttn += (data-bkgd)
            sum_pttn[sum_pttn<1] = 0
            sum_pttn /= (i+1)
            if isinstance(self.pttn_roi_sum_window,type(None)): 
                self.pttn_roi_sum_window = Plot2D()
            vmin = float(self.pttn_vmin_input.text())
            vmax = float(self.pttn_vmax_input.text())
            colormap = setup_colormap(data,vmin=vmin,vmax=vmax)
            self.pttn_roi_sum_window.addImage(sum_pttn,colormap=colormap)
            self.pttn_roi_sum_window.setYAxisInverted(flag=True)
            self.pttn_roi_sum_window.setKeepDataAspectRatio(True)                      
            self.pttn_roi_sum_window.show()
            print(time()-t)
        except Exception as e:
            print(e)
            pass 
        
        
