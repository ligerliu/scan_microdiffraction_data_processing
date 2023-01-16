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
    
class Iq_analysis(QWidget):
    # slow reduntantly load proc data feels like
    def __init__(self):#,*args,**kwargs):
        super().__init__()
        self.Iq_ana_UI()
        self.fn = None
        self.bkgd = None
        #self.vertex = []
        self.show()
    
    def Iq_ana_UI(self):
        self.setGeometry(500,40,900,650)
        self.setWindowTitle('Iq analysis')
        
        load_file_widget = QWidget(self)
        load_file_widget.setGeometry(20,20,180,100)
        load_file_widget.setWindowTitle('load proc hdf')
        fileDialogBttn = QPushButton('open proc',load_file_widget)
        fileDialogBttn.setGeometry(10,20,80,20)
        fileDialogBttn.clicked.connect(self.jl_Iq_load)
        
        self.fileInputWin = QLineEdit('',load_file_widget)
        self.fileInputWin.setGeometry(10,50,160,20)
        self.fileInputWin.textEdited.connect(self.line_edit_input)
        
        roi_setup    = QWidget(self)
        roi_setup.setGeometry(10,120,180,320)
        qmin_label = QLabel('qmin',roi_setup)
        qmin_label.setGeometry(10,10,60,20)
        self.qmin = QLineEdit('',roi_setup)
        self.qmin.setAlignment(Qt.AlignRight)
        self.qmin.setGeometry(80,10,100,20)
        qmax_label = QLabel('qmax',roi_setup)
        qmax_label.setGeometry(10,40,60,20)
        self.qmax = QLineEdit('',roi_setup)
        self.qmax.setAlignment(Qt.AlignRight)
        self.qmax.setGeometry(80,40,100,20)
        
        self.low_int_thrhd_label = QLabel('low thrhd',roi_setup)
        self.low_int_thrhd_label.setGeometry(10,130,60,20)
        self.low_int_thrhd_line = QLineEdit('0',roi_setup)
        self.low_int_thrhd_line.setGeometry(80,130,100,20)
        self.high_int_thrhd_label = QLabel('high thrhd',roi_setup)
        self.high_int_thrhd_label.setGeometry(10,160,60,20)
        self.high_int_thrhd_line = QLineEdit('1e6',roi_setup)
        self.high_int_thrhd_line.setGeometry(80,160,100,20)
        
        # this is function generate intensity of map of certain q range
        self.roi_sum_bttn = QPushButton('roi sum',roi_setup)
        self.roi_sum_bttn.move(10,190)
        self.roi_sum_bttn.clicked.connect(self.Iq_roi_sum)
        
        # this is ave of roi within intensity map generated above
        self.Iq_roi_ave_bttn = QPushButton('Iq roi ave',roi_setup)
        self.Iq_roi_ave_bttn.move(100,190)
        self.Iq_roi_ave_bttn.clicked.connect(self.Iq_roi_ave)
        self.ave_window = None
        
        # need to check the txt loading function        
        self.bkgd_load_bttn = QPushButton('bkgd load',roi_setup)
        self.bkgd_load_bttn.move(10,220)
        self.bkgd_load_bttn.clicked.connect(self.load_bkgd)
        
        self.bkgd_sub_box = QCheckBox('bkgd subtract',roi_setup)
        self.bkgd_sub_box.move(10,250)
        
        self.keep_current_Iq = QCheckBox('keep curve',roi_setup)
        self.keep_current_Iq.move(10,280)
        
        self.show() 
    
    
    def load_bkgd(self):
        try:
            self.bkgd_fn,_ = QFileDialog.getOpenFileName(self,"open file","","")
            bkgd_ = np.loadtxt(self.bkgd_fn)
            self.bkgd = bkgd_[:,1]
        except:
            pass
             
    def Iq_roi_sum(self):
        try:
            qmin = float(self.qmin.text())
            qmax = float(self.qmax.text())
            low_thrd = float(self.low_int_thrhd_line.text())
            high_thrd = self.high_int_thrhd_line.text()
            if high_thrd != '':
                high_thrd = float(high_thrd)
            else:
                high_thrd = np.inf
                #Iq[Iq>high_thrd] = np.nan
            if self.bkgd_sub_box.isChecked():
                self.roi = Iq_roi_map_2dmap(self.fn,q=self.q,
                                qmin=qmin,qmax=qmax,
                                low_thrd=low_thrd,high_thrd=high_thrd,
                                bkgd=self.bkgd)
            else:
                self.roi = Iq_roi_map_2dmap(self.fn,q=self.q,
                                qmin=qmin,qmax=qmax,
                                low_thrd=low_thrd,high_thrd=high_thrd)
            roi = np.copy(self.roi)
            
            roi_map = QWidget(self)
            self.Iq_roi_map = PlotWindow(roi_map)#,position=True)
            roi_map.setGeometry(220,20,680,640)
            colormap = setup_colormap(roi)
            self.Iq_roi_map.addImage(roi,colormap=colormap,replace=True)
            self.Iq_roi_map.setYAxisInverted(flag=True)
            toolBar = qt.QToolBar()
            self.Iq_roi_map.addToolBar(
            qt.Qt.BottomToolBarArea,
            toolBar)
            position = PositionInfo(plot= \
            self.Iq_roi_map,
            converters=[('X', 
            lambda x,y: int(np.round(x))),
            ('Y',
            lambda x,y: int(np.round(y))),
            ('value',
            lambda x,y: roi[int(np.round(y)),int(np.round(x))]),
            ])
            toolBar.addWidget(position)
            
            self.subwindow = None
            self.Iq_roi_map.sigPlotSignal.connect(self.roi_map_clicked)
            self.Iq_roi_map.sigPlotSignal.connect(self.roi_map_polygon)
            roi_map.show()
            #self.subwindow.show()
         
        except Exception as e:
            print(e)
            pass
                
        
    def load_proc_hdf(self):    
        if self.fn != '':
            try:
                #self.Iq_h5,_ = QFileDialog.getOpenFileName(self,"load Iq","","")
                self.q = load_proc_dataset(self.fn,'q',proc_type='integrate1d')
                self.qmin.setText(str(np.round(np.min(self.q),decimals=4)))
                self.qmax.setText(str(np.round(np.max(self.q),decimals=4)))
            except Exception as e:
                print(e)
        else:
            print(f"data path: {self.fn} is not correct")
            pass
         
    def line_edit_input(self):
        fn_line_input = self.fileInputWin.text()
        if isinstance(self.fn,type(None)):
            self.fn = file_line_input  
            self.load_proc_hdf()
        elif self.fn != fn_line_input:
            self.fn = file_line_input  
            self.load_proc_hdf()
        
    def jl_Iq_load(self):
        self.fn,_ = QFileDialog.getOpenFileName(self,"open file","","")
        self.fileInputWin.setText(self.fn)
        self.load_proc_hdf()
    
    def roi_map_polygon(self,event):          
        try:
            #self.mask_roi = self.Iq_roi_map.getSelectionMask()         
            self.vertex = []
            mask = (self.roi*0+1).astype(bool)
            if event['event'] == 'drawingFinished':
                coord = np.array(event['points'])
                coord = coord.astype(np.int)
                self.vertex.append(coord)
                #print(vertex)
                self.mask_roi = mask_making(mask,self.vertex)
        except Exception as e:
                print(e)
                pass
    
    def Iq_roi_ave(self):
        try:
            self.mask_roi = self.mask_roi.astype(bool)
            #print(qphi.shape,qphi.dtype,self.mask_roi.shape,self.mask_roi.dtype)
            scan_shape = load_proc_Iq_size(self.fn,proc_type='integrate1d')
            rr,cc = np.mgrid[0:scan_shape[0],0:scan_shape[1]]
            r = rr[self.mask_roi]
            c = cc[self.mask_roi]
            Iq = []
            for _ in range(len(r)):
                Iq.append(load_proc_single_Iq(self.fn,r[_],c[_]))
            Iq = np.array(Iq)
            Iq_ave = np.nanmean(Iq,axis=0)
            if self.bkgd_sub_box.isChecked():
                #print(qphi_ave.dtype,self.bkgd.dtype)
                if not isinstance(self.bkgd,type(None)):
                    Iq_ave = Iq_ave - self.bkgd
            
            if isinstance(self.subwindow,type(None)):
                self.subwindow = Plot1D()
                self.Iq_old = []
            else:
                if not self.keep_current_Iq.isChecked():
                    self.subwindow.clear()
                    self.Iq_old = []
            if self.keep_current_Iq.isChecked():
                self.Iq_old.append(Iq_ave)
                for i in range(len(self.Iq_old)):
                    self.subwindow.addCurve(self.q,self.Iq_old[i],
                            legend=f"{i}",xlabel='Q',ylabel='I (a.u.)')
            else:
                self.Iq_old = []
                self.subwindow.addCurve(self.q,Iq_ave,
                xlabel='Q',ylabel='I (a.u.)')
            self.subwindow.show()
        except Exception as e:
            print(e)
            pass
     
    def roi_map_clicked(self,event):
        widget = self.Iq_roi_map.getWidgetHandle()
        if event['event'] == 'mouseClicked':
            position = widget.mapFromGlobal(qt.QCursor.pos())
            xPixel,yPixel = position.x(),position.y()
            dataPos = self.Iq_roi_map.pixelToData(xPixel,yPixel,check=True)
            try:
                col,row = (int(dataPos[0]),int(dataPos[1]))
                print('\n pttn position: x-{}, y-{}'.format(col,row))
                
                if self.bkgd_sub_box.isChecked():
                    if isinstance(self.bkgd,type(None)):
                        #data = np.copy(self.qphi[row,col])
                        data = np.copy(load_proc_single_Iq(self.fn,row,col,
                                            proc_type="integrate1d"))
                    else:
                        bkgd = np.copy(self.bkgd).astype(np.float)
                        data = np.copy(load_proc_single_Iq(self.fn,row,col,
                                            proc_type="integrate1d"))
                        data = data - bkgd
                else:
                    data = np.copy(load_proc_single_Iq(self.fn,row,col,
                                            proc_type="integrate1d"))
                data[np.isnan(data)] = 0
                data[np.isinf(data)] = 0
            
                if isinstance(self.subwindow,type(None)):
                    self.subwindow = Plot1D()
                    self.Iq_old = []
                else:
                    if not self.keep_current_Iq.isChecked():
                        self.subwindow.clear()
                        self.Iq_old = []
                if self.keep_current_Iq.isChecked():
                    self.Iq_old.append(data)
                    for i in range(len(self.Iq_old)):
                        self.subwindow.addCurve(self.q,self.Iq_old[i],
                                legend=f"{i}",xlabel='Q',ylabel='I (a.u.)')
                else:
                    self.Iq_old = []
                    self.subwindow.addCurve(self.q,data,legend='0',
                                xlabel='Q',ylabel='I (a.u.)')
                self.subwindow.show()
            except Exception as e:
                print(e)
        else:
            pass
        
        
