import sys
sys.path.append('../xs_proc')
from h5_data_search import *
from xs_data_proc import *
from proc_data_ana import *     
from visual_func import * 
from gui_scan_xrd_analysis import *
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
    
class qphi_analysis(QWidget):
    # slow reduntantly load proc data feels like
    def __init__(self):#,*args,**kwargs):
        super().__init__()
        self.qphi_ana_UI()
        self.fn = None
        self.bkgd = None
        #self.vertex = []
        self.show()
    
    def qphi_ana_UI(self):
        self.setGeometry(500,40,900,650)
        self.setWindowTitle('qphi analysis')
        
        load_file_widget = QWidget(self)
        load_file_widget.setGeometry(20,20,180,100)
        load_file_widget.setWindowTitle('load proc hdf')
        fileDialogBttn = QPushButton('open proc',load_file_widget)
        fileDialogBttn.setGeometry(10,20,80,20)
        fileDialogBttn.clicked.connect(self.jl_qphi_load)
        
        self.fileInputWin = QLineEdit('',load_file_widget)
        self.fileInputWin.setGeometry(10,50,160,20)
        self.fileInputWin.textEdited.connect(self.line_edit_input)
        
        roi_setup    = QWidget(self)
        roi_setup.setGeometry(10,120,180,420)
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
        amin_label = QLabel('amin',roi_setup)
        amin_label.setGeometry(10,70,60,20)
        self.amin = QLineEdit('',roi_setup)
        self.amin.setAlignment(Qt.AlignLeft)
        self.amin.setGeometry(80,70,100,20)
        amax_label = QLabel('amax',roi_setup)
        amax_label.setGeometry(10,100,60,20)
        self.amax = QLineEdit('',roi_setup)
        self.amax.setAlignment(Qt.AlignLeft)
        self.amax.setGeometry(80,100,100,20)
        self.low_int_thrhd_label = QLabel('low thrhd',roi_setup)
        self.low_int_thrhd_label.setGeometry(10,130,60,20)
        self.low_int_thrhd_line = QLineEdit('0',roi_setup)
        self.low_int_thrhd_line.setGeometry(80,130,100,20)
        self.high_int_thrhd_label = QLabel('high thrhd',roi_setup)
        self.high_int_thrhd_label.setGeometry(10,160,60,20)
        self.high_int_thrhd_line = QLineEdit('1e6',roi_setup)
        self.high_int_thrhd_line.setGeometry(80,160,100,20)
        
        
        self.roi_sum_bttn = QPushButton('roi sum',roi_setup)
        self.roi_sum_bttn.move(10,190)
        self.roi_sum_bttn.clicked.connect(self.qphi_roi_sum)
        
        self.qphi_roi_ave_bttn = QPushButton('qphi roi ave',roi_setup)
        self.qphi_roi_ave_bttn.move(100,190)
        self.qphi_roi_ave_bttn.clicked.connect(self.qphi_roi_ave)
        self.ave_window = None
        
        self.bkgd_load_bttn = QPushButton('bkgd load',roi_setup)
        self.bkgd_load_bttn.move(10,220)
        self.bkgd_load_bttn.clicked.connect(self.load_bkgd)
        
        self.bkgd_sub_box = QCheckBox('bkgd subtract',roi_setup)
        self.bkgd_sub_box.move(10,250)
        
        pttn_vmin_input_ui = QWidget(roi_setup)
        pttn_vmin_input_ui.setGeometry(20,310,200,20)
        pttn_vmin_label = QLabel('pttn_vmin',pttn_vmin_input_ui)
        pttn_vmin_label.setGeometry(0,0,80,20)
        self.pttn_vmin_input = QLineEdit('0',pttn_vmin_input_ui)
        self.pttn_vmin_input.setGeometry(100,0,80,20)
        
        pttn_vmax_input_ui = QWidget(roi_setup)
        pttn_vmax_input_ui.setGeometry(20,330,200,20)
        pttn_vmax_label = QLabel('pttn_vmax',pttn_vmax_input_ui)
        pttn_vmax_label.setGeometry(0,0,80,20)
        self.pttn_vmax_input = QLineEdit('1',pttn_vmax_input_ui)
        self.pttn_vmax_input.setGeometry(100,0,80,20) 
        
        self.pttn_log_box = QCheckBox('pttn log',roi_setup)
        self.pttn_log_box.move(10,360)
        
        self.pttn_zoom_box = QCheckBox('pttn zoom',roi_setup)
        self.pttn_zoom_box.move(10,390)
        
        self.ori_map_bttn = QPushButton('ori map',roi_setup)
        self.ori_map_bttn.move(10,280)
        self.ori_map_bttn.clicked.connect(self.ori_map_cal)
        
        self.show() 
    
    def ori_map_cal(self):
        try:
            qmin = float(self.qmin.text())
            qmax = float(self.qmax.text())
            low_thrd = float(self.low_int_thrhd_line.text())
            high_thrd = float(self.high_int_thrhd_line.text())
            t = time()
            if self.bkgd_sub_box.isChecked():
                ori_mat,wid_mat = ori_determ2d_para(
                                self.fn,self.a,self.q,qmin,qmax,
                                ll_thrhd = low_thrd, hl_thrhd = high_thrd,
                                bkgd=self.bkgd)
            else:
                ori_mat,wid_mat = ori_determ2d_para(
                                self.fn,self.a,self.q,qmin,qmax,
                                ll_thrhd = low_thrd, hl_thrhd = high_thrd)
            print('\n process took: {}'.format(time()-t))
            ori_mat[ori_mat<0] = 180 + ori_mat[ori_mat<0]
            ori_mat[ori_mat>90] = 180 - ori_mat[ori_mat>90]
            self.ori_map_widget = QWidget()
            self.ori_map_widget.setGeometry(500,700,1300,500)
            ori_map = Plot2D(self.ori_map_widget)#PlotWindow(position=True)
            ori_map.setGeometry(0,10,640,480)
            colormap={'name':'jet','normalization':'linear',
                      'autoscale':False,'vmin':0,'vamx':90}
            #colormap = setup_colormap(wid_mat)
            ori_map.addImage(ori_mat,legend='ori_map',colormap=colormap)
            #ori_map.setDefaultColormap(colormap)
            ori_map.setYAxisInverted(flag=True)
                        
            wid_map = Plot2D(self.ori_map_widget)
            wid_map.setGeometry(660,10,640,480)
            colormap = setup_colormap(1/wid_mat)
            wid_map.addImage(1/wid_mat,legend='wid_map',colormap=colormap)
            wid_map.setYAxisInverted(flag=True)
            
            self.ori_map_widget.show()
        except Exception as e:
            print(e)
            pass
    
    def load_bkgd(self):
        try:
            self.bkgd_fn,_ = QFileDialog.getOpenFileName(self,"open file","","")
            self.bkgd = fabio.open(self.bkgd_fn).data.astype(np.float)
        except:
            pass
             
    def qphi_roi_sum(self):
        try:
            qmin = float(self.qmin.text())
            qmax = float(self.qmax.text())
            amin = float(self.amin.text())
            amax = float(self.amax.text())
            #qphi = np.copy(self.qphi)
            low_thrd = float(self.low_int_thrhd_line.text())
            #qphi[qphi<low_thrd] = np.nan
            high_thrd = float(self.high_int_thrhd_line.text())
            if high_thrd != '':
                high_thrd = float(high_thrd)
            else:
                high_thrd = np.inf
                #qphi[qphi>high_thrd] = np.nan
            if self.bkgd_sub_box.isChecked():
                if not isinstance(self.bkgd,type(None)):
                    roi = qphi_roi_sum_2dmap(self.fn,q=self.q,azi=self.a,
                           qmin=qmin,qmax=qmax,amin=amin,amax=amax,
                           low_thrd=low_thrd,high_thrd=high_thrd,bkgd=self.bkgd,
                           proc_h5_list=self.proc_h5_list) 
                else:
                    roi = qphi_roi_sum_2dmap(self.fn,q=self.q,azi=self.a,
                           qmin=qmin,qmax=qmax,amin=amin,amax=amax,
                           low_thrd=low_thrd,high_thrd=high_thrd,
                           proc_h5_list=self.proc_h5_list) 
            else:
                roi = qphi_roi_sum_2dmap(self.fn,q=self.q,azi=self.a,
                           qmin=qmin,qmax=qmax,amin=amin,amax=amax,
                           low_thrd=low_thrd,high_thrd=high_thrd,
                           proc_h5_list=self.proc_h5_list) 
            self.roi = np.copy(roi)
             
            roi_map = QWidget(self)
            self.qphi_roi_map = PlotWindow(roi_map)#,position=True)
            roi_map.setGeometry(220,20,680,640)
            colormap = setup_colormap(roi)
            self.qphi_roi_map.addImage(roi,colormap=colormap,replace=True)
            self.qphi_roi_map.setYAxisInverted(flag=True)
            toolBar = qt.QToolBar()
            self.qphi_roi_map.addToolBar(
            qt.Qt.BottomToolBarArea,
            toolBar)
            position = PositionInfo(plot= \
            self.qphi_roi_map,
            converters=[('X', 
            lambda x,y: int(np.round(x))),
            ('Y',
            lambda x,y: int(np.round(y))),
            ('value',
            lambda x,y: roi[int(np.round(y)),int(np.round(x))]),
            ])
            toolBar.addWidget(position)
            
            self.subwindow = None
            self.qphi_roi_map.sigPlotSignal.connect(self.roi_map_clicked)
            self.qphi_roi_map.sigPlotSignal.connect(self.roi_map_polygon)
            roi_map.show()
            #self.subwindow.show()
         
        except Exception as e:
            print(e)
            pass
                
        
    def load_proc_hdf(self):    
        if self.fn != '':
            try:
                self.q = load_proc_dataset(self.fn,'q')
                self.a = load_proc_dataset(self.fn,'angle')
                try:
                    self.proc_h5_list = load_proc_dataset(self.fn,'proc_h5_list')
                except:
                    self.proc_h5_list = load_proc_dataset(self.fn,'proc_h5_list_1')
                self.qmin.setText(str(np.round(np.min(self.q),decimals=4)))
                self.qmax.setText(str(np.round(np.max(self.q),decimals=4)))
                self.amin.setText(str(np.round(np.min(self.a),decimals=4)))
                self.amax.setText(str(np.round(np.max(self.a),decimals=4)))
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
        
    def jl_qphi_load(self):
        self.fn,_ = QFileDialog.getOpenFileName(self,"open file","","")
        self.fileInputWin.setText(self.fn)
        self.load_proc_hdf()
    
    def roi_map_polygon(self,event):          
        try:
            #this will select all selected mask, need to select most recent
            #self.mask_roi = self.qphi_roi_map.getSelectionMask()         
            self.vertex = []
            mask = (self.roi*0+1).astype(bool)
            if event['event'] == 'drawingFinished':
                #for _ in event:
                #    print(_)
                coord = np.array(event['points'])
                coord = coord.astype(np.int)
                self.vertex.append(coord)
                #print(self.vertex)
                self.mask_roi = mask_making(mask,self.vertex)
        except Exception as e:
                print(e)
                pass
    
    def qphi_roi_ave(self):
        try:
            # have no idea why the dtype change to uint8 instance boolean here
            # thus add one line to ensure the mask is boolean type
            self.mask_roi = self.mask_roi.astype(bool)
            #print(qphi.shape,qphi.dtype,self.mask_roi.shape,self.mask_roi.dtype)
            scan_shape = load_proc_qphi_size(self.fn)
            rr,cc = np.mgrid[0:scan_shape[0],0:scan_shape[1]] 
            r = rr[self.mask_roi]
            c = cc[self.mask_roi]
            qphi = []
            for _ in range(len(r)):
                qphi.append(load_proc_single_qphi(self.fn,r[_],c[_],
                                       proc_h5_list=self.proc_h5_list))
            qphi = np.array(qphi)
            qphi[qphi==0] = np.nan
            qphi_ave = np.nanmean(qphi,axis=0).astype(float)
            if isinstance(self.ave_window,type(None)):
                self.ave_window = Plot2D(self.ave_window)
            self.ave_window.show()
            if self.bkgd_sub_box.isChecked():
                #print(qphi_ave.dtype,self.bkgd.dtype)
                qphi_ave = qphi_ave - self.bkgd.astype(np.float)
            vmin = float(self.pttn_vmin_input.text())
            vmax = float(self.pttn_vmax_input.text())
            #colormap = setup_colormap(data,vmin=vmin,vmax=vmax)
            if self.pttn_log_box.isChecked():
                colormap = setup_colormap(qphi_ave,vmin=vmin,vmax=vmax,normalization='log')
            else:
                colormap = setup_colormap(qphi_ave,vmin=vmin,vmax=vmax)

            #colormap = setup_colormap(qphi_ave)
            # need out xAxis and yAxis info this simple way doesn't work
            #if self.pttn_zoom_box.isChecked():
            #    click_reset = False
            #else:
            #    click_reset = True
            
            xlabel = r'$\rm{Q\,\,(\AA^{-1})}$'
            ylabel = r'$\rm{\phi\,\,(^{o})}$'
            self.ave_window.addImage(qphi_ave,
                colormap=colormap,
                origin=(self.q[0],self.a[0]),
                replace=True,
                #resetzoom=False,
                xlabel = xlabel,
                ylabel = ylabel,
                scale=((self.q[-1]-self.q[0])/len(self.q),
                       (self.a[-1]-self.a[0])/len(self.a))
                )
        except Exception as e:
            print(e)
            pass
     
    def roi_map_clicked(self,event):
        widget = self.qphi_roi_map.getWidgetHandle()
        if event['event'] == 'mouseClicked':
            if isinstance(self.subwindow,type(None)):
                #self.subwindow = QWidget()
                #self.subwindow.setGeometry(1260,620,900,760)
                #self.subwindow2 = PlotWindow(self.subwindow,position=True) 
                #self.subwindow2 = ImageView(self.subwindow)
                #self.subwindow.setGeometry(1260,620,680,660)
                self.subwindow = Plot2D(self.subwindow)
                self.subwindow.move(5,5)
            self.subwindow.show()
            position = widget.mapFromGlobal(qt.QCursor.pos())
            xPixel,yPixel = position.x(),position.y()
            dataPos = self.qphi_roi_map.pixelToData(xPixel,yPixel,check=True)
            try:
                col,row = (int(dataPos[0]),int(dataPos[1]))
                print('\n pttn position: x-{} y-{}'.format(col,row))
                
                if self.bkgd_sub_box.isChecked():
                    if isinstance(self.bkgd,type(None)):
                        data = np.copy(load_proc_single_qphi(self.fn,row,col,
                                    proc_h5_list=self.proc_h5_list,
                                    proc_type="integrate2d"))
                    else:
                        bkgd = np.copy(self.bkgd).astype(np.float)
                        data = np.copy(load_proc_single_qphi(self.fn,row,col,
                                    proc_h5_list=self.proc_h5_list,
                                    proc_type="integrate2d"))
                        data = data - bkgd
                else:
                    data = np.copy(load_proc_single_qphi(self.fn,row,col,
                                proc_h5_list=self.proc_h5_list,
                                proc_type="integrate2d"))
                data[np.isnan(data)] = 0
                data[np.isinf(data)] = 0
                vmin = float(self.pttn_vmin_input.text())
                vmax = float(self.pttn_vmax_input.text())
                if self.pttn_log_box.isChecked():
                    colormap = setup_colormap(data,vmin=vmin,vmax=vmax,normalization='log')
                else:
                    colormap = setup_colormap(data,vmin=vmin,vmax=vmax)
    
                #colormap = setup_colormap(data)
                if self.pttn_zoom_box.isChecked():
                    click_reset = False
                else:
                    click_reset = True
                xlabel = r'$\rm{Q\,\,(\AA^{-1})}$'
                ylabel = r'$\rm{\phi\,\,(^{o})}$'
                self.subwindow.addImage(data,
                origin=(self.q[0],self.a[0]),
                xlabel = xlabel,
                ylabel = ylabel,
                scale=((self.q[-1]-self.q[0])/len(self.q),
                       (self.a[-1]-self.a[0])/len(self.a)),
                colormap = colormap,
                resetzoom=click_reset,
                )
                #self.subwindow2.setColormap(colormap)
                self.subwindow.setYAxisInverted(flag=True)
                
                #print(self.subwindow.getXAxis.getLimits())
                #toolBar = qt.QToolBar()
                #self.subwindow.addToolBar(
                #qt.Qt.BottomToolBarArea,
                #toolBar)
                #position = PositionInfo(plot= \
                #self.subwindow,
                #converters=[('X', lambda x,y: x),
                #('Y', lambda x,y: y)])
                #toolBar.addWidget(position)
            except Exception as e:
                print(e)
        else:
            pass
        
        
