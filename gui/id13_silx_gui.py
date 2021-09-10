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

from roi_sum_UI import roi_sum
from Iq_calculate_UI import Iq_calculate
from qphi_calculate_UI import qphi_calculate
from qphi_analysis_UI import qphi_analysis
from XRF_load_UI import XRF_UI
from colormap_setup import setup_colormap
                

class input_window(QMainWindow):
    def __init__(self):
        super().__init__()
        #self.data_path = 'entry_0000/measurement/data'
        self.data_path = 'entry/data/data'
        self.input_UI()
    
    def input_UI(self):
        self.setGeometry(100,100,1280,480)
        self.setWindowTitle('jl_beta_ver')#fuck')
        
        path_load = QWidget(self)
        path_load.setGeometry(20,20,180,30)
        path_load_bttn = QPushButton('path load',path_load)
        path_load_bttn.clicked.connect(self.load_path)
            
        path = QWidget(self)
        path.setGeometry(20,60,180,20)
        path_name = QLabel('path',path)
        #self.path_input = QLineEdit('',path)
        #self.path_input = QLineEdit('/data/id13/inhouse12/DATAPOLICY_I12_1/eh2/inhouse/blc12867/id13',path)
        self.path_input = QLineEdit('/data/visitor/sc5005/id13',path)
        self.proposal_path = self.path_input.text()
        self.path_input.move(60,0)
        self.path_input.textEdited.connect(self.path_set)

        prop_id = QWidget(self) 
        prop_id.setGeometry(20,90,180,20)
        id_name = QLabel('proposal',prop_id)
        #self.id_input   = QLineEdit('',prop_id)
        #self.id_input   = QLineEdit('blc12867',prop_id)
        self.id_input   = QLineEdit('sc5005',prop_id)
        self.proposal_id = self.id_input.text()
        self.id_input.move(60,0)
        self.id_input.textEdited.connect(self.id_set)

        h5_data_path = QWidget(self) 
        h5_data_path.setGeometry(20,120,180,20)
        h5_data_name = QLabel('data path',h5_data_path)
        self.h5_data_path_input   = QLineEdit(self.data_path,h5_data_path)
        self.data_h5path = self.h5_data_path_input.text()
        self.h5_data_path_input.move(60,0)
        self.h5_data_path_input.textEdited.connect(self.h5_data_path_set)

       
        self.tree = QTreeWidget(self)
        self.tree.setGeometry(220,20,360,340)
        header = QTreeWidgetItem(['File'])
        self.tree.setHeaderItem(header)
        
        self.msgbox = QListWidget(self)
        self.msgbox.insertItem(0,'scan info:\n')
        self.msgbox.setGeometry(600,20,360,160)
        
        self.proc_msg = QTextBrowser(self)
        self.proc_msg.setGeometry(600,200,360,160)
        
        info_enter = QPushButton("enter",self)
        info_enter.setGeometry(20,380,60,20)
        info_enter.clicked.connect(self.load_id13_h5)
        
        proc_widget = QWidget(self)
        proc_widget.setGeometry(1000,20,260,320)
        roi_sum_proc = QPushButton("roi sum",proc_widget)
        roi_sum_proc.setGeometry(10,10,60,20)
        roi_sum_proc.clicked.connect(self.roi_sum_proc)
        
        qphi_cal = QPushButton("qphi cal",proc_widget)
        qphi_cal.setGeometry(80,10,60,20)
        qphi_cal.clicked.connect(self.qphi_cal_proc)
        
        qphi_ana_proc = QPushButton("qphi ana",proc_widget)
        qphi_ana_proc.setGeometry(10,40,60,20)
        qphi_ana_proc.clicked.connect(self.qphi_ana_proc)
       
        pttn_sum_proc = QPushButton("total sum",proc_widget)
        pttn_sum_proc.setGeometry(80,40,60,20)
        pttn_sum_proc.clicked.connect(self.total_sum_proc)
         
        Iq_cal = QPushButton("Iq cal",proc_widget)
        Iq_cal.setGeometry(10,70,60,20)
        Iq_cal.clicked.connect(self.Iq_cal_proc)
        
        XRF_load = QPushButton("load XRF",proc_widget)
        XRF_load.setGeometry(10,100,60,20)
        XRF_load.clicked.connect(self.XRF_load)
        
        self.show()
    
    def total_sum_proc(self):
        try:
            if hasattr(self.scan_obj,'data_h5'):
                t = time()
                h5_path_list = self.scan_obj.data_h5
                for i in range(len(h5_path_list)):
                    h5_path = h5_path_list[i]
                    name = os.path.split(h5_path)[1][:-3]
                    with h5py.File(h5_path,'r') as f:
                        total_num = f[self.data_path].shape[0]
                    idx_list = chunk_idx(total_num,2000)
                    if total_num >= 2000:
                        res =  parallel_func(partial_sum,12,idx_list,
                                             h5_path=h5_path,data_path=self.data_path,
                                             )
                        print("\nsum processing time:\n%5.2f sec" %(time()-t))
                        data = res[0]*0
                        for __ in res:
                            data += __
                        data /= len(res)
                    else:
                        print("\nsum processing time:\n%5.2f sec" %(time()-t))
                        with h5py.File(h5_path,"r") as f:
                            data = np.nanmean(f[self.data_path][:],axis=0)
                    data = data.astype(np.uint32)
                    img = fabio.edfimage.EdfImage(data=data)
                    #
                    save_path = QFileDialog.getExistingDirectory()
                    os.chdir(save_path)
                    if not os.path.exists("edf_sum_file"):
                        os.mkdir("edf_sum_file")
                    os.chdir("edf_sum_file")
                    img.write("{}.edf".format(name))
                    os.chdir(save_path)
        except Exception as e:
            print(e)
            pass    


    def load_path(self):
        try:
            path = QFileDialog.getExistingDirectory()
            self.path_input.setText(path)
            proposal_id = path.split('/')[-2]
            self.id_input.setText(proposal_id)
            self.path_set()
            self.id_set()
        except Exception as e:
            print(e)
            pass    
    
    def path_set(self):
        self.proposal_path = self.path_input.text()

    def id_set(self):
        self.proposal_id   = self.id_input.text()
        
    def h5_data_path_set(self):
        self.data_h5path   = self.h5_data_path_input.text()
    
    def load_id13_h5(self):
        try:
            
            self.scan_obj = id13_h5_search(
                                self.proposal_path,
                                self.proposal_id,
                                data_h5path = self.data_h5path)
            self.add_file()
        except Exception as e:
            print(e)
            pass
    
    def add_file(self):
        try:
            path = os.path.join(self.proposal_path,
                   '{}_id13.h5'.format(self.proposal_id))
            with h5py.File(path,'r') as f:
                self.f = f
                self.filename = self.f.filename.split('/')[-1]
                self.tree_root = QTreeWidgetItem(self.tree,
                                    [self.filename])
                self.add_branch(self.tree_root,self.f)
                self.tree.itemClicked.connect(self.onItemClicked)
        except Exception as e:
            print(e)
            pass
     
    def add_branch(self,tree_root,h5file):
        for _ in h5file:
            try:
                branch = QTreeWidgetItem([str(_)])
                tree_root.addChild(branch)
            except Exception as e:
                print(path,'\n\n',e,'\n\n')
                pass
    
    
    @pyqtSlot(QTreeWidgetItem,int)
    def onItemClicked(self,item):
        self.sample_keyword = item.text(0)
        try:
            self.scan_obj.keyword_search(self.sample_keyword)
            self.msgbox.clear()
            self.proc_msg.clear()
            self.msgbox.insertItem(0,'scan info:')
            for _,__ in enumerate(self.scan_obj.data_info.keys()):
                msg = str(__) + ':  ' + \
                      str(self.scan_obj.data_info[__])
                self.msgbox.insertItem((_+1),msg)
            self.total_pttns,self.scan_shape,self.idx_list = \
                    scan_info(self.scan_obj)
            self.h5_list,self.path_idx,self.pttn_idx = \
                    scan_h5_data_info(self.scan_obj,
                               self.scan_shape,self.idx_list)
            #self.update_info(qphi_calculate)
            #self.update_info(Iq_calculate)           
            #print(self.h5_list,self.scan_obj._data_name,self.scan_obj.data_h5) 
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass

    def update_info(self,obj):
        obj.h5_list    = self.h5_list
        obj.path_idx   = self.path_idx
        obj.pttn_idx   = self.pttn_idx
        obj.scan_shape = self.scan_shape
        obj.obj        = self.scan_obj 
    
    def roi_sum_proc(self):
        try:
            self.roi_sum_win = roi_sum(self.scan_obj)
                        #self.h5_list,
                        #self.path_idx,self.pttn_idx,
                        #self.scan_shape)
            self.roi_sum_win.show()
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass
    
    def qphi_cal_proc(self):
        try:
            self.qphi_cal_win = qphi_calculate(self.scan_obj)
            self.qphi_cal_win.show()
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass
    
    def Iq_cal_proc(self):
        try:
            self.Iq_cal_win = Iq_calculate(self.scan_obj)
            self.Iq_cal_win.show()
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass
    
    def XRF_load(self):
        try:
            self.XRF_load_win = XRF_UI(self.scan_obj)
            self.XRF_load_win.show()
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass
    
    def qphi_ana_proc(self):
        try:
            self.qphi_ana_win = qphi_analysis()

            self.qphi_ana_win.show()
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass 

         
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = input_window()
    sys.exit(app.exec_())
        
        
