import sys
sys.path.append('../xs_proc')
from h5_data_search import *
from xs_data_proc import *
from proc_data_ana import *     
from visual_func import * 
from multi_scan_qphi_proc import * 
from multi_scan_Iq_proc import * 
import os
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

def setup_colormap(data):
    data = data.astype(np.float)
    vmin = np.nanmedian(data) - np.nanstd(data)*3
    if vmin < np.nanmin(data):
        vmin = np.nanmin(data)
    vmax = np.nanmedian(data) + np.nanstd(data)*3
    if vmax > np.nanmax(data)
        vmax = np.nanmax(data)
    colormap = Colormap(name='gray',
                normalization='linear',
                vmin = vmin,
                vmax = vmax
                )
    return colormap
                

class input_window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.data_path = 'entry_0000/measurement/data'
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
        
    def load_id13_h5(self):
        try:
            
            self.scan_obj = id13_h5_search(
                                self.proposal_path,
                                self.proposal_id)
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
                        #self.h5_list,
                        #self.path_idx,self.pttn_idx,
                        #self.scan_shape,self.scan_obj)
            self.qphi_cal_win.show()
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass
    
    def Iq_cal_proc(self):
        try:
            self.Iq_cal_win = Iq_calculate(self.scan_obj)
                        #self.h5_list,
                        #self.path_idx,self.pttn_idx,
                        #self.scan_shape,self.scan_obj)
            self.Iq_cal_win.show()
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass
    
    def qphi_ana_proc(self):
        try:
            self.qphi_ana_win = qphi_analysis()
                        #self.h5_list,
                        #self.path_idx,self.pttn_idx,
                        #self.scan_shape)

            self.qphi_ana_win.show()
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass 

            
class roi_sum(QWidget):
    def __init__(self,obj):#h5_list,path_idx,pttn_idx,scan_shape,*args,**kwargs):
        super().__init__()
        #self.h5_list  = h5_list
        #self.path_idx = path_idx
        #self.pttn_idx = pttn_idx
        #self.scan_shape = scan_shape
        self.obj = obj
        self.roi_sum_UI()
        self.show()
        
    def roi_sum_UI(self):
        self.setGeometry(500,40,900,650)
        self.setWindowTitle('roi_sum')
        
        up_corner = QWidget(self)
        up_corner.setGeometry(20,20,180,20)
        up_corner_label = QLabel('up corner',up_corner)
        up_corner_label.setGeometry(0,0,60,20)
        up_corner_input = QLineEdit('0,0',up_corner)
        self.ul_corner_set(up_corner_input.text())
        up_corner_input.setGeometry(80,0,80,20)
        up_corner_input.textChanged[str].connect(self.ul_corner_set)
        
        dr_corner = QWidget(self)
        dr_corner.setGeometry(20,60,180,20)
        dr_corner_label = QLabel('down corner',dr_corner)
        dr_corner_label.setGeometry(0,0,60,20)
        dr_corner_input = QLineEdit('-1,-1',dr_corner)
        self.dr_corner_set(dr_corner_input.text())
        dr_corner_input.setGeometry(80,0,80,20)
        dr_corner_input.textChanged[str].connect(self.dr_corner_set)
        
        proc_button = QPushButton('process',self)
        proc_button.setGeometry(20,120,80,20)
        proc_button.clicked.connect(self.roi_sum_process)
        
        pttn_roi_sum_button = QPushButton('roi sum pttn',self)
        pttn_roi_sum_button.setGeometry(20,150,80,20)
        pttn_roi_sum_button.clicked.connect(self.pttn_roi_sum_process)
        self.pttn_roi_sum_window = None
         
        self.bkgd_load_bttn = QPushButton('bkgd load',self)
        self.bkgd_load_bttn.move(10,220)
        self.bkgd_load_bttn.clicked.connect(self.load_bkgd)
        
        self.bkgd_sub_box = QCheckBox('bkgd subtract',self)
        self.bkgd_sub_box.move(10,250)
        
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
            res = parallel_func(scan_pttn_roi_sum,12,
                                np.arange(len(self.path_idx.flatten())),
                                h5_list  = self.h5_list,
                                path_idx = self.path_idx.flatten(),
                                pttn_idx = self.pttn_idx.flatten(),
                                data_path = 'entry_0000/measurement/data',
                                left_top  = self.ul,
                                right_bottom = self.dr,
                                )
            self.roi_map = np.array(res).reshape(self.scan_shape)    
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
            #self.subwindow1 = PlotWindow(position=True)
            #self.subwindow1.show()
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
                    colormap = setup_colormap(data)
                    self.subwindow1.addImage(data,colormap=colormap)
                    self.subwindow1.setYAxisInverted(flag=True)
                      
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
            if isinstance(self.pttn_roi_sum_window,type(None)): 
                self.pttn_roi_sum_window = Plot2D()
            colormap = setup_colormap(sum_pttn)
            self.pttn_roi_sum_window.addImage(sum_pttn,colormap=colormap)
            self.pttn_roi_sum_window.show()
            print(time()-t)
        except Exception as e:
            print(e)
            pass 
    
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
        data_path_label = QLabel("data_path:",Iq_shape)
        data_path_label.setGeometry(10,70,60,20)
        self.data_path_input = QLineEdit("entry_0000/measurement/data",Iq_shape)
        self.data_path_input.setGeometry(80,70,120,20)
        save_path_label = QLabel("save_path:",Iq_shape)
        save_path_label.setGeometry(10,100,60,20)
        self.save_path_input = QLineEdit("/data/id13/inhouse12/jiliang/code_v3/example/sc5005/",Iq_shape)
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
        
    def load_Iq(self):
        try:
            self.Iq_h5,_ = QFileDialog.getOpenFileName(self,"load Iq","","")
            self.Iq = load_proc_dataset(self.Iq_h5,'Iq',proc_type="integrate1d")
            self.q  = load_proc_dataset(self.Iq_h5,'q',proc_type="integrate1d")
        except Exception as e:
            print(e)
            pass
        
    def Iq_map_show(self):
        Iq_show = Plot2D(self)#PlotWindow(self)#,position=True)
        Iq_show.setGeometry(240,20,640,480)
        data = np.copy(self.Iq)
        data = data.reshape((self.Iq.shape[0]*self.Iq.shape[1],
                             self.Iq.shape[2]))

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
            self.data_path = str(self.data_path_input.text())
            num_core = int(self.core_num.text())
            
            res = parallel_func(scan_calculate_Iq,
                           num_core,
                           np.arange(len(self.path_idx.flatten())),
                           h5_list   = self.h5_list,
                           path_idx  = self.path_idx.flatten(),
                           pttn_idx  = self.pttn_idx.flatten(),
                           data_path = self.data_path,
                           pyfai_obj = self.ai,
                           mask      = self.mask,
                           q_npts    = self.q_npts,
                           **kwargs
                           )
            q    = res[0][0]
            q   /= 10
            scan_shape = self.scan_shape
            Iq   = np.zeros((scan_shape[0],
                             scan_shape[1],
                             len(res[0][1]),
                            ))
            for _ in range(len(res)):
                #if len(idx_list) > 1:
                #    i1 = int((_+i*t1.single_h5_shape[0])/scan_shape[1])
                #    i2 = int((_+i*t1.single_h5_shape[0])%scan_shape[1])
                #else:
                i1 = int(_/scan_shape[1])
                i2 = int(_%scan_shape[1])
                Iq[i1,i2,:]   = res[_][1]
            name = self.obj._data_name[0].split('.')[0]
            save_path =self.save_path_input.text()
            print(time()-t)
            save_Iq_as_h5(self.obj,save_path,name,
                    q    = q,
                    Iq   = Iq,
                    path_idx = self.path_idx,
                    pttn_idx = self.pttn_idx,
                    h5_path_list = self.h5_list)
            print(time()-t)#,'\n\nfuck complete')
            self.q  = q
            self.Iq = Iq 
        except Exception as e:
            print(e)
            pass


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
        mask_load.setGeometry(20,120,240,80)
        mask_file_load = QPushButton('load mask',mask_load)
        mask_file_load.clicked.connect(self.load_mask_window)
        #if isinstance(self.mask_file,type(None)):
        #    print("mask is needed, at current develop version")
        #    return
        #else:
        #    self.mask = np.load(self.mask_file)['mask']
        
        qphi_shape = QWidget(self)
        qphi_shape.setGeometry(20,240,240,120)
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
        data_path_label = QLabel("data_path:",qphi_shape)
        data_path_label.setGeometry(10,100,60,20)
        self.data_path_input = QLineEdit("entry_0000/measurement/data",qphi_shape)
        self.data_path_input.setGeometry(80,100,120,20)
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
            self.data_path = str(self.data_path_input.text())
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
        #self.bkgd_sub_box.stateChanged.connect(self.qphi_roi_sum)
        
        self.ori_map_bttn = QPushButton('ori map',roi_setup)
        self.ori_map_bttn.move(10,280)
        self.ori_map_bttn.clicked.connect(self.ori_map_cal)
        
        self.show() 
    
    def ori_map_cal(self):
        try:
            qmin = float(self.qmin.text())
            qmax = float(self.qmax.text())
            qphi = np.copy(self.qphi).astype(np.float)
            if self.bkgd_sub_box.isChecked():
                if not isinstance(self.bkgd,type(None)):
                    bkgd = np.copy(self.bkgd).reshape(1,1,qphi.shape[-2],qphi.shape[-1])
                    qphi = qphi - np.tile(bkgd,(qphi.shape[0],qphi.shape[1],1,1))
            ori_mat,wid_mat = ori_determ2d(qphi,self.a,self.q,qmin,qmax,
                                ll_thrhd = 0, hl_thrhd = np.inf)
            #
            self.ori_map_widget = QWidget()
            self.ori_map_widget.setGeometry(500,700,1300,500)
            ori_map = Plot2D(self.ori_map_widget)#PlotWindow(position=True)
            ori_map.setGeometry(0,10,640,480)
            colormap={'name':'gray','normalization':'linear',
                                        'autoscale':False,'vmin':0,'vamx':180}
            ori_map.addImage(ori_mat,legend='ori_map',colormap=colormap)
            #ori_map.setDefaultColormap(colormap)
            
            wid_map = Plot2D(self.ori_map_widget)
            wid_map.setGeometry(660,10,640,480)
            colormap = setup_colormap(wid_mat)
            wid_map.addImage(wid_mat,legend='wid_map',colormap=colormap)
            self.ori_map_widget.show()
        except Exception as e:
            print(e)
            pass
    
    def load_bkgd(self):
        try:
            self.bkgd_fn,_ = QFileDialog.getOpenFileName(self,"open file","","")
            self.bkgd = fabio.open(self.bkgd_fn).data
        except:
            pass
             
    def qphi_roi_sum(self):
        try:
            qmin = float(self.qmin.text())
            qmax = float(self.qmax.text())
            amin = float(self.amin.text())
            amax = float(self.amax.text())
            qphi = np.copy(self.qphi)
            low_thrd = float(self.low_int_thrhd_line.text())
            qphi[qphi<low_thrd] = np.nan
            high_thrd = self.high_int_thrhd_line.text()
            if high_thrd != '':
                high_thrd = float(high_thrd)
                qphi[qphi>high_thrd] = np.nan
            roi = sum_roi_2dmap(qphi,q=self.q,a=self.a,
                                qmin=qmin,qmax=qmax,
                                amin=amin,amax=amax)
            self.roi = np.copy(roi)
            #self.polygon_choose = QtCheckBox("polygon vertx",self)
            #self.polygon_choose.setGeometry(80,160,100,20)
            #self.polygon_choose.isChecked.connect(self.add_vertex)
            
            roi_map = QWidget(self)
            self.qphi_roi_map = PlotWindow(roi_map)#,position=True)
            roi_map.setGeometry(220,20,680,640)
            colormap = setup_colormap(roi)
            self.qphi_roi_map.addImage(roi,colormap=colormap,replace=True)
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
                self.qphi = load_proc_dataset(self.fn,'map_qphi')
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
            self.mask_roi = self.qphi_roi_map.getSelectionMask()         
            self.vertex = []
            mask = (self.roi*0+1).astype(bool)
            if event['event'] == 'drawingFinished':
                coord = np.array(event['points'])
                coord = coord.astype(np.int)
                self.vertex.append(coord)
                #print(vertex)
                self.mask_roi = mask_making(mask,self.vertex)
                # did't find the mask signal event, drawing signal emit before mask update
                #self.mask_roi = self.qphi_roi_map.getSelectionMask()
                #self.mask_roi = self.qphi_roi_map._maskToolsDockWidget.getSelectionMask()
                #self.mask_roi = self.mask_roi.astype(bool)
                #print(self.mask_roi.dtype)
                # for test correctness of mask
                #self.mask_display = Plot2D()
                #self.mask_display.addImage(self.mask_roi)
                #self.mask_display.show()
        except Exception as e:
                print(e)
                pass
    
    def qphi_roi_ave(self):
        try:
            qphi = np.copy(self.qphi)
            qphi[qphi==0] = np.nan
            # have no idea why the dtype change to uint8 instance boolean here
            # thus add one line to ensure the mask is boolean type
            self.mask_roi = self.mask_roi.astype(bool)
            #print(qphi.shape,qphi.dtype,self.mask_roi.shape,self.mask_roi.dtype)
            qphi_ave = np.nanmean(qphi[self.mask_roi],axis=0).astype(float)
            if isinstance(self.ave_window,type(None)):
                self.ave_window = Plot2D(self.ave_window)
            self.ave_window.show()
            if self.bkgd_sub_box.isChecked():
                #print(qphi_ave.dtype,self.bkgd.dtype)
                qphi_ave = qphi_ave - self.bkgd.astype(np.float)
            colormap = setup_colormap(qphi_ave)
            self.ave_window.addImage(qphi_ave,
                colormap=colormap,
                origin=(0,-180),
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
                self.subwindow = QWidget()
                self.subwindow.setGeometry(1260,620,900,760)
                #self.subwindow2 = PlotWindow(self.subwindow,position=True) 
                self.subwindow2 = ImageView(self.subwindow)
                self.subwindow2.move(5,5)
            self.subwindow.show()
            position = widget.mapFromGlobal(qt.QCursor.pos())
            xPixel,yPixel = position.x(),position.y()
            dataPos = self.qphi_roi_map.pixelToData(xPixel,yPixel,check=True)
            try:
                col,row = (int(dataPos[0]),int(dataPos[1]))
                if self.bkgd_sub_box.isChecked():
                    if isinstance(self.bkgd,type(None)):
                        data = np.copy(self.qphi[row,col])
                    else:
                        bkgd = np.copy(self.bkgd).astype(np.float)
                        data = np.copy(self.qphi[row,col])\
                                 - bkgd
                else:
                    data = np.copy(self.qphi[row,col])
                #self.subwindow2.addImage(data,
                colormap = setup_colormap(data)
                self.subwindow2.setImage(data,
                origin=(0,-180),
                scale=((self.q[-1]-self.q[0])/len(self.q),
                       (self.a[-1]-self.a[0])/len(self.a))
                )
                self.subwindow2.setColormap(colormap)
                self.subwindow2.setYAxisInverted(flag=True)
                
                toolBar = qt.QToolBar()
                self.subwindow2.addToolBar(
                qt.Qt.BottomToolBarArea,
                toolBar)
                position = PositionInfo(plot= \
                self.subwindow2,
                converters=[('X', lambda x,y: x),
                ('Y', lambda x,y: y)])
                toolBar.addWidget(position)
                #self.subwindow2.setLimits(self.q[0],
                #            self.q[-1],self.a[0],self.a[-1])
                
                # this is for raw x,y input without rescale
                #toolBar = qt.QToolBar()
                #
                #self.subwindow2.addToolBar(
                #qt.Qt.BottomToolBarArea,
                #toolBar)
                #
                #position = PositionInfo(plot= \
                #self.subwindow2,
                #converters=[('Q', 
                #lambda x,y: self.q[int(np.round(x))]),
                #('Azimuth',
                #lambda x,y: self.a[int(np.round(y))])])
                #
                #toolBar.addWidget(position)
            except Exception as e:
                print(e)
        else:
            pass
#class PlotCanvas(FigureCanvas):
#    def __init(self,data):
#        super().__init__()
         
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = input_window()
    sys.exit(app.exec_())
        
        
