import sys
sys.path.append('../xs_proc')
from h5_data_search import *
from xs_data_proc import *
from proc_data_ana import *       
import os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import h5py
import numpy as np
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import file_tree as ft
from time import time

from silx.gui.plot.PlotWidget import PlotWidget
from silx.gui import qt
from silx.gui.plot.PlotWindow import PlotWindow 
from silx.gui.plot.PlotTools import PositionInfo

class input_window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.data_path = 'entry_0000/measurement/data'
        self.input_UI()
    
    def input_UI(self):
        self.setGeometry(100,100,1280,480)
        self.setWindowTitle('fuck')
        
        path = QWidget(self)
        path.setGeometry(20,20,180,20)
        path_name = QLabel('path',path)
        #self.path_input = QLineEdit('',path)
        self.path_input = QLineEdit('/data/id13/inhouse12/DATAPOLICY_I12_1/eh2/inhouse/blc12867/id13',path)
        self.proposal_path = self.path_input.text()
        self.path_input.move(60,0)
        self.path_input.textEdited.connect(self.path_set)

        prop_id = QWidget(self) 
        prop_id.setGeometry(20,50,180,20)
        id_name = QLabel('proposal',prop_id)
        #self.id_input   = QLineEdit('',prop_id)
        self.id_input   = QLineEdit('blc12867',prop_id)
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
        self.show()
    
    def path_set(self):
        self.proposal_path = self.path_input.text()

    def id_set(self):
        self.proposal_id   = self.id_input.text()
        
    def load_id13_h5(self):
        try:
            #print('\n\n',self.proposal_path,
            #      '\n\n',self.proposal_id)
            
            self.scan_obj = id13_h5_search(
                                self.proposal_path,
                                self.proposal_id)
            self.add_file()
        except Exception as e:
            print(e)
            pass
    
    def add_file(self):
        try:
            #print(1)
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
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass

    def roi_sum_proc(self):
        try:
            self.roi_sum_win = roi_sum(self.h5_list,
                        self.path_idx,self.pttn_idx,
                        self.scan_shape)
            self.roi_sum_win.show()
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass

            
class roi_sum(QWidget):
    def __init__(self,h5_list,path_idx,pttn_idx,scan_shape):
        super().__init__()
        self.h5_list  = h5_list
        self.path_idx = path_idx
        self.pttn_idx = pttn_idx
        self.scan_shape = scan_shape
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
        proc_button.setGeometry(20,120,60,20)
        proc_button.clicked.connect(self.roi_sum_process)
        
        
        
        self.show()
         
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
            self.roi_show = PlotWindow(self,position=True)
            self.roi_show.addImage(self.roi_map)
            #print(self.roi_show.sigPlotSignal.emit('mouseMoved'))
            #from silx.gui.plot.ImageView import ImageView
            #self.roi_show = ImageView(self)
            #self.roi_show.setImage(self.roi_map)
            self.position = PositionInfo(plot=self.roi_show)
            self.position._plotEvent({'event':'mouseClicked'})
            self.roi_show.sigPlotSignal.emit({'event':'mouseClicked'})
            self.roi_show.sigPlotSignal.connect(self.roi_map_clicked)
            self.roi_show.move(200,20)
            self.subwindow1 = PlotWindow(position=True)
            self.subwindow1.show()
            self.roi_show.show()
            print(time() - t)
        except Exception as e:
            print(str(e))
            pass
    
    def roi_map_clicked(self):
        widget = self.roi_show.getWidgetHandle()
        position = widget.mapFromGlobal(qt.QCursor.pos())
        xPixel,yPixel = position.x(),position.y()
        dataPos = self.roi_show.pixelToData(xPixel,yPixel,check=True)
        try:
            col,row = (int(dataPos[0]),int(dataPos[1]))
            with h5py.File(self.h5_list[self.path_idx[row,col]],'r') as f:
                data = f['entry_0000/measurement/data'][self.pttn_idx[row,col]]
                self.subwindow1.addImage(data)
        except Exception as e:
            print(e)
        #print(self.roi_show.pixelToData(12,30))
    #def pttn_plot(self):
    #    c,r = cl
                            
#class PlotCanvas(FigureCanvas):
#    def __init(self,data):
#        super().__init__()
         
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = input_window()
    sys.exit(app.exec_())
        
        
