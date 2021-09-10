import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import sys
import h5py
import numpy as np
import fabio

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from time import time

from load_id13_position_info import scan_pos_info

class input_window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.input_UI()

    def input_UI(self):
        self.setGeometry(100,100,980,480)
        self.setWindowTitle('display frame pos info of scan')

        file_load = QWidget(self)
        file_load.setGeometry(20,20,180,30)
        file_load_bttn = QPushButton('load h5',file_load)
        file_load_bttn.clicked.connect(self.load_file)

        frame_num_win= QWidget(self)
        frame_num_win.setGeometry(20,60,180,20)
        frame_num_label = QLabel('Frame No.:',frame_num_win)
        self.frame_num_input = QLineEdit('0',frame_num_win)
        self.frame_num = int(self.frame_num_input.text())
        self.frame_num_input.move(60,0)
        self.frame_num_input.textEdited.connect(self.frame_num_set)

        self.tree = QTreeWidget(self)
        self.tree.setGeometry(220,20,360,340)
        header = QTreeWidgetItem(['File'])
        self.tree.setHeaderItem(header)
        
        self.msgbox = QListWidget(self)
        self.msgbox.insertItem(0,'scan info:\n')
        self.msgbox.setGeometry(600,20,360,160)

        self.proc_msg = QTextBrowser(self)
        self.proc_msg.setGeometry(600,200,360,160)
        
        info_update = QPushButton("update",self)
        info_update.setGeometry(20,340,60,20)
        info_update.clicked.connect(self.update_id13_h5)
        
        info_enter = QPushButton("enter",self)
        info_enter.setGeometry(20,380,60,20)
        info_enter.clicked.connect(self.load_id13_h5)
        self.show()

    def load_file(self):
        try:
            self.h5_fn,_ = QFileDialog.getOpenFileName(self,"load_file","","*.h5")
        except Exception as e:
            print(e)
            pass

    def frame_num_set(self):
        try:
            if self.frame_num_input.text() == '':
                self.frame_num = 0
            self.frame_num = int(self.frame_num_input.text())
            self.scan_obj = scan_pos_info(
                                self.h5_fn,
                                self.frame_num)
            self.scan_obj.keyword_search(self.sample_keyword)
            self.msgbox.clear()
            self.proc_msg.clear()
            self.msgbox.insertItem(0,'scan info:')
            for _,__ in enumerate(self.scan_obj.data_info.keys()):
                msg = str(__) + ':  ' + \
                      str(self.scan_obj.data_info[__])
                self.msgbox.insertItem((_+1),msg)
            
            pos_msg = ('\n####################\n'+
                       f'frame number: {self.frame_num}\n'+
                       "{}, {}, {:.4f}, {:.4f}".format(self.scan_obj.mot1, 
                                                        self.scan_obj.mot2, 
                                                        self.scan_obj.pos1[self.frame_num],
                                                        self.scan_obj.pos2[self.frame_num])
                      +'\n####################\n')
            self.proc_msg.setText(pos_msg)
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass
        
    def load_id13_h5(self):
        try:
            self.add_file()
        except Exception as e:
            print(e)
            pass

    def update_id13_h5(self):
        try:
            self.tree.clear()
            self.add_file()
        except Exception as e:
            print(e)
            pass
    
    def add_file(self):
        try:
            with h5py.File(self.h5_fn,'r') as f:
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
            self.scan_obj = scan_pos_info(
                                self.h5_fn,
                                self.frame_num)
            self.scan_obj.keyword_search(self.sample_keyword)
            self.msgbox.clear()
            self.proc_msg.clear()
            self.msgbox.insertItem(0,'scan info:')
            for _,__ in enumerate(self.scan_obj.data_info.keys()):
                msg = str(__) + ':  ' + \
                      str(self.scan_obj.data_info[__])
                self.msgbox.insertItem((_+1),msg)
            
            pos_msg = ('\n####################\n'+
                       f'frame number: {self.frame_num}\n'+
                       "{}, {}, {:.4f}, {:.4f}".format(self.scan_obj.mot1, 
                                                        self.scan_obj.mot2, 
                                                        self.scan_obj.pos1[self.frame_num],
                                                        self.scan_obj.pos2[self.frame_num])
                      +'\n####################\n')
            self.proc_msg.setText(pos_msg)
        except Exception as e:
            self.proc_msg.setText(str(e))
            pass

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = input_window()
    sys.exit(app.exec_())

                                       
