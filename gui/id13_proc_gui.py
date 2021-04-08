import sys
sys.path.append('../xs_proc')
from h5_data_search import *
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

class input_window(QWidget):
    def __init__(self):
        super().__init__()
        self.input_UI()
    
    def input_UI(self):
        self.setGeometry(100,100,1280,480)
        self.setWindowTitle('fuck')
        
        path = QWidget(self)
        path.setGeometry(20,20,180,20)
        path_name = QLabel('path',path)
        self.path_input = QLineEdit('',path)
        self.path_input.move(60,0)
        self.path_input.textEdited.connect(self.path_set)

        prop_id = QWidget(self) 
        prop_id.setGeometry(20,50,180,20)
        id_name = QLabel('proposal',prop_id)
        self.id_input   = QLineEdit('',prop_id)
        self.id_input.move(60,0)
        self.id_input.textEdited.connect(self.id_set)
       
        self.tree = QTreeWidget(self)
        self.tree.setGeometry(220,20,360,320)
        header = QTreeWidgetItem(['File'])
        self.tree.setHeaderItem(header)
        
        info_enter = QPushButton("enter",self)
        info_enter.setGeometry(20,380,60,20)
        info_enter.clicked.connect(self.load_id13_h5)
        
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
            print(1)
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
        self.scan_obj.keyword_search(self.sample_keyword)
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = input_window()
    sys.exit(app.exec_())
        
        
