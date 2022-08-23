import sys
from PyQt5.QtWidgets import (QPushButton, QHBoxLayout, QApplication, QWidget, QTreeView, QTreeWidget, QVBoxLayout, QTreeWidgetItem,QFileDialog)
from PyQt5.QtGui import QStandardItemModel,QStandardItem
from PyQt5.QtCore import Qt, pyqtSlot,QSize
import h5py
import numpy as np

class tree(QWidget):
    def __init__(self):
        # inital inherited QWidget
        super().__init__()
        self.title = 'Tree of h5 data'
        self.left = 200
        self.top = 20
        self.width = 320
        self.height = 280
        
        self.setWindowTitle(self.title)
        self.setGeometry(self.left,self.top,self.width,self.height)
        
        #self.datalayout= QVBoxLayout()
        self.tree = QTreeWidget()
        header = QTreeWidgetItem(['File'])
        self.tree.setHeaderItem(header)
        #self.datalayout.addWidget(self.tree)
        #path,fn = self.h5_
        
    def clear(self):
        self.tree.clear()
            
    def add_file(self,h5file):  

        self.h5_file_path = h5file
        self.f = h5py.File(h5file,'r')
        self.filename = self.f.filename.split('/')[-1]
        
        self.tree_root = QTreeWidgetItem(self.tree,[self.filename,self.h5_file_path,'/'])
        #self.tree.setColumnWidth(0,250)
        #self.tree.setColumnWidth(1,0)
        #self.tree.setColumnWidth(2,0)
        self.data_path = self.tree_root.text(2)
        
        self.add_branch(self.tree_root,self.f)
        self.tree.itemClicked.connect(self.onItemClicked)
        #self.setLayout(self.datalayout)
        #self.show()

    def add_branch(self,tree_root,h5file):
        for _ in h5file.keys():
            try:
                data_path = tree_root.text(2) + _ + '/'
                branch = QTreeWidgetItem([str(_),
                                          str(self.h5_file_path),
                                          data_path,
                                            ])
                tree_root.addChild(branch)
                #self.data_path = tree_root.text(2)
                #if  isinstance(h5file[_],h5py.Group):
                #    self.add_branch(branch,h5file[_])
            except Exception as e:
                print(e)
                pass 
    
    @pyqtSlot(QTreeWidgetItem,int)
    def onItemClicked(self,item):
        self.data_path = item.text(2)
        self.keywords = item.text(0)
        #p = self.h5_file_path.split('/id13/')[]
        #proposal_num  = p.split('/')[-1]
        #proposal_path = p + '/id13/'
        #print('\n\n',self.h5_file_path,
        #      proposal_path,'\n\n',proposal_num)
        
        #with h5py.File(self.h5_file_path,'r') as f:
        #    try:
        #        print(f[self.data_path])
        #    except Exception as e:
        #        print(e)
        #        pass

class titledTable():
    def __init__(self,title):
        self.title = QLabel(title)
        self.table = QTableWidget()
        self.table.setShowGrid(True)
        
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.title)
        self.layout.addWidget(self.table)
        
    def clear():
        self.table.setRowCount(0)
        self.table.setColumnCount(0)
        self.clear()
    
    def set_item(self,row,col,item):
        if isinstance(item, str):
            self.table.setItem(row, col, QTableWidgetItem(item))
        else:
            print("Type Error: Item Must Be a String")
    
    def num_cols(self,values):
        value_shape = np.shape(values)
        numcols = 1
        
        if len(value_shape) > 1:    
            numcols = value_shape[1]    
        
#if __name__ == '__main__':
#    app = QApplication(sys.argv)
#    Tree = tree()
#    Tree.show()
#    sys.exit(app.exec_())
