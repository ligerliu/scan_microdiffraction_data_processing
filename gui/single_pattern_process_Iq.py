import sys
import numpy as np
import pyFAI
import fabio
import h5py

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

from silx.gui.plot import *
from silx.gui.plot.PlotWidget import PlotWidget
from silx.gui.plot.PlotWindow import PlotWindow
from silx.gui import qt

from silx.gui.plot.PlotTools import PositionInfo
from silx.gui.plot.ImageView import ImageView
from silx.gui.colors import Colormap



'''
need to load the data and need to load the mask and poni for integration process
'''

class single_pattern_window(QWidget):
    def __init__(self):
        super().__init__()
        self.data  = None
        self.poni  = None
        self.mask  = None
        self.bkgd  = None
        self.mainwindow()

    def mainwindow(self):
        self.setGeometry(200,200,1000,520)
        self.setWindowTitle('process single pattern')
        
        self.file_load = QWidget(self)
        self.file_load.setGeometry(20,20,180,30)
        self.file_load_bttn = QPushButton('load data',self.file_load)
        self.file_load_bttn.clicked.connect(self.load_data)

        self.poni_load = QWidget(self)
        self.poni_load.setGeometry(20,60,180,30)
        self.poni_load_bttn = QPushButton('load poni',self.poni_load)
        self.poni_load_bttn.clicked.connect(self.load_poni)
        
        self.mask_load = QWidget(self)
        self.mask_load.setGeometry(20,100,180,30)
        self.mask_load_bttn = QPushButton('load mask',self.mask_load)
        self.mask_load_bttn.clicked.connect(self.load_mask)
        
        self.bkgd_load  = QWidget(self)
        self.bkgd_load.setGeometry(20,140,180,30)
        self.bkgd_load_bttn = QPushButton('load bkgd',self.bkgd_load)
        self.bkgd_load_bttn.clicked.connect(self.load_bkgd)
        
        self.bkgd_sub_box = QCheckBox('bkgd subtract',self)
        self.bkgd_sub_box.move(20,180)
       
        self.update_pattern = QWidget(self)
        self.update_pattern.setGeometry(20,220,180,30)
        self.update_pattern_bttn = QPushButton('show img',self.update_pattern)
        self.update_pattern_bttn.clicked.connect(self.show_pattern) 
        
        self.cal_qphi_pattern = QWidget(self)
        self.cal_qphi_pattern.setGeometry(20,260,180,30)
        self.cal_qphi_pattern_bttn = QPushButton('cal qphi',self.cal_qphi_pattern)
        self.cal_qphi_pattern_bttn.clicked.connect(self.cal_qphi) 
        
        self.cal_Iq_pattern = QWidget(self)
        self.cal_Iq_pattern.setGeometry(20,300,180,30)
        self.cal_Iq_pattern_bttn = QPushButton('cal Iq',self.cal_Iq_pattern)
        self.cal_Iq_pattern_bttn.clicked.connect(self.cal_Iq) 
        
        self.pttn_vmin_label = QLabel("vmin:",self)
        self.pttn_vmin_label.setGeometry(10,420,60,20)
        self.pttn_vmin_input = QLineEdit("0",self)
        self.pttn_vmin_input.setGeometry(80,420,120,20)
        self.pttn_vmax_label = QLabel("vmax:",self)
        self.pttn_vmax_label.setGeometry(10,450,60,20)
        self.pttn_vmax_input = QLineEdit("1",self)
        self.pttn_vmax_input.setGeometry(80,450,120,20)
        
        self.pyfai = None      
        self.show()
    
    def show_pattern(self):
        try:
            self.pattern_display = None
            self.pattern_display = Plot2D(self)
            self.pattern_display.clear()
            self.pattern_display.setGeometry(220,20,700,500)
            if self.bkgd_sub_box.isChecked():
                if (not isinstance(self.bkgd,type(None))) and \
                   (not isinstance(self.data,type(None))):
                    try:
                        img = self.data - self.bkgd
                    except Exception as e:
                        print(Exception)
                        pass
                else:   
                        img = self.data
            else:
                img = self.data
            
            self.pttn_vmin = float(self.pttn_vmin_input.text())
            self.pttn_vmax = float(self.pttn_vmax_input.text())
            colormap = Colormap(name = 'inferno', normalization = 'linear',
                                vmin = self.pttn_vmin, vmax = self.pttn_vmax)

            self.pattern_display.addImage(img,colormap=colormap)
            self.pattern_display.setYAxisInverted(flag=False)
            self.pattern_display.setKeepDataAspectRatio(flag=True)
            if not isinstance(self.pyfai,type(None)):
                toolBar = qt.QToolBar()
                self.pattern_display.addToolBar(qt.Qt.BottomToolBarArea,toolBar)
                position = PositionInfo(plot = self.pattern_display,
                    converters = [('Q', lambda x,y: self.pyfai.qArray()[int(y),int(x)]/10),
                                 ('Azimuth', lambda x,y: self.pyfai.chiArray()[int(y),int(x)])])
                toolBar.addWidget(position)
            
            self.pattern_display.show()
        except Exception as e:
            print(e)
            pass
    
    def load_data(self):
        try:
            self.data_fn,_ = QFileDialog.getOpenFileName(self,"load_data","","")
            self.data = fabio.open(self.data_fn).data
            #self.show_pattern()
        except Exception as e:
            print(e)
            pass
    
    def load_poni(self):
        try:
            self.poni,_ = QFileDialog.getOpenFileName(self,"load_poni","","")
            if not isinstance(self.poni,type(None)): 
                self.pyfai  = pyFAI.load(self.poni)
            else:
                print("poni is need for integration")
        except Exception as e:
            print(e)
            pass

    def load_mask(self):
        try:
            self.mask_fn,_ = QFileDialog.getOpenFileName(self,"load_mask","","")
            if 'npz' in self.mask_fn:
                self.mask = np.load(self.mask_fn)['mask']
            else:
                self.mask = fabio.open(self.mask_fn).data
        except Exception as e:
            print(e)
            pass
        
    def load_bkgd(self):
        try:
            self.bkgd_fn,_ = QFileDialog.getOpenFileName(self,"load_bkgd","","")
            self.bkgd = fabio.open(self.bkgd_fn).data
        except Exception as e:
            print(e)
            pass
    
    def cal_qphi(self):
        self.qphi_win = QWidget()
        self.qphi_win.setGeometry(400,400,870,500)
        self.qphi_q_label = QLabel("q_npts:",self.qphi_win)
        self.qphi_q_label.setGeometry(10,10,60,20)
        self.qphi_q_shape = QLineEdit("720",self.qphi_win)
        self.qphi_q_shape.setGeometry(80,10,120,20)
        self.qphi_a_label = QLabel("a_npts:",self.qphi_win)
        self.qphi_a_label.setGeometry(10,40,60,20)
        self.qphi_a_shape = QLineEdit("360",self.qphi_win)
        self.qphi_a_shape.setGeometry(80,40,120,20)
        
        
        self.qphi_vmin_label = QLabel("vmin:",self.qphi_win)
        self.qphi_vmin_label.setGeometry(10,120,60,20)
        self.qphi_vmin_input = QLineEdit("0",self.qphi_win)
        self.qphi_vmin_input.setGeometry(80,120,120,20)
        self.qphi_vmax_label = QLabel("vmax:",self.qphi_win)
        self.qphi_vmax_label.setGeometry(10,150,60,20)
        self.qphi_vmax_input = QLineEdit("1",self.qphi_win)
        self.qphi_vmax_input.setGeometry(80,150,120,20)
        
        qphi_process_bttn = QPushButton("qphi map",self.qphi_win)
        qphi_process_bttn.setGeometry(10,70,60,20)
        qphi_process_bttn.clicked.connect(self.show_qphi)
        
        self.sub_qphi_win = None 
        self.qphi_win.show()
    
    def show_qphi(self):
        try:
            if  isinstance(self.sub_qphi_win,type(None)):
                self.sub_qphi_win = Plot2D(self.qphi_win)
                self.sub_qphi_win.setGeometry(220,20,640,480)
            else:
                self.sub_qphi_win.clear()
            self.sub_qphi_win.show()
            self.q_npts = int(self.qphi_q_shape.text()) 
            self.a_npts = int(self.qphi_a_shape.text())
            qphi,self.q,self.a = self.calculate_Iqphi(q_npts = self.q_npts,
                                            a_npts = self.a_npts,
                                            )

            self.qphi_vmin = float(self.qphi_vmin_input.text())
            self.qphi_vmax = float(self.qphi_vmax_input.text())
            colormap = Colormap(name = 'gray', normalization = 'linear',
                                vmin = self.qphi_vmin, vmax = self.qphi_vmax)
            self.sub_qphi_win.addImage(qphi,colormap=colormap,
                                       origin=(0,-180),
                                       scale=((self.q[-1]-self.q[0])/len(self.q),
                                       (self.a[-1]-self.a[0])/len(self.a)))

        except Exception as e:
            print(e)
            pass
                
    def calculate_Iqphi(self,
                       method  = None,
                       q_npts=300,
                       a_npts=180,
                       q_unit="q_A^-1",
                       **kwargs):
        if isinstance(method,type(None)):
            method = 'csr'
        else:
            method = method
        if self.bkgd_sub_box.isChecked():
            if (not isinstance(self.bkgd,type(None))) & \
               (not isinstance(self.data,type(None))):
                try:
                    data = self.data - self.bkgd
                except Exception as e:
                    print(e)
                    pass
        else:
            data = self.data
        if hasattr(self,'pyfai'):
            if isinstance(self.mask,type(None)):
                qphi,q,a = self.pyfai.integrate2d(
                                        data,
                                        npt_rad  = q_npts,
                                        npt_azim = a_npts,
                                        unit     = q_unit,
                                        method   = method,
                                        **kwargs)
            else:
                try:
                    qphi,q,a = self.pyfai.integrate2d(
                                        data,
                                        npt_rad  = q_npts,
                                        npt_azim = a_npts,
                                        unit     = q_unit,
                                        method   = method,
                                        **kwargs)
                except Exception as e:
                    print(e)
                    pass
            return qphi,q,a
        else:
            print('poni is needed for pyfai calculation')
            pass



    def cal_Iq(self):
        self.Iq_win = QWidget()
        self.Iq_win.setGeometry(400,400,870,500)
        self.Iq_q_label = QLabel("q_npts:",self.Iq_win)
        self.Iq_q_label.setGeometry(10,10,60,20)
        self.Iq_q_shape = QLineEdit("720",self.Iq_win)
        self.Iq_q_shape.setGeometry(80,10,120,20)
        
        Iq_process_bttn = QPushButton("Iq cal",self.Iq_win)
        Iq_process_bttn.setGeometry(10,70,60,20)
        Iq_process_bttn.clicked.connect(self.show_Iq)
        
        self.keep_current_Iq = QCheckBox('keep curve',self.Iq_win)
        self.keep_current_Iq.move(10,100)

        self.Iq_hmsk_label = QLabel("hint mask:",self.Iq_win)
        self.Iq_hmsk_label.setGeometry(10,150,60,20)
        self.Iq_hmsk_shape = QLineEdit("65000",self.Iq_win)
        self.Iq_hmsk_shape.setGeometry(80,150,120,20)

        self.sub_Iq_win = None 
        self.Iq_win.show()
    
    def show_Iq(self):
        try:
            if  isinstance(self.sub_Iq_win,type(None)):
                self.sub_Iq_win = Plot1D(self.Iq_win)
                self.sub_Iq_win.setGeometry(220,20,640,480)
                self.Iq_old = []
            else:
                if not self.keep_current_Iq.isChecked():
                    self.sub_Iq_win.clear()
                    self.Iq_old = []
            self.q_npts = int(self.Iq_q_shape.text()) 
            self.hmsk = int(self.Iq_hmsk_shape.text())
            self.q,Iq = self.calculate_Iq(q_npts = self.q_npts,hmsk=self.hmsk)
            
            if self.keep_current_Iq.isChecked():
                self.Iq_old.append(Iq)
                for i in range(len(self.Iq_old)):
                    #print(self.Iq_old[i])
                    self.sub_Iq_win.addCurve(self.q,self.Iq_old[i],legend=f'{i}',xlabel='Q',ylabel='I (a.u.)')
            else:
                self.Iq_old = []
                self.sub_Iq_win.addCurve(self.q,Iq,replace=True)
            self.sub_Iq_win.show()

        except Exception as e:
            print(e)
            pass

    
    def calculate_Iq(self,
                     q_npts = 300,
                     hmsk   = 65000):
        if self.bkgd_sub_box.isChecked():
            if (not isinstance(self.bkgd,type(None))) & \
               (not isinstance(self.data,type(None))):
                try:
                    data = self.data - self.bkgd
                except Exception as e:
                    print(e)
                    pass
        else:
            data = self.data
        self.mask = (self.mask + (data>hmsk))
        if hasattr(self,'pyfai'):
            if isinstance(self.mask,type(None)):
                q,I  = self.pyfai.integrate1d(data,
                                              q_npts,
                                             )
                q /= 10
            else:
                try:
                    q,I = self.pyfai.integrate1d(data,
                                         q_npts,
                                         mask = self.mask)                    
                    q /= 10
                except Exception as e:
                    print(e)
                    pass
            return q,I
        else:
            print('poni is needed for pyfai calculation')
            pass
        return q,I

    
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = single_pattern_window()
    sys.exit(app.exec_()) 

