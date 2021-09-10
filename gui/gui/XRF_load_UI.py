import sys
sys.path.append('../xs_proc')
from auto_load_XRF import *

import h5py
import numpy as np
import fabio

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
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
    if vmax > np.nanmax(data):
        vmax = np.nanmax(data)
    colormap = Colormap(name='gray',
                normalization='linear',
                vmin = vmin,
                vmax = vmax
                )
    return colormap

class XRF_UI(QWidget):
    def __init__(self,obj):
        super().__init__()
        self.obj = obj
        self.XRF_UI()
        self.show()
    
    def XRF_UI(self):
        self.setGeometry(500,40,900,650)
        self.setWindowTitle('load XRF')
        
        energy_set = QWidget(self)
        energy_set.setGeometry(20,240,240,120)
        energy_max_label = QLabel("Emax:",energy_set)
        energy_max_label.setGeometry(10,10,60,20)
        self.energy_max_value = QLineEdit("40000",energy_set)
        self.energy_max_value.setGeometry(80,10,120,20)
        
        energy_min_label = QLabel("Emin:",energy_set)
        energy_min_label.setGeometry(10,40,60,20)
        self.energy_min_value = QLineEdit("0",energy_set)
        self.energy_min_value.setGeometry(80,40,120,20)
        
        process_button = QPushButton("load XRF",self)
        process_button.setGeometry(20,420,60,20)
        process_button.clicked.connect(self.load_XRF_roi)
    
    
    def load_XRF_roi(self):
        # currently still in form of ID13 customization, need to be generized
        try:
            t = time()
            self.Emax = int(self.energy_max_value.text())
            self.Emin = int(self.energy_min_value.text())
            print(self.Emin,self.Emax)
            self.xrf_map = auto_load_xrf(self.obj,
                                    element=None,
                                    energy_roi = (self.Emin,self.Emax)
                                    )

            print(np.max(self.xrf_map))
            self.xrf_roi_show = PlotWindow(self)#,position=True)
            colormap = setup_colormap(self.xrf_map)
            self.xrf_roi_show.addImage(self.xrf_map,colormap=colormap)

            toolBar = qt.QToolBar()
            self.xrf_roi_show.addToolBar(
            qt.Qt.BottomToolBarArea,
            toolBar)
            position = PositionInfo(plot= \
            self.xrf_roi_show,
            converters=[('X', 
            lambda x,y: int(np.round(x))),
            ('Y',
            lambda x,y: int(np.round(y)))])
            toolBar.addWidget(position)

            self.position = PositionInfo(plot=self.xrf_roi_show)
            #self.xrf_roi_show.sigPlotSignal.connect(self.roi_map_clicked)
            #self.xrf_roi_show.sigPlotSignal.connect(self.roi_map_polygon)
            self.xrf_roi_show.move(250,20)
            self.subwindow1 = None
            self.xrf_roi_show.show()
            print(time() - t)
        except Exception as e:
            print(e)
            pass
