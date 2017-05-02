######################################################################
# Author: Udaysankar Nair, Department of Atmospheric Science,        #
#         University of Alabama in Huntsville                        #
#                                                                    #
# Python class WrfDataset, which is used for manipulating WRF output #
# files.  This class developed as a part of ESS680 course in Spring  #
# of 2016.                                                           #
######################################################################


from __future__ import print_function

import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

import numpy as np
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backends import qt_compat
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import netCDF4
import PyGEOMET.utils.wrf_functions as wrf
import PyGEOMET.utils.SkewTobj as SkewT
import PyGEOMET.utils.DerivedVar as wrf_dvar
import PyGEOMET.datasets.GOESDataset as GOES
import PyGEOMET.datasets.WRFDataset as WrfDataset
import PyGEOMET.datasets.NCARDataset as NCAR
import PyGEOMET.datasets.MERRADataset as MERRA
import PyGEOMET.datasets.GOESClassDataset as GOESClass
import PyGEOMET.datasets.RadarDataset as NEXRAD
import PyGEOMET.datasets.CMAQDataset as CmaqDataset
import PyGEOMET.datasets.METDataset as METDataset
import PyGEOMET.utils.LayoutFormat as Layout
import time
import os
import warnings

#This is used to suppress plotting errors associated with bad lat/long data
warnings.filterwarnings('ignore', 'invalid*', RuntimeWarning)

###############################################################################
####                        Begin CanvasWidget() Object                    ####
###############################################################################

class CanvasWidget(QWidget):
    def __init__(self, parent, plotNumber):
        QWidget.__init__(self)
        self.fig = Figure((4,1.5), dpi=100)
        #self.fig = Figure((4,1.5), tight_layout())
        #self.fig = plt.figure(figsize=(8,4))
        self.setMinimumSize(500,500)
        self.canvas = FigureCanvas(self.fig)
        self.setMouseTracking(False)
        self.canvas.setParent(parent)
        self.canvas.setMinimumSize(self.canvas.size())
        boxTitleString = 'Plot # ' + str(plotNumber)
        self.gbox = QGroupBox(boxTitleString,parent)
        self.gbox.setStyleSheet(Layout.QGroupBox2())

        plotLayout = QVBoxLayout()
        plotLayout.addWidget(self.canvas)
        self.mpl_toolbar = NavigationToolbar(self.canvas, parent)
        plotLayout.addWidget(self.mpl_toolbar)

        self.gbox.setLayout(plotLayout)

        widgetLayout = QVBoxLayout()
        widgetLayout.addWidget(self.gbox)
        self.setLayout(widgetLayout)

    def setPlot(self,pltobj):
        pltobj.setFigure(self.fig)
        self.plotObj = pltobj
    
    def drawPlot(self):

        #self.plotObj.figure.subplots_adjust(left=0.05,right=0.95,
        #                                    bottom=0.05,top=0.95)
        #Horizontal Plot
        if self.plotObj.currentPType == 0:
            self.plotObj.plot()
            if self.plotObj.appobj.eom == True:
                if self.plotObj.cid is None:
                    self.plotObj.cid = self.canvas.mpl_connect('button_press_event', self.onclick)
                if self.plotObj.col is not None and self.plotObj.row is not None:
                    self.plotObj.getTable(self.plotObj.col,self.plotObj.row)
            
            #self.plotObj.figure.tight_layout()
           # self.plotObj.figure.subplots_adjust(left=0.08,right=0.92,
           #                                     bottom=0.08,top=0.92)
            self.canvas.draw()
        
        #Vertical Cross section
        if self.plotObj.currentPType == 1:
            self.plotObj.plotCS()
            #self.plotObj.figure.tight_layout()
            self.canvas.draw()

        #Skew T
        if self.plotObj.currentPType == 2:
            if self.plotObj.replot2d == True:
                self.plotObj.plot()
                #self.plotObj.figure.tight_layout()
                self.canvas.draw()
            if self.plotObj.cid is None: 
                self.plotObj.cid = self.canvas.mpl_connect('button_press_event', self.onclick)
            if self.plotObj.col is not None and self.plotObj.row is not None:
                self.plotObj.plotSkewT(self.plotObj.col,self.plotObj.row)
                plt.ion()
                plt.show()  

        #Vertical Profile
        if self.plotObj.currentPType == 3:
            if self.plotObj.replot2d == True:
                self.plotObj.plot()
                #self.plotObj.figure.tight_layout()
                self.canvas.draw()
            if self.plotObj.cid is None:
                self.plotObj.cid = self.canvas.mpl_connect('button_press_event', self.onclick)
            if self.plotObj.col is not None and self.plotObj.row is not None:
                self.plotObj.plotVertPro(self.plotObj.col,self.plotObj.row)
                plt.ion()
                plt.show()

        #Time Series
        if self.plotObj.currentPType == 4:
            if self.plotObj.replot2d == True:
                self.plotObj.plot()
                #self.plotObj.figure.tight_layout()
                self.canvas.draw() 
            if self.plotObj.cid is None:
                self.plotObj.cid = self.canvas.mpl_connect('button_press_event', self.onclick)
            if self.plotObj.col is not None and self.plotObj.row is not None:
                if self.plotObj.newvar == True:
                    self.plotObj.createTimeSeries(self.plotObj.col,self.plotObj.row)
                    self.plotObj.newvar = False
                self.plotObj.plotTimeSeries()
                plt.ion()
                plt.show()

        #Difference Plot
        if self.plotObj.currentPType == 5:
            self.plotObj.plotDifference()
            self.canvas.draw()

    def mouseMoveEvent(self, event):
        if event.buttons() == Qt.NoButton:
            print("Simple mouse motion",event.buttons())
        elif event.buttons() == Qt.LeftButton:
            print ("Left click drag")
        elif event.buttons() == Qt.RightButton:
            print("Right click drag")  

    def onclick(self,event):
        ix, iy = event.xdata, event.ydata
        lon, lat  = self.plotObj.dataSet.map[self.plotObj.currentGrid-1](
                    self.plotObj.dataSet.glons[self.plotObj.currentGrid-1],
                    self.plotObj.dataSet.glats[self.plotObj.currentGrid-1])
        clon,clat = self.plotObj.dataSet.map[self.plotObj.currentGrid-1](
                    self.plotObj.dataSet.lon0[self.plotObj.currentGrid-1],
                    self.plotObj.dataSet.lat0[self.plotObj.currentGrid-1])
        tmp = np.argsort(np.abs(lat[:,int(self.plotObj.dataSet.nx[self.plotObj.currentGrid -1]/2)] - clat))[0]
        tmp2 = np.argsort(np.abs(lon[tmp,:] - clon))[0]
        row = np.argsort(np.abs(lat[:,tmp2] - iy))[0]
        col  = np.argsort(np.abs(lon[row,:] - ix))[0]
        


        print( col, row )
        coords = [ix, iy]
       
        if len(coords) == 2:
            self.plotObj.col = col
            self.plotObj.row = row
            self.plotObj.replot2d = False
            self.plotObj.newvar = True
            self.drawPlot()
        else:
            self.plotObj.col = None
            self.plotObj.row = None
        
###############################################################################
####                         End CanvasWidget() Object                     ####
###############################################################################

###############################################################################
####                         Begin PlotSlab() Object                       ####
###############################################################################

class PlotSlab:   
    def __init__(self, fig=None,dset=None,AppWid=None):
        self.figure = None
        self.setFigure(fig)
        self.level = 0
        self.x1 = 0
        self.x2 = 0
        self.y1 = 0
        self.y2 = 0
        self.varTitle = 'Variable'
        self.x = None
        self.y = None
        self.zval = None
        #Index 1 defaults to first loaded dataset
        self.currentDset = len(dset)-1
        self.dSet = dset
        self.diffdata = self.dSet[-1]
        self.dataSet = self.dSet[-1]
        self.CBar = None
        self.pltFxn = None
        self.currentVar = None
        self.currentdVar = None
        self.derivedVar = False
        self.currentGrid = 1
        self.currentLevel = 0
        self.currentTime = 0
        self.currentPType = 0
        self.currentOrient = 0
        self.nz = None
        self.plotCount = 0
        self.CSPlotCount = 0
        self.selectLevel = None
        self.varUpdate = 1
        self.u10 = None
        self.v10 = None
        self.var = None
        self.var2 = None
        self.updateListView = True
        self.ref_pt = self.dataSet.lat0[self.currentGrid-1]
        self.minpress = 100.
        self.replot2d = True
        self.cid = None
        self.col = None
        self.row = None
        self.newvar = True
        self.colormin = None
        self.colormax = None
        self.ncontours = None
        self.colorbox = None
        self.vslicebox = None
        self.orientList = ['xz','yz','custom']
        self.appobj = AppWid
        self.appobj.filltype = "contourf"
        self.appobj.cmap = 'jet'
        self.appobj.plotbarbs = False
        self.appobj.plotvectors = False
        self.appobj.plotcontour2 = False
        self.vplot = None
        self.vertControl = None
        self.pControl = None
        self.diffControl = None
        self.selectedParcel = 0
        self.skewParcel = 'SB'
        self.appobj.background = None
        self.appobj.recallProjection = True
        self.appobj.changeColor = False
        self.ColorBar = None
        self.extend = 'both'
        self.appobj.cs = None
        self.appobj.cs2 = None
        self.appobj.cs2label = None
        self.appobj.barbs = None
        self.appobj.vectors = None
        self.appobj.vectorkey = None
        self.appobj.resolution = 'l'
        self.appobj.domain_average = None
        self.coasts = None
        self.countries = None
        self.states = None
        self.cmap = self.appobj.cmap
        self.appobj.eom = None
        self.colorlock = False
        #self.cPIndex = self.appobj.plotControlTabs.currentIndex()

    def setFigure(self, fig):
        if(fig != None):
            self.figure = fig

    def setConnection(self,fxn=None,plotnum=None):
        if(fxn != None):
            self.pltFxn = fxn
        if(plotnum != None):
            self.pNum = plotnum
           
    def readField(self):
        #print("read plot tab", self.appobj.plotControlTabs.currentIndex())
        #print("PIndex:",self.cPIndex)
        #print("Grid Num", self.dataSet.currentGrid-1)
        #if (self.appobj.plotControlTabs.currentIndex() != self.cPIndex):
        #    self.dataSet.updateData()
        #    print("in")
        if not self.derivedVar:
            if self.currentVar != None:
                #self.currentVar = np.where(np.array(self.dataSet.variableList) == "T2")[0][0]
                #self.nz = 1
                self.var = self.dataSet.readNCVariable(self.dataSet.variableList[self.currentVar],
                    barbs=self.appobj.plotbarbs, vectors = self.appobj.plotvectors,
                    contour2=self.appobj.plotcontour2)
                self.varTitle = self.dataSet.description +' ('+self.dataSet.units
                self.varTitle = self.varTitle +') \n'+self.dataSet.getTime()
                #print(self.varTitle)
                #print(self.var.shape)
                #account for unstaggering of the grids in 3-dimensional variables
                if self.appobj.dname == 'WRF' or self.appobj.dname == 'MET':
                    if len(self.var.shape) == 3:
                        if self.var.shape[0] == self.dataSet.nz[self.dataSet.currentGrid-1]:
                            self.var = wrf.unstaggerZ(self.var)
                        if self.var.shape[1] == self.dataSet.ny[self.dataSet.currentGrid-1]:
                            self.var = wrf.unstaggerY(self.var)
                        if self.var.shape[2] == self. dataSet.nx[self.dataSet.currentGrid-1]:
                            self.var = wrf.unstaggerX(self.var) 
                        self.varTitle = self.varTitle + ', Level=' + str(self.dataSet.levelList[self.currentLevel])
                    else:
                        if self.var.shape[0] == self.dataSet.ny[self.dataSet.currentGrid-1]:
                            plotself.var = wrf.unstaggerY(self.var)
                        if self.var.shape[1] == self.dataSet.nx[self.dataSet.currentGrid-1]:
                            self.var = wrf.unstaggerX(self.var)
                self.diffvar = self.var
            else:
                pass
        else:
            if self.currentdVar != None:
                t0 = time.clock()
                t1 = time.time()
                dvar = wrf_dvar.WRFDerivedVar(dset = self.dataSet, 
                                              var = self.dataSet.dvarlist[self.currentdVar],
                                              ptype = self.currentPType)
                print( time.clock() - t0,"seconds process time plots")
                print( time.time() - t1,"seconds wall time plots")

                self.var = dvar.var
                self.var2 = dvar.var2
                self.varTitle = dvar.varTitle
                if len(self.var.shape) == 3:
                    self.varTitle = self.varTitle + ', Level=' + str(self.dataSet.levelList[self.currentLevel])
                if len(self.var.shape) == 3:
                    self.u10 = dvar.u10[self.currentLevel]
                    self.v10 = dvar.v10[self.currentLevel]
                    self.var2 = dvar.var2[self.currentLevel]
                else:
                    self.u10 = dvar.u10
                    self.v10 = dvar.v10
                    self.var2 = dvar.var2

                self.diffvar = self.var
            else:
                pass
        #self.cPIndex = self.appobj.plotControlTabs.currentIndex() 
        #self.cPIndex = 1 
        #print("After pindex:", self.cPIndex)
 
    def readDiffField(self):

        if not self.derivedVar:
            self.diffvar = self.diffdata.readNCVariable(self.diffdata.variableList[self.currentVar],
                barbs=self.appobj.plotbarbs, vectors = self.appobj.plotvectors,
                contour2=self.appobj.plotcontour2)
            #self.varTitle = self.dataSet.description +' ('+self.dataSet.units
            #self.varTitle = self.varTitle +') \n'+self.dataSet.getTime()
            #print(self.varTitle)
            #print(self.var.shape)
            #account for unstaggering of the grids in 3-dimensional variables
            if self.appobj.dname == 'WRF' or self.appobj.dname == 'MET':
                if len(self.diffvar.shape) == 3:
                    if self.diffvar.shape[0] == self.diffdata.nz[self.diffdata.currentGrid-1]:
                        self.diffvar = wrf.unstaggerZ(self.diffvar)
                    if self.diffvar.shape[1] == self.diffdata.ny[self.diffdata.currentGrid-1]:
                        self.diffvar = wrf.unstaggerY(self.diffvar)
                    if self.diffvar.shape[2] == self. diffdata.nx[self.diffdata.currentGrid-1]:
                        self.diffvar = wrf.unstaggerX(self.diffvar)
                    #self.varTitle = self.varTitle + ', Level=' + str(self.dataSet.levelList[self.currentLevel])
                else:
                    if self.diffvar.shape[0] == self.diffdata.ny[self.diffdata.currentGrid-1]:
                        self.diffvar = wrf.unstaggerY(self.diffvar)
                    if self.diffvar.shape[1] == self.diffdata.nx[self.diffdata.currentGrid-1]:
                        self.diffvar = wrf.unstaggerX(self.diffvar)
        else:
            #Call wrf_dvar class
            #t0 = time.clock()
            #t1 = time.time()
            dvar = wrf_dvar.WRFDerivedVar(dset = self.diffdata,
                                            var = self.dataSet.dvarlist[self.currentdVar],
                                            ptype = self.currentPType)

            self.diffvar = dvar.var

        self.varTitle = self.varTitle +'\n'+self.vartitle_ext

    def getControlBar(self):

        self.appobj.cbar = self.dataSet.getDsetControlBar(self)
        self.dataSet.qscroll.setStyleSheet(Layout.QScrollArea())
        self.optionTabs = QTabWidget()
        self.optionTabs.setMaximumHeight(self.appobj.screeny*.8*.4)        
        self.appobj.tabLayout.addWidget(self.appobj.cbar)
 
    def lockColorBar(self):
        self.colorlock = True
  
    def unlockColorBar(self):
        self.colorlock = False

    def controlColorBar(self):

        boxTitleString = 'Colorbar Control'
        self.colorbox = QGroupBox()
        self.colorbox.setStyleSheet(Layout.QGroupBox())

        #font = QFont();
        #font.setBold(True);
        #self.colorbox.setFont(font);
        colorboxLayout = QVBoxLayout()
        self.colorbox.setLayout(colorboxLayout)

        radio = QWidget()
        radiolayout = QHBoxLayout()
        radio.setLayout(radiolayout)
        selectLock = QRadioButton("Lock Color Range")
        selectLock.setStyleSheet(Layout.QCheckBox())
        selectLock.clicked.connect(lambda:self.lockColorBar())
        unlock = QRadioButton("Unlock Color Range")
        unlock.setStyleSheet(Layout.QCheckBox())
        unlock.clicked.connect(lambda:self.unlockColorBar())
        unlock.setChecked(True)
        radiolayout.addWidget(selectLock)
        radiolayout.addWidget(unlock)        

        #Maximum
        maxLabel = QLabel()
        maxLabel.setText('Maximum Value')
        self.selectMax = QLineEdit()
        self.selectMax.setStyleSheet(Layout.QLineEdit())
        self.selectMax.setText(str(self.colormax))
        #self.selectMax.editingFinished.connect(self.enterCValue)

        #Minimum
        minLabel = QLabel()
        minLabel.setText('Minimum Value')
        self.selectMin = QLineEdit()
        self.selectMin.setStyleSheet(Layout.QLineEdit())
        self.selectMin.setText(str(self.colormin))
        #self.selectMin.editingFinished.connect(self.enterCValue)   

        #Number of contours
        contoursLabel = QLabel()
        contoursLabel.setText('Number of Contours')
        self.selectcontours = QLineEdit()
        self.selectcontours.setStyleSheet(Layout.QLineEdit())
        self.selectcontours.setText(str(self.ncontours))
        #self.selectcontours.editingFinished.connect(self.enterCValue)     

        #Replot button
        self.PlotButton = QPushButton('Replot')
        self.PlotButton.setStyleSheet(Layout.QPushButton3())
        self.PlotButton.resize(self.PlotButton.minimumSizeHint())
        self.PlotButton.clicked.connect(self.enterCValue)

        #Connect widgets
        colorboxLayout.addWidget(radio)
        colorboxLayout.addWidget(maxLabel)
        colorboxLayout.addWidget(self.selectMax)
        colorboxLayout.addWidget(minLabel)
        colorboxLayout.addWidget(self.selectMin)
        colorboxLayout.addWidget(contoursLabel)
        colorboxLayout.addWidget(self.selectcontours)
        colorboxLayout.addWidget(self.PlotButton)
        #self.qscrollLayout.addWidget(self.colorbox)
        self.optionTabs.addTab(self.colorbox,boxTitleString)
        self.optionTabs.setCurrentIndex(self.optionTabs.count()-1)
        self.tabbingLayout.addWidget(self.optionTabs)

    def verticalSliceControl(self):
        boxTitleString = 'Vertical Slice Control'
        self.vslicebox = QGroupBox()      
        self.vslicebox.setStyleSheet(Layout.QGroupBox())

        vsliceboxLayout = QVBoxLayout()
        self.vslicebox.setLayout(vsliceboxLayout)
        orientLabel = QLabel()
        orientLabel.setText('Cross-Section Orientation')
        selectorient = QComboBox()
        selectorient.setStyleSheet(Layout.QComboBox())
        selectorient.addItems(self.orientList)
        selectorient.currentIndexChanged.connect(self.selectionChangeOrient)
        self.refboxLabel = QLabel()
        self.refboxLabel.setText('Lat:')
        self.refbox = QLineEdit()
        self.refbox.setStyleSheet(Layout.QLineEdit())
        self.refbox.setText(str(self.ref_pt))
        self.refbox.editingFinished.connect(self.enterPress)
        pboxLabel = QLabel()
        pboxLabel.setText('PTOP (hPa)')
        self.pbox = QLineEdit()
        self.pbox.setStyleSheet(Layout.QLineEdit())
        self.pbox.setText(str(self.minpress))
        self.pbox.editingFinished.connect(self.enterPress)
        vsliceboxLayout.addWidget(orientLabel)
        vsliceboxLayout.addWidget(selectorient)
        vsliceboxLayout.addWidget(self.refboxLabel)
        vsliceboxLayout.addWidget(self.refbox)
        vsliceboxLayout.addWidget(pboxLabel)
        vsliceboxLayout.addWidget(self.pbox)
        self.qscrollLayout.addWidget(self.vslicebox)
        self.optionTabs.addTab(self.vslicebox,boxTitleString)
        self.optionTabs.setCurrentIndex(self.optionTabs.count()-1)
        self.tabbingLayout.addWidget(self.optionTabs)
 
    def plot(self):
        #print(self.appobj.eom)
        if self.appobj.changeColor == True:
            self.cmap = self.appobj.cmap      

        #Set background map
        alpha = 1
        if self.appobj.clear:
            self.figure.clear()
            self.appobj.clear = False
            self.appobj.cs = None
            self.coasts = None
            self.states = None
            self.countries = None
            self.ColorBar = None
            self.appobj.domain_average = None
        if self.appobj.background is not None:
            pre = 'self.dataSet.map['
            gridnum = str(self.currentGrid-1)+']'
            end = self.appobj.background+'(ax=self.appobj.axes1['+str(self.pNum-1)+'])'
            exec(pre+gridnum+end,None,None)
            #Make variable transparent
            alpha=0.75

        if self.appobj.recallProjection == True:
            #print("axes", len(self.appobj.axes1))
            #print("Pnum",self.pNum)
            #if len(self.appobj.axes1) >= self.pNum:
            #print("plot tab", self.appobj.plotControlTabs.currentIndex())
            if len(self.appobj.axes1) >= self.appobj.plotCount:
                 self.appobj.axes1[self.pNum-1] = self.figure.add_subplot(111)
                 self.dataSet.resolution = self.appobj.resolution
                 self.dataSet.setProjection(self.currentGrid,axs=self.appobj.axes1[self.pNum-1])

            else:
                 self.appobj.axes1.append(self.figure.add_subplot(111))
                 self.dataSet.resolution = self.appobj.resolution
                 self.dataSet.setProjection(self.currentGrid,axs=self.appobj.axes1[self.pNum-1])
            #print("axes-after", len(self.appobj.axes1))
        xb, yb = self.dataSet.map[self.currentGrid-1](
                 self.dataSet.glons[self.currentGrid-1],
                 self.dataSet.glats[self.currentGrid-1])

        #Set NEXRAD settings
        if self.dataSet.dsetname == 'NEXRAD Radar':
            self.nz = 1
            self.appobj.filltype = 'pcolormesh'
            self.cmap = self.dataSet.cmap
            self.colormin = self.dataSet.range[0]
            self.colormax = self.dataSet.range[1]

        if(self.currentVar != None or self.currentdVar != None):
            if(self.nz == 1):
                pltfld = self.var
            else:
                pltfld = self.var[self.currentLevel]
            davg = np.nanmean(pltfld)
            time = self.dataSet.timeObj.strftime("%Y%m%d_%H%M%S")
            if not self.derivedVar:
                varname = self.dataSet.variableList[self.currentVar]
            else:
                varname = self.dataSet.dvarlist[self.currentdVar]
            #change the default plot name
            self.figure.canvas.get_default_filename = lambda: (varname + '_' + time + '.png')

            #Initialize colorbar range and increment
            if self.colormin is None:
                #self.colormin = np.floor(np.nanmin(pltfld))
                self.colormin = np.nanmin(pltfld)
                if self.colormin != self.colormin:
                    self.colormin = 0.
            if self.colormax is None:
                #self.colormax = np.ceil(np.nanmax(pltfld))
                self.colormax = np.nanmax(pltfld)
                if self.colormax != self.colormax:
                    self.colormax = 1.

            if self.ncontours is None:
                self.ncontours = 11.
            if self.colormax == self.colormin:
               if self.colormax == 0:
                   self.colormax = 1.
               else:
                   self.colormax = np.abs(self.colormax) * 2
               
            if self.colorbox is None:
                self.controlColorBar()

            if self.appobj.cs2 != None:
                for coll in self.appobj.cs2.collections:
                    coll.remove()
                for labels in self.appobj.cs2label:
                    labels.remove()
            if self.appobj.barbs != None:
                self.appobj.barbs.remove()
                self.appobj.vectors2.remove()
            if self.appobj.vectors != None:
                self.appobj.vectors.remove()
            if self.appobj.vectorkey != None:
                self.appobj.vectorkey.remove()
            if self.coasts != None and self.countries != None and self.states != None:
                self.coasts.remove()
                self.countries.remove()
                self.states.remove()

            if self.appobj.filltype == "contourf":
                if self.appobj.cs is not None:
                    for coll in self.appobj.cs.collections:
                        coll.remove()
                #print(self.dataSet.glons[self.currentGrid-1].shape)
                #print(pltfld.shape)
                self.appobj.cs = self.dataSet.map[self.currentGrid-1].contourf(
                                  self.dataSet.glons[self.currentGrid-1],
                                  self.dataSet.glats[self.currentGrid-1],
                                  pltfld,
                                  levels=np.linspace(self.colormin,self.colormax,self.ncontours),
                                  latlon=True,extend=self.extend,
                                  cmap=plt.cm.get_cmap(str(self.cmap)),
                                  alpha=alpha, ax=self.appobj.axes1[self.pNum-1])
                divider = make_axes_locatable(self.appobj.axes1[self.pNum-1])
                cax = divider.append_axes("right", size="5%", pad=0.1)
                if self.ColorBar != None:
                    self.ColorBar.remove()
                self.ColorBar = self.figure.colorbar(self.appobj.cs,cax=cax, extend = self.extend)
            elif self.appobj.filltype == "pcolormesh":
                lvls = np.linspace(self.colormin,self.colormax,self.ncontours)
                pltfld = np.ma.masked_less_equal(pltfld,self.colormin)
                norm = matplotlib.colors.Normalize(vmin=np.amin(lvls),vmax=np.amax(lvls))
                #if (len(np.where(pltfld <= self.colormin)[0]) > 1):
                if (self.derivedVar == True):
                    if (self.dataSet.dvarlist[self.currentdVar] == 'refl' or
                        np.ma.count_masked(pltfld).any()):
                        self.appobj.axes1[self.pNum-1] = None
                        self.appobj.domain_average = None
                        self.ColorBar = None
                        self.figure.clear()
                        if len(self.appobj.axes1) >= self.pNum:
                            self.appobj.axes1[self.pNum-1] = self.figure.add_subplot(111)
                            self.dataSet.resolution = self.appobj.resolution
                            self.dataSet.setProjection(self.currentGrid,axs=self.appobj.axes1[self.pNum-1])
                        else:
                            self.appobj.axes1.append(self.figure.add_subplot(111))
                            self.dataSet.resolution = self.appobj.resolution
                            self.dataSet.setProjection(self.currentGrid,axs=self.appobj.axes1[self.pNum-1])
                    xb, yb = self.dataSet.map[self.currentGrid-1](
                                              self.dataSet.glons[self.currentGrid-1],
                                              self.dataSet.glats[self.currentGrid-1])
                else:
                    if (self.dataSet.variableList[self.currentVar] == 'REFL_10CM' or
                        self.dataSet.dsetname == 'NEXRAD Radar' or 
                        np.ma.count_masked(pltfld).any()):
                        self.appobj.axes1[self.pNum-1] = None
                        self.ColorBar = None
                        self.figure.clear()
                        if len(self.appobj.axes1) >= self.pNum:
                            self.appobj.axes1[self.pNum-1] = self.figure.add_subplot(111)
                            self.dataSet.resolution = self.appobj.resolution
                            self.dataSet.setProjection(self.currentGrid,axs=self.appobj.axes1[self.pNum-1])
                        else:
                            self.appobj.axes1.append(self.figure.add_subplot(111))
                            self.dataSet.resolution = self.appobj.resolution
                            self.dataSet.setProjection(self.currentGrid,axs=self.appobj.axes1[self.pNum-1])
                    xb, yb = self.dataSet.map[self.currentGrid-1](
                                              self.dataSet.glons[self.currentGrid-1],
                                              self.dataSet.glats[self.currentGrid-1]) 
                self.appobj.cs = self.dataSet.map[self.currentGrid-1].pcolormesh(
                                  self.dataSet.glons[self.currentGrid-1],
                                  self.dataSet.glats[self.currentGrid-1],
                                  pltfld,
                                  norm=norm,
                                  latlon=True,cmap=plt.cm.get_cmap(str(self.cmap)),
                                  alpha=alpha, ax=self.appobj.axes1[self.pNum-1])
                #print(self.appobj.cs.QuadMesh)
                divider = make_axes_locatable(self.appobj.axes1[self.pNum-1])
                cax = divider.append_axes("right", size="5%", pad=0.1)
                if self.ColorBar != None:
                    self.ColorBar.remove()                
                self.ColorBar = self.figure.colorbar(self.appobj.cs,cax=cax)
            self.ColorBar.ax.tick_params(labelsize=9)
            self.appobj.axes1[self.pNum-1].set_title(self.varTitle,fontsize = 10)
            self.displayValuesXYZ(pltfld)

            if self.appobj.domain_average != None:
                self.appobj.domain_average.remove()
            #self.appobj.domain_average = self.appobj.axes1[self.pNum-1].text(0.95, -0.08,
            #     ("Domain Average: " + str(davg)),
            #     verticalalignment='bottom',horizontalalignment='right',
            #     transform = self.appobj.axes1[self.pNum-1].transAxes,
            #     color='k',fontsize=12)

            if self.appobj.plotbarbs == True or self.appobj.plotvectors == True:
                self.readField()

                if not self.derivedVar:
                    if self.nz != 1:
                        self.u10 = self.dataSet.u10[self.currentLevel]
                        self.v10 = self.dataSet.v10[self.currentLevel]
                    else:
                        self.u10 = self.dataSet.u10
                        self.v10 = self.dataSet.v10

                if (self.dataSet.dsetname == "MERRA" or self.dataSet.dsetname == 'NCEP/NCAR Reanalysis II'):
                    if len(self.u10.shape) == 3:
                        self.u10 = self.dataSet.u10[self.currentLevel]
                    if len(self.v10.shape) == 3:
                        self.v10 = self.dataSet.v10[self.currentLevel]
                    if self.dataSet.dsetname == 'MERRA':
                        interval = 20
                    if self.dataSet.dsetname == 'NCEP/NCAR Reanalysis II':
                        interval = 8

                elif (self.dataSet.dx[self.currentGrid-1] >= 15000):
                    if (self.dataSet.ny[self.currentGrid-1] < 250. and self.dataSet.nx[self.currentGrid-1] < 250.):
                        interval = 3
                    else:
                        interval = 5
                else:
                    if (self.dataSet.ny[self.currentGrid-1] < 250. and self.dataSet.nx[self.currentGrid-1] < 250.):
                        interval = 5
                    else:
                        interval = 10
                maxspeed = np.amax((self.u10[::interval,::interval]**2+self.v10[::interval,::interval]**2)**(1./2.))
                if self.appobj.plotbarbs == True:
                    self.appobj.barbs = self.appobj.axes1[self.pNum-1].barbs(
                        xb[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval], 
                        yb[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval], 
                        self.u10[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval], 
                        self.v10[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval], 
                        length = 5, sizes={"emptybarb":0.0}, barbcolor='k', flagcolor='k',linewidth=0.50, pivot='middle')
 
                    unorm = self.u10/np.sqrt(np.power(self.u10,2) + np.power(self.v10,2))
                    vnorm = self.v10/np.sqrt(np.power(self.u10,2) + np.power(self.v10,2))
                    self.appobj.vectors2 = self.appobj.axes1[self.pNum-1].quiver(xb[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                        yb[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                        unorm[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                        vnorm[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                        pivot='middle', headwidth=0, headlength=0, headaxislength=0,scale=60,width=0.0015)

                if self.appobj.plotvectors == True:
                        self.appobj.vectors = self.appobj.axes1[self.pNum-1].quiver(
                            xb[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval], 
                            yb[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                            self.u10[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                            self.v10[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval], 
                            color='k',units='inches',pivot='mid',scale=maxspeed*4)
                        self.appobj.vectorkey = self.appobj.axes1[self.pNum-1].quiverkey(self.appobj.vectors, 0.05, -0.08, 
                                                int(maxspeed*.75), str(int(maxspeed*.75))+' $m s^{-1}$',labelpos='E')
#                    else:
#                        maxspeed = np.amax((self.u10[points]**2+self.v10[points]**2)**(1./2.))
#                        print(maxspeed)
#                        self.appobj.vectors = self.appobj.axes1[self.pNum-1].quiver(
#                            xb[points], yb[points], self.u10[points],
#                            self.v10[points], color='k',units='inches',pivot='mid',scale=maxspeed*4)
#                        self.appobj.vectorkey = self.appobj.axes1[self.pNum-1].quiverkey(self.appobj.vectors, 0.05, -0.08, 
#                                                int(maxspeed*.75), str(int(maxspeed*.75))+' $m s^{-1}$',labelpos='E')

            if self.appobj.plotcontour2 == True:
                self.readField()
                if not self.derivedVar:
                    self.var2 = pltfld
                self.appobj.cs2 = self.appobj.axes1[self.pNum-1].contour(
                                xb, yb, self.var2,
                                colors='k',linewidths=1.5)
                if np.abs(self.var2).max() <= 1:
                    fmt = '%8.6f'
                else:
                    fmt = '%d'
                self.appobj.cs2label = self.appobj.cs2.clabel(inline=1, 
                                       font=10, fmt=fmt)
        #print(self.currentGrid,self.pNum)
        self.coasts = self.dataSet.map[self.currentGrid-1].drawcoastlines(ax=self.appobj.axes1[self.pNum-1])
        self.countries = self.dataSet.map[self.currentGrid-1].drawcountries(ax=self.appobj.axes1[self.pNum-1])
        self.states = self.dataSet.map[self.currentGrid-1].drawstates(ax=self.appobj.axes1[self.pNum-1])
        # draw parallels
        if self.dataSet.projectionType == "robin" or self.dataSet.projectionType == "geos" or \
           self.dataSet.projectionType == "ortho" or self.dataSet.projectionType == "aeqd":
            parallels = np.arange(-90.,90.,15)
            meridians = np.arange(0.,360., 30)
        else:
            plim = max(int((np.nanmax(np.array(np.abs(self.dataSet.glats[self.currentGrid-1])))
                           -np.nanmin(np.array(np.abs(self.dataSet.glats[self.currentGrid-1]))))/10.),1) 
            parallels = np.arange(-90.,90.,plim)
            mlim = max(int((np.nanmax(np.array(np.abs(self.dataSet.glons[self.currentGrid-1])))
                           -np.nanmin(np.array(np.abs(self.dataSet.glons[self.currentGrid-1]))))/10.),1)
            meridians = np.arange(0.,360.,mlim)
        
        #Can't put labels on Geostationary, Orthographic or Azimuthal Equidistant basemaps 	
        if self.dataSet.projectionType == "ortho" or self.dataSet.projectionType == "aeqd" or self.dataSet.projectionType == "geos":
            self.dataSet.map[self.currentGrid-1].drawparallels(parallels,ax=self.appobj.axes1[self.pNum-1])
            # draw meridians
            self.dataSet.map[self.currentGrid-1].drawmeridians(meridians,ax=self.appobj.axes1[self.pNum-1])
        else:	
            self.dataSet.map[self.currentGrid-1].drawparallels(parallels,labels=[1,0,0,0],fontsize=6,ax=self.appobj.axes1[self.pNum-1])
            # draw meridians
            self.dataSet.map[self.currentGrid-1].drawmeridians(meridians,labels=[0,0,0,1],fontsize=6,ax=self.appobj.axes1[self.pNum-1])

        if self.dataSet.dsetname == 'NEXRAD Radar':
            self.dataSet.map[self.currentGrid-1].plot(self.dataSet.lon0[0],
                    self.dataSet.lat0[0], marker='o',markersize=4,color='black',latlon=True, ax = self.appobj.axes1[self.pNum-1])
        self.appobj.recallProjection = False
        self.appobj.changeColor = False
        line, = self.appobj.axes1[self.pNum-1].plot([0], [0])  # empty line
        self.dataSelector = DataSelector(line)
 
    def plotCS(self):
        if self.appobj.changeColor == True:
            self.cmap = self.appobj.cmap
        if self.appobj.recallProjection == True:
            if len(self.appobj.axes1) >= self.pNum:
                 self.appobj.axes1[self.pNum-1] = self.figure.add_subplot(111)
            else:
                 self.appobj.axes1.append(self.figure.add_subplot(111))
        else:
            #self.appobj.axes1.remove(self.appobj.axes1[self.pNum-1])
            self.appobj.axes1[self.pNum-1] = None
            self.figure.clear()
            #self.appobj.axes1.append(self.figure.add_subplot(111))
            self.appobj.axes1[self.pNum-1] = self.figure.add_subplot(111)

        if(self.nz < 3):
            self.error3DVar()
        else:
            if (self.dataSet.dsetname == 'CMAQ'):
                var = self.var
                dims = var.shape 
                if self.orientList[self.currentOrient] == 'xz':
                    diff = abs(self.dataSet.glats[self.dataSet.currentGrid-1][:,0] - self.ref_pt)
                    ind = np.argsort(diff)
                    horiz = np.tile(np.sum(self.dataSet.glons[self.dataSet.currentGrid-1][ind[0:4],:],axis=0)/4.,(dims[0],1))
                    horiz2 = np.tile(np.sum(self.dataSet.glats[self.dataSet.currentGrid-1][ind[0:4],:],axis=0)/4.,(dims[0],1))
                    pvar = np.squeeze(np.sum(var[:,ind[0:4],:],axis=1)/4.)
                    plevs = np.swapaxes(np.tile(np.linspace(1,dims[0],dims[0]),(dims[2],1)),0,1)

                if self.orientList[self.currentOrient] == 'yz':
                    diff = abs(self.dataSet.glons[self.dataSet.currentGrid-1][0,:]-self.ref_pt)
                    ind = np.argsort(diff)
                    horiz = np.tile(np.sum(self.dataSet.glats[self.dataSet.currentGrid-1][:,ind[0:4]],axis=1)/4.,(dims[0],1))
                    horiz2 =np.tile(np.sum(self.dataSet.glons[self.dataSet.currentGrid-1][:,ind[0:4]],axis=1)/4.,(dims[0],1))
                    pvar = np.squeeze(np.sum(var[:,:,ind[0:4]],axis=2)/4.)
                    plevs = np.swapaxes(np.tile(np.linspace(1,dims[0],dims[0]),(dims[1],1)),0,1)               
 
                #Initialize colorbar range and increment
                if self.colormin is None:
                    #self.colormin = np.floor(pvar.min())
                    self.colormin = pvar.min()
                    if self.colormin != self.colormin:
                        self.colormin = 0.
                if self.colormax is None:
                    #self.colormax = np.ceil(pvar.max())
                    self.colormax = pvar.max()
                    if self.colormax != self.colormax:
                        self.colormax = 1.
                if self.ncontours is None:
                    self.ncontours = 11.
                if self.colormax == self.colormin:
                    if self.colormax == 0:
                        self.colormax = 1.
                    else:
                        self.colormax = self.colormax * 2.

                if self.colorbox is None:
                    self.controlColorBar()

                if self.vslicebox is None and self.currentPType == 1:
                    self.verticalSliceControl()

                if self.appobj.filltype == "contourf":
                    self.appobj.cs = self.appobj.axes1[self.pNum-1].contourf(horiz,
                          plevs, pvar,
                          levels=np.linspace(self.colormin,self.colormax,
                          self.ncontours), extend=self.extend,
                          cmap=plt.cm.get_cmap(str(self.cmap)))
                elif self.appobj.filltype == "pcolormesh":
                    lvls = np.linspace(self.colormin,self.colormax,self.ncontours)
                    pvar = np.ma.masked_less_equal(pvar,self.colormin)
                    norm = matplotlib.colors.Normalize(vmin=np.amin(lvls),vmax=np.amax(lvls))
                    if (self.dataSet.variableList[self.currentVar] == 'REFL_10CM' or
                        np.ma.count_masked(pvar).any()):
                        self.appobj.axes1[self.pNum-1] = None
                        self.ColorBar = None
                        self.figure.clear()
                    self.appobj.cs = self.appobj.axes1[self.pNum-1].pcolormesh(
                                      horiz, plevs, pvar, norm=norm,
                                      cmap=plt.cm.get_cmap(str(self.cmap)))

                if self.ColorBar != None:
                    self.ColorBar.remove()
                self.Colorbar = self.figure.colorbar(self.appobj.cs, ax=self.appobj.axes1[self.pNum-1], 
                                                     orientation='horizontal', pad=0.06)
                self.Colorbar.ax.tick_params(labelsize=9)
                ax1 = self.appobj.axes1[self.pNum-1]
                ax1.set_ylabel("Model Level")
                ax1.set_ylim(1,dims[0])
                ax1.set_title(self.varTitle,fontsize = 10)

                axins = inset_axes(ax1,width="30%",height="30%",loc=1,borderpad=0)
                map2 = self.dataSet.map[self.currentGrid-1]
                map2.ax=axins
                map2.etopo()
                map2.drawcoastlines()
                map2.drawcountries()
                map2.drawstates()
                # draw parallels.
                plim = max(int((abs(self.dataSet.glats[self.currentGrid-1]).max()-abs(self.dataSet.glats[self.currentGrid-1]).min())/2.5),0.5)
                parallels = np.arange(-90.,90.,plim)
                map2.drawparallels(parallels,labels=[0,1,0,0],fontsize=6)
                # draw meridians
                mlim = max(int((abs(self.dataSet.glons[self.currentGrid-1]).max()-abs(self.dataSet.glons[self.currentGrid-1]).min())/2.5),0.5)
                meridians = np.arange(0.,360.,mlim)
                map2.drawmeridians(meridians,labels=[0,0,1,0],fontsize=6)
                if self.orientList[self.currentOrient] == 'xz':
                    xx, yy = map2(horiz[0,:], horiz2[0,:])
                if self.orientList[self.currentOrient] == 'yz':
                    xx, yy = map2(horiz2[0,:],horiz[0,:])

                map2.plot(xx,yy,color='red',linewidth=2)
                self.appobj.recallProjection = False
                #self.displayValuesXYZ(pvar)
            else: 
                var = self.var
                self.u10, self.v10, w, press, height = self.dataSet.getVertVars()
                var2 = w
                dims = var.shape
                if self.orientList[self.currentOrient] == 'xz':
                    diff = abs(self.dataSet.glats[self.dataSet.currentGrid-1][:,0] - self.ref_pt)
                    ind = np.argsort(diff)
                    horiz = np.tile(np.sum(self.dataSet.glons[self.dataSet.currentGrid-1][ind[0:4],:],axis=0)/4.,(dims[0],1))
                    horiz2 =np.tile(np.sum(self.dataSet.glats[self.dataSet.currentGrid-1][ind[0:4],:],axis=0)/4.,(dims[0],1))
                    pvar = np.squeeze(np.sum(var[:,ind[0:4],:],axis=1)/4.)
                    if self.appobj.dname != 'MET':
                        pvar2 = np.squeeze(np.sum(var2[:,ind[0:4],:],axis=1)/4.)
                        wwind = np.squeeze(np.sum(w[:,ind[0:4],:],axis=1)/4.)
                    hgt = np.squeeze(np.sum(height[:,ind[0:4],:],axis=1)/4.)
                    plevs = np.squeeze(np.sum(press[:,ind[0:4],:],axis=1)/4.)
                    hwind = np.squeeze(np.sum(self.u10[:,ind[0:4],:],axis=1)/4.)
        
                    #add the wind vectors in knots
#                    hwind*=1.943844492574
#                    wwind*=1.943844492574
        
                    xx = np.arange(0,self.u10.shape[2],int(self.dataSet.nx[self.currentGrid-1]/15.))
                    yy = np.arange(0,press.shape[0],int(self.dataSet.nz[self.currentGrid-1]/20.))
                    points = np.meshgrid(yy,xx)
                    hdim = self.dataSet.nx[self.currentGrid-1]
    
                if self.orientList[self.currentOrient] == 'yz':
                    diff = abs(self.dataSet.glons[self.dataSet.currentGrid-1][0,:]-self.ref_pt)
                    ind = np.argsort(diff)
                    horiz = np.tile(np.sum(self.dataSet.glats[self.dataSet.currentGrid-1][:,ind[0:4]],axis=1)/4.,(dims[0],1))
                    horiz2 =np.tile(np.sum(self.dataSet.glons[self.dataSet.currentGrid-1][:,ind[0:4]],axis=1)/4.,(dims[0],1))
                    pvar = np.squeeze(np.sum(var[:,:,ind[0:4]],axis=2)/4.)
                    if self.appobj.dname != 'MET':
                        pvar2 = np.squeeze(np.sum(var2[:,:,ind[0:4]],axis=2)/4.)
                        wwind = np.squeeze(np.sum(w[:,:,ind[0:4]],axis=2)/4.)
                    hgt = np.squeeze(np.sum(height[:,:,ind[0:4]],axis=2)/4.)
                    plevs = np.squeeze(np.sum(press[:,:,ind[0:4]],axis=2)/4.)
                    hwind = np.squeeze(np.sum(self.v10[:,:,ind[0:4]],axis=2)/4.)
        
                    #add the wind vectors
#                    hwind*=1.943844492574
#                    wwind*=1.943844492574
        
                    xx = np.arange(0,self.v10.shape[1],int(self.dataSet.nx[self.currentGrid-1]/15.))
                    yy = np.arange(0,press.shape[0],int(self.dataSet.nz[self.currentGrid-1]/20.))
                    points = np.meshgrid(yy,xx)    
                    hdim = self.dataSet.ny[self.currentGrid-1]    

                #Initialize colorbar range and increment
                if self.colormin is None:
                    #self.colormin = np.floor(pvar.min())
                    self.colormin = pvar.min()
                    if self.colormin != self.colormin:
                        self.colormin = 0.
                if self.colormax is None:
                    #self.colormax = np.ceil(pvar.max())
                    self.colormax = pvar.max()
                    if self.colormax != self.colormax:
                        self.colormax = 1.
                if self.ncontours is None:
                    self.ncontours = 11.
                if self.colormax == self.colormin:
                    if self.colormax == 0:
                        self.colormax = 1.
                    else:
                        self.colormax = self.colormax * 2.

                if self.colorbox is None:
                    self.controlColorBar()

                if self.vslicebox is None and self.currentPType == 1:
                    self.verticalSliceControl()

                if self.appobj.filltype == "contourf":
                    self.appobj.cs = self.appobj.axes1[self.pNum-1].contourf(horiz,
                          plevs, pvar, 
                          levels=np.linspace(self.colormin,self.colormax,
                          self.ncontours), extend=self.extend,
                          cmap=plt.cm.get_cmap(str(self.cmap)))
                elif self.appobj.filltype == "pcolormesh":
                    lvls = np.linspace(self.colormin,self.colormax,self.ncontours)
                    pvar = np.ma.masked_less_equal(pvar,self.colormin)
                    norm = matplotlib.colors.Normalize(vmin=np.amin(lvls),vmax=np.amax(lvls))
                    if (self.dataSet.variableList[self.currentVar] == 'REFL_10CM' or
                        np.ma.count_masked(pvar).any()):
                        self.appobj.axes1[self.pNum-1] = None
                        self.ColorBar = None
                        self.figure.clear()
                        if len(self.appobj.axes1) >= self.pNum:
                            self.appobj.axes1[self.pNum-1] = self.figure.add_subplot(111)
                        else:
                            self.appobj.axes1.append(self.figure.add_subplot(111))

                    self.appobj.cs = self.appobj.axes1[self.pNum-1].pcolormesh(
                                      horiz, plevs, pvar, norm=norm,
                                      cmap=plt.cm.get_cmap(str(self.cmap)))
                if self.ColorBar != None:
                    self.ColorBar.remove()
                self.Colorbar = self.figure.colorbar(self.appobj.cs, ax=self.appobj.axes1[self.pNum-1], 
                                                     orientation='horizontal', pad=0.1)
                self.Colorbar.ax.tick_params(labelsize=9)
                if (self.appobj.plotbarbs == True or self.appobj.plotvectors == True) and self.appobj.dname != 'MET':
                    if (self.dataSet.nz[self.currentGrid-1] < 60.):
                        interval = 2
                    else:
                        interval = 4
                    if (self.dataSet.ny[self.currentGrid-1] < 250. and self.dataSet.nx[self.currentGrid-1] < 250.):
                        interval2 = 10
                    else:
                        interval2 = 20
                    maxspeed = np.amax((hwind[::interval,::interval2]**2+wwind[::interval,::interval2]**2)**(1./2.))

                    if self.appobj.plotbarbs == True:
                        self.appobj.barbs = self.appobj.axes1[self.pNum-1].barbs(
                            horiz[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            plevs[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            hwind[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            wwind[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            length = 5, sizes={"emptybarb":0.0},barbcolor='k', flagcolor='k',linewidth=0.50, pivot='middle')

                        hnorm = hwind/np.sqrt(np.power(hwind,2) + np.power(wwind,2))
                        wnorm = wwind/np.sqrt(np.power(hwind,2) + np.power(wwind,2))

                        self.appobj.vectors2 = self.appobj.axes1[self.pNum-1].quiver(
                            horiz[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            plevs[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            hnorm[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            wnorm[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            pivot='middle', headwidth=0, headlength=0, headaxislength=0,scale=60,width=0.0015)

                    if self.appobj.plotvectors == True:
                        self.appobj.vectors = self.appobj.axes1[self.pNum-1].quiver(
                            horiz[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            plevs[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            hwind[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            wwind[int(interval/2.):self.dataSet.nz[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval2:hdim-1-int(interval2/2.):interval2],
                            color='k',units='inches',pivot='mid',scale=maxspeed*4)
                        self.appobj.vectorkey = self.appobj.axes1[self.pNum-1].quiverkey(
                                       self.appobj.vectors, 0.05, -0.08,
                                       int(maxspeed*.75), str(int(maxspeed*.75))+' $m s^{-1}$',labelpos='E')

                if self.appobj.plotcontour2 == True and self.appobj.dname != 'MET':
                    self.appobj.cs2 = self.appobj.axes1[self.pNum-1].contour(horiz, plevs,
                    pvar2, colors='k', linewidths=1.5)
                    if np.abs(pvar2).max() <= 1:
                        fmt = '%6.4f'
                    else:
                        fmt = '%d'
                    self.appobj.cs2label = self.appobj.cs2.clabel(inline=1, font=10, fmt=fmt)
                ax1 = self.appobj.axes1[self.pNum-1]
                ax1.set_ylabel("Pressure [hPa]")
                ax1.set_yscale('log')
                ax1.set_ylim(plevs.max()+1, self.minpress)
                print(plevs.max(),self.minpress)
                subs = [1,2,3,5,7,8.5]
                loc = matplotlib.ticker.LogLocator(base=10., subs=subs)
                ax1.yaxis.set_major_locator(loc)
                fmt = matplotlib.ticker.FormatStrFormatter("%g")
                ax1.yaxis.set_major_formatter(fmt)
                ax1.set_title(self.varTitle,fontsize = 10)
                ax2 = self.appobj.axes1[self.pNum-1].twinx()
                ax2.set_ylabel("Altitude [km]")
                if self.dataSet.dsetname == "MERRA":
                    if self.dataSet.grid[5] == 'P':
                        min_ind = 0
                    else:
                        min_ind = -1
                else:
                    min_ind = 0
                diff = abs(plevs[:,0] - self.minpress)
                ind = np.where(diff == diff.min())
                if len(ind[0]) > 1:
                    ind = ind[0]
                ax2.set_ylim(hgt.min(), hgt[ind[0],0])
                ax2.set_xlim(horiz[min_ind,:].min(),horiz[min_ind,:].max())
                ax2.fill_between(horiz[min_ind,:],0, hgt[min_ind,:], facecolor='peru')
                ax2.plot(horiz[min_ind,:],hgt[min_ind,:],color='black')
                fmt = matplotlib.ticker.FormatStrFormatter("%g")
                ax2.yaxis.set_major_formatter(fmt)
                axins = inset_axes(ax1,width="30%",height="30%",loc=1,borderpad=0)
                map2 = self.dataSet.map[self.currentGrid-1]
                map2.ax=axins
                map2.etopo()
                map2.drawcoastlines()
                map2.drawcountries()
                map2.drawstates()
                # draw parallels.
                plim = max(int((abs(self.dataSet.glats[self.currentGrid-1]).max()-abs(self.dataSet.glats[self.currentGrid-1]).min())/2.5),0.5)
                parallels = np.arange(-90.,90.,plim)
                map2.drawparallels(parallels,labels=[0,1,0,0],fontsize=6)
                # draw meridians
                mlim = max(int((abs(self.dataSet.glons[self.currentGrid-1]).max()-abs(self.dataSet.glons[self.currentGrid-1]).min())/2.5),0.5)
                meridians = np.arange(0.,360.,mlim)
                map2.drawmeridians(meridians,labels=[0,0,1,0],fontsize=6)
                if self.orientList[self.currentOrient] == 'xz':
                    xx, yy = map2(horiz[0,:], horiz2[0,:])
                if self.orientList[self.currentOrient] == 'yz':
                    xx, yy = map2(horiz2[0,:],horiz[0,:])
    
                map2.plot(xx,yy,color='red',linewidth=2)
                self.appobj.changeColor = False
                self.appobj.recallProjection = False

    def plotSkewT(self,i,j):

        #The WRF MET files have different pressure, temp and humidity variables      
        if self.appobj.dname == 'MET':
            #Get Variables
            p = np.squeeze(self.dataSet.readNCVariable('PRES'))/100. #Convert from pascal
            t = np.squeeze(self.dataSet.readNCVariable('TT'))-273.15 #Convert from K
            rh = np.squeeze(self.dataSet.readNCVariable('RH'))/100. #Convert from % 
            
            #Make sure RH values are valid (0-100%)
            rh[rh > 1.] = 1.
            rh[rh < 0.0] = 0.0
            #Calculate vapor pressure from temperature and RH
            qv = 6.11*np.exp((17.67*t)/(243.5+t))*rh
            self.skewt = SkewT.SkewTobj(temp = t[:,j,i], pressure = p[:,j,i], qv = qv[:,j,i], parcel = self.skewParcel)
                       

        else:
            #WRF
            #Get Variables 
            p = np.squeeze(self.dataSet.readNCVariable('P'))
            pb = np.squeeze(self.dataSet.readNCVariable('PB'))
            theta_prime = np.squeeze(self.dataSet.readNCVariable('T'))
            qv = np.squeeze(self.dataSet.readNCVariable('QVAPOR'))
        
            #Construct full fields        
            pressure = (p + pb) * 0.01
            temp = (theta_prime + 300.) * (pressure/1000.) ** (287.04/1004.)

            self.skewt = SkewT.SkewTobj(temp = temp[:,j,i], pressure = pressure[:,j,i], qv = qv[:,j,i], parcel = self.skewParcel)

        #Create Figure
        self.skewt.fig.canvas.set_window_title('Skew-T at i='+str(i)+' j='+str(j))
        self.skewt.ax.set_title(self.Plist[self.selectedParcel]+' Profile')        
        self.skewt.ax.set_xlabel('Temperature [$^o$C] \n'+self.dataSet.getTime())

    def plotVertPro(self,i,j):
        
        if self.vplot is not None:
            #Create new figure that will pop-up when called
            fig = plt.figure(figsize=(8, 7))
            fig.canvas.set_window_title('Vertical Profile at i='+str(i)+' j='+str(j))
            ax = fig.add_subplot(111)
            ax.set_xlabel(self.VertvarTitle)
            ax.set_ylabel("Pressure [hPa]")
            ax.set_yscale('log')
            ax2 = ax.twinx()
            ax2.set_ylabel("Altitude [km]")

            #Read in variables
            ph = np.squeeze(self.dataSet.readNCVariable('PH'))
            phb = np.squeeze(self.dataSet.readNCVariable('PHB'))
            p = np.squeeze(self.dataSet.readNCVariable('P'))
            pb = np.squeeze(self.dataSet.readNCVariable('PB'))
                
            #Create full fields at input location
            press = (p[:,j,i] + pb[:,j,i])/100.
            height = wrf.unstaggerZ((ph + phb)/9.81)/1000.
            height = height[:,j,i]
        
            #Get vertical profile of input field
            dims = self.vplot.shape
            if dims[0] != len(press):
                dvar = wrf.unstaggerZ(self.vplot)
                var = dvar[:,j,i] 
            else:
                var = self.vplot[:,j,i]       

            #Create plot
            ax.plot(var, press, 'k', linewidth=2.0)
            ax.xaxis.grid(True)
            #ax.set_ylim(press.max()*1.02,press.min())
            ax.set_ylim(1000., press.min())
            #mhgt = height.min()*np.log((press.max()*1.02)/press.max())
            mhgt = height.min()*np.log(1000./press.max())
            ax2.set_ylim(mhgt, height.max())
            subs = [0.5,1,2,3,4,5,6,7,8,9,10]
            loc = matplotlib.ticker.LogLocator(base=10., subs=subs)
            ax.yaxis.set_major_locator(loc)
            fmt = matplotlib.ticker.FormatStrFormatter("%g")
            ax.yaxis.set_major_formatter(fmt)
            ax.margins(0.02)
            ax.yaxis.grid(True)
        else:
           self.errorSelectVVar()

    def plotDifference(self):
        #Change cmap to blue-white-red'
        #icmap = plt.cm.get_cmap('bwr')
        icmap = matplotlib.cm.bwr
        #self.appobj.cmap = 'bwr'
        if self.appobj.recallProjection == True:
            if len(self.appobj.axes1) >= self.pNum:
                 self.appobj.axes1[self.pNum-1] = self.figure.add_subplot(111)
                 self.dataSet.resolution = self.appobj.resolution
                 self.dataSet.setProjection(self.currentGrid,axs=self.appobj.axes1[self.pNum-1])

            else:
                 self.appobj.axes1.append(self.figure.add_subplot(111))
                 self.dataSet.resolution = self.appobj.resolution
                 self.dataSet.setProjection(self.currentGrid,axs=self.appobj.axes1[self.pNum-1])
        xb, yb = self.dataSet.map[self.currentGrid-1](
                 self.dataSet.glons[self.currentGrid-1],
                 self.dataSet.glats[self.currentGrid-1])

        #Set background map
        alpha = 1
        if self.appobj.background is not None:
            pre = 'self.dataSet.map['
            gridnum = str(self.currentGrid-1)+']'
            end = self.appobj.background+'(ax=self.appobj.axes1['+str(self.pNum-1)+'])'
            exec(pre+gridnum+end)
            #Make variable transparent
            alpha=0.75
        if(self.currentVar != None or self.currentdVar != None):
            if(self.nz == 1):
                pltfld = self.var - self.diffvar
            else:
                pltfld = self.var[self.currentLevel] - self.diffvar[self.currentLevel]
            
            #Initialize colorbar range and increment
            if self.colormin is None:
                #self.colormin = np.floor(np.nanmin(pltfld))
                self.colormin = np.nanmin(pltfld)
                if self.colormin != self.colormin or self.colormin >= 0:
                    self.colormin = -1.
            if self.colormax is None:
                #self.colormax = np.ceil(np.nanmax(pltfld))
                self.colormax = np.nanmax(pltfld)
                if self.colormax != self.colormax or self.colormax <= 0:
                    self.colormax = 1.

            if self.ncontours is None:
                self.ncontours = 11.
            
            if self.colormax == self.colormin:
                if self.colormax == 0:
                    self.colormax = 1.
                else:
                    self.colormax = np.abs(self.colormax) * 2
            #Checks to make scale look better
            if self.colormin >= 0:
                self.colormin = -1.*self.colormax
            if self.colormax <=0:
                self.colormax = -1.*self.colormin

            if self.colorbox is None:
                self.controlColorBar()
            else:
                #This isnt the greatest way to handle the initial change but it works for now
                self.selectMin.setText(str(self.colormin))
                self.selectMax.setText(str(self.colormax))
                self.selectcontours.setText(str(self.ncontours))             
   
            if self.appobj.cs2 != None:
                for coll in self.appobj.cs2.collections:
                    coll.remove()
                for labels in self.appobj.cs2label:
                    labels.remove()
            if self.appobj.barbs != None:
                self.appobj.barbs.remove()
            if self.appobj.vectors != None:
                self.appobj.vectors.remove()
            if self.appobj.vectorkey != None:
                self.appobj.vectorkey.remove()
            if self.coasts != None and self.countries != None and self.states != None:
                self.coasts.remove()
                self.countries.remove()
                self.states.remove()
            
            #Colormap normalization
            midpoint = 1 - self.colormax/(self.colormax+abs(self.colormin))
            self.shiftedColorMap(cmap=icmap, midpoint=midpoint, name='shifted')
            self.cmap = self.newmap


            if self.appobj.filltype == "contourf":
                if self.appobj.cs is not None:
                    for coll in self.appobj.cs.collections:
                        coll.remove()
                self.appobj.cs = self.dataSet.map[self.currentGrid-1].contourf(
                                  self.dataSet.glons[self.currentGrid-1],
                                  self.dataSet.glats[self.currentGrid-1],
                                  pltfld,
                                  levels=np.linspace(self.colormin,self.colormax,self.ncontours),
                                  latlon=True,extend=self.extend,
                                  cmap=self.cmap,
                                  #cmap=plt.cm.get_cmap(str(self.appobj.cmap)),
                                  alpha=alpha, ax=self.appobj.axes1[self.pNum-1])
                divider = make_axes_locatable(self.appobj.axes1[self.pNum-1])
                cax = divider.append_axes("right", size="5%", pad=0.1)
                if self.ColorBar != None:
                    self.ColorBar.remove()
                self.ColorBar = self.figure.colorbar(self.appobj.cs,cax=cax, extend = self.extend)
            elif self.appobj.filltype == "pcolormesh":
                lvls = np.linspace(self.colormin,self.colormax,self.ncontours)
                pltfld = np.ma.masked_less_equal(pltfld,self.colormin)
                norm = matplotlib.colors.Normalize(vmin=np.amin(lvls),vmax=np.amax(lvls))
                self.appobj.cs = self.dataSet.map[self.currentGrid-1].pcolormesh(
                                  self.dataSet.glons[self.currentGrid-1],
                                  self.dataSet.glats[self.currentGrid-1],
                                  pltfld,
                                  norm=norm,
                                  latlon=True,cmap=self.cmap,
                                  alpha=alpha, ax=self.appobj.axes1[self.pNum-1])
                divider = make_axes_locatable(self.appobj.axes1[self.pNum-1])
                cax = divider.append_axes("right", size="5%", pad=0.1)
                if self.ColorBar != None:
                    self.ColorBar.remove()
                self.ColorBar = self.figure.colorbar(self.appobj.cs,cax=cax)
            self.ColorBar.ax.tick_params(labelsize=9)
            self.appobj.axes1[self.pNum-1].set_title(self.varTitle,fontsize = 10)
            self.displayValuesXYZ(pltfld)

            if self.appobj.plotbarbs == True or self.appobj.plotvectors == True:
                self.readField()

                if not self.derivedVar:
                    if self.nz != 1:
                        self.u10 = self.dataSet.u10[self.currentLevel]
                        self.v10 = self.dataSet.v10[self.currentLevel]
                    else:
                        self.u10 = self.dataSet.u10
                        self.v10 = self.dataSet.v10

                if (self.dataSet.dx[self.currentGrid-1] >= 15000):
                    if (self.dataSet.ny[self.currentGrid-1] < 250. and self.dataSet.nx[self.currentGrid-1] < 250.):
                        interval = 3
                    else:
                        interval = 5
                else:
                    if (self.dataSet.ny[self.currentGrid-1] < 250. and self.dataSet.nx[self.currentGrid-1] < 250.):
                        interval = 5
                    else:
                        interval = 10
                maxspeed = np.amax((self.u10[::interval,::interval]**2+self.v10[::interval,::interval]**2)**(1./2.))

                if self.appobj.plotbarbs == True:
                    self.appobj.barbs = self.appobj.axes1[self.pNum-1].barbs(
                        xb[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                        yb[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                        self.u10[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                        self.v10[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                        length = 5, barbcolor='k', flagcolor='k',linewidth=0.50)
                if self.appobj.plotvectors == True:
                        self.appobj.vectors = self.appobj.axes1[self.pNum-1].quiver(
                            xb[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                            yb[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                            self.u10[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                            self.v10[int(interval/2.):self.dataSet.ny[self.currentGrid-1]-1-int(interval/2.):int(interval/2.),interval:self.dataSet.nx[self.currentGrid-1]-1-int(interval/2.):interval],
                            color='k',units='inches',pivot='mid',scale=maxspeed*4)
                        self.appobj.vectorkey = self.appobj.axes1[self.pNum-1].quiverkey(self.appobj.vectors, 0.05, -0.08,
                                                int(maxspeed*.75), str(int(maxspeed*.75))+' $m s^{-1}$',labelpos='E')

            if self.appobj.plotcontour2 == True:
                self.readField()
                if not self.derivedVar:
                    self.var2 = pltfld
                self.appobj.cs2 = self.appobj.axes1[self.pNum-1].contour(
                                xb, yb, self.var2,
                                colors='k',linewidths=1.5)
                if np.abs(self.var2).max() <= 1:
                    fmt = '%8.6f'
                else:
                    fmt = '%d'
                self.appobj.cs2label = self.appobj.cs2.clabel(inline=1,
                                       font=10, fmt=fmt)

        self.coasts = self.dataSet.map[self.currentGrid-1].drawcoastlines(ax=self.appobj.axes1[self.pNum-1])
        self.countries = self.dataSet.map[self.currentGrid-1].drawcountries(ax=self.appobj.axes1[self.pNum-1])
        self.states = self.dataSet.map[self.currentGrid-1].drawstates(ax=self.appobj.axes1[self.pNum-1])
        # draw parallels
        if self.dataSet.projectionType == "robin" or self.dataSet.projectionType == "geos" or \
           self.dataSet.projectionType == "ortho" or self.dataSet.projectionType == "aeqd":
            parallels = np.arange(-90.,90.,15)
            meridians = np.arange(0.,360., 30)
        else:
            plim = max(int((np.nanmax(np.array(np.abs(self.dataSet.glats[self.currentGrid-1])))
                           -np.nanmin(np.array(np.abs(self.dataSet.glats[self.currentGrid-1]))))/10.),1)
            parallels = np.arange(-90.,90.,plim)
            mlim = max(int((np.nanmax(np.array(np.abs(self.dataSet.glons[self.currentGrid-1])))
                           -np.nanmin(np.array(np.abs(self.dataSet.glons[self.currentGrid-1]))))/10.),1)
            meridians = np.arange(0.,360.,mlim)

        #Can't put labels on Geostationary, Orthographic or Azimuthal Equidistant basemaps
        if self.dataSet.projectionType == "ortho" or self.dataSet.projectionType == "aeqd" or self.dataSet.projectionType == "geos":
            self.dataSet.map[self.currentGrid-1].drawparallels(parallels,ax=self.appobj.axes1[self.pNum-1])
            # draw meridians
            self.dataSet.map[self.currentGrid-1].drawmeridians(meridians,ax=self.appobj.axes1[self.pNum-1])
        else:
            self.dataSet.map[self.currentGrid-1].drawparallels(parallels,labels=[1,0,0,0],fontsize=6,ax=self.appobj.axes1[self.pNum-1])
            # draw meridians
            self.dataSet.map[self.currentGrid-1].drawmeridians(meridians,labels=[0,0,0,1],fontsize=6,ax=self.appobj.axes1[self.pNum-1])

        self.appobj.recallProjection = False
        line, = self.appobj.axes1[self.pNum-1].plot([0], [0])  # empty line
        self.dataSelector = DataSelector(line)

    def plotTimeSeries(self):
        
        #Create new figure that will pop-up when called
        fig = plt.figure(figsize=(12, 8))
        fig.canvas.set_window_title('Time Series Plot')
        ax = fig.add_subplot(111)
        ax.set_ylabel(self.varTitle.split("\n")[0])
        ax.set_xlabel('Date-Time')        

        #Create a point for each time
        ind = matplotlib.dates.date2num(self.time_series)
        #Plot
        if len(self.var_series.shape) > 1:
            #3D Variable
            ax.plot_date(ind,self.var_series[:,self.currentLevel],'bo-')
        else:
            #2D Variable
            ax.plot_date(ind,self.var_series,'bo-')
        
        #Create date-time tick marks
        plt.xticks(fontsize=10)
        #Make conditional spacing options
        time_diff = self.time_series[-1] - self.time_series[0]
        tdiff_seconds = time_diff.days*24.*3600. + time_diff.seconds
        if (tdiff_seconds >= 10.*24.*3600.):
            ax.xaxis.set_major_locator(matplotlib.dates.DayLocator(interval=2))
            ax.xaxis.set_minor_locator(matplotlib.dates.HourLocator(interval=12))
            ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d-%Y %H:%M'))
        elif (tdiff_seconds >= 5.*24.*3600.):
            ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=12))
            ax.xaxis.set_minor_locator(matplotlib.dates.HourLocator(interval=6))
            ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d-%Y %H:%M'))
        elif (tdiff_seconds >= 3.*24.*3600.):
            ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=6))
            ax.xaxis.set_minor_locator(matplotlib.dates.HourLocator(interval=3))
            ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d-%Y %H:%M'))
        elif (tdiff_seconds >= 0.5*24.*3600.):
            ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
            ax.xaxis.set_minor_locator(matplotlib.dates.HourLocator(interval=1))
            ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d-%Y %H:%M'))
        else:
            ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=15))
            ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(interval=5))
            ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d-%Y %H:%M'))


        plt.tight_layout()
        plt.margins(0.01)
        plt.grid(True)
        fig.autofmt_xdate()

    def createTimeSeries(self,i,j):
        #t0 = time.clock()
        #t1 = time.time()
        QApplication.setOverrideCursor(Qt.WaitCursor) 
        #Create list to hold values at selected point and time
        var_series = []
        self.time_series = []
        append_var = var_series.append
        append_time = self.time_series.append
        numfiles = self.dataSet.getNumFiles()
        ntimes = self.dataSet.ntimes
        input_time_index = self.dataSet.currentTimeIndex
        #If the number of dimensions is less than 3 then the selected variable is 2D
        if len(self.var.shape) < 3:
            for ii in range(numfiles*ntimes):
                timeIndex = self.dataSet.setTimeIndex(ii)
                self.dataSet.getTime()
                append_time(self.dataSet.timeObj)#.strftime("%Y-%m-%d %H:%M:%S"))
                self.readField()
                append_var(self.var[j,i])
        else:
            for ii in range(numfiles*ntimes):
                timeIndex = self.dataSet.setTimeIndex(ii)
                self.dataSet.getTime()
                append_time(self.dataSet.timeObj)
                self.readField()
                append_var(self.var[:,j,i])        
        self.dataSet.setTimeIndex(input_time_index)
        self.var_series = np.array(var_series) 
        QApplication.restoreOverrideCursor()
        #print(time.clock() - t0, "Seconds process time")
        #print(time.time() - t1, "Seconds wall time")

    def displayValuesXYZ(self,var):
        self.appobj.axes1[self.pNum-1].imshow(var, interpolation='nearest')
        #numrows, numcols = var.shape
        def format_coord(x, y):
            if self.dataSet.dsetname == "GOES Class":
                clon,clat = self.dataSet.map[self.currentGrid-1](self.dataSet.lon0[self.currentGrid-1],self.dataSet.lat0[self.currentGrid-1])
                #Get accurate lat lon points with nan lat longs in the dataset
                lon, lat  = self.dataSet.map[self.currentGrid-1](self.dataSet.maplon[self.currentGrid-1], self.dataSet.maplat[self.currentGrid-1])
                tmp = np.argsort(np.abs(lat[:,int(self.dataSet.nx[self.currentGrid -1]/2)] - clat))[0]
                tmp2 = np.argsort(np.abs(lon[tmp,:] - clon))[0]
                row = np.argsort(np.abs(lat[:,tmp2] - y))[0]
                col  = np.argsort(np.abs(lon[row,:] -x))[0]
                x1 = self.dataSet.maplon[self.currentGrid-1][row,col]
                y1 = self.dataSet.maplat[self.currentGrid-1][row,col]
                if x1 >= 180.:
                    x1 = -360. + x1
                #Now for the data that is on original lon/lat grid
                lon, lat  = self.dataSet.map[self.currentGrid-1](self.dataSet.glons[self.currentGrid-1], self.dataSet.glats[self.currentGrid-1])
                tmp = np.argsort(np.abs(lat[:,int(self.dataSet.nx[self.currentGrid -1]/2)] - clat))[0]
                tmp2 = np.argsort(np.abs(lon[tmp,:] - clon))[0]
                row = np.argsort(np.abs(lat[:,tmp2] - y))[0]
                col  = np.argsort(np.abs(lon[row,:] -x))[0]
                z = var[row,col]
            else:
                lon, lat  = self.dataSet.map[self.currentGrid-1](self.dataSet.glons[self.currentGrid-1], self.dataSet.glats[self.currentGrid-1])
                clon,clat = self.dataSet.map[self.currentGrid-1](self.dataSet.lon0[self.currentGrid-1],self.dataSet.lat0[self.currentGrid-1])
                tmp = np.argsort(np.abs(lat[:,int(self.dataSet.nx[self.currentGrid -1]/2)] - clat))[0]
                tmp2 = np.argsort(np.abs(lon[tmp,:] - clon))[0]
                row = np.argsort(np.abs(lat[:,tmp2] - y))[0]
                col  = np.argsort(np.abs(lon[row,:] -x))[0]
                x1 = self.dataSet.glons[self.currentGrid-1][row,col]
                y1 = self.dataSet.glats[self.currentGrid-1][row,col]
                if x1 >= 180.:
                    x1 = -360. + x1
                z = var[row,col]

            #if self.currentPType == 2 or self.currentPType == 3 or self.currentPType == 4:
            #   return 'i=%6.2f, j=%6.2f, value=%10.6f'%(col, row, z)
            #else:
            return 'lon=%6.2f, lat=%6.2f, i=%6.2f, j=%6.2f, value=%10.6f'%(x1, y1, col, row, z)
        if self.dataSet.dsetname != 'NEXRAD Radar':
            self.appobj.axes1[self.pNum-1].format_coord = format_coord

    def getTable(self,col,row):
        self.table = QTableView()
        self.table.resizeColumnsToContents()
        self.table.resizeRowsToContents()
        self.table.setWindowTitle("Equations of Motion")
        #self.table.setStyle(QStyleFactory.create("Windows"))
        self.model = QStandardItemModel()
        row1 = ['X-Comp', 'Y-Comp', 'Z-Comp']
        col1 = ['PGF', 'Cor','Grav']
        QStandardItemModel.setHorizontalHeaderLabels(self.model,col1)
        QStandardItemModel.setVerticalHeaderLabels(self.model,row1)
        if(self.nz == 1):
            var = self.var
        else:
            var = self.var[self.currentLevel]

        pgfx, pgfy, pgfz = self.dataSet.calcPGF(col=col,row=row,level=self.currentLevel)
        corx, cory, corz = self.dataSet.calcCOR(col=col,row=row,level=self.currentLevel)
        self.model.setItem(0,0, QStandardItem("{:.3E}".format(pgfx)))
        self.model.setItem(1,0, QStandardItem("{:.3E}".format(pgfy)))
        self.model.setItem(2,0, QStandardItem("{:.3E}".format(pgfz)))
        self.model.setItem(0,1, QStandardItem("{:.3E}".format(corx)))
        self.model.setItem(1,1, QStandardItem("{:.3E}".format(cory)))
        self.model.setItem(2,1, QStandardItem("{:.3E}".format(corz)))
        self.model.setItem(0,2, QStandardItem("N/A"))
        self.model.setItem(1,2, QStandardItem("N/A"))
        self.model.setItem(2,2, QStandardItem("9.81"))

        self.table.setModel(self.model)
        self.table.show()

    def selectionChangePlot(self,i):
        if self.cid is not None:
            self.cid = None
            self.row = None
            self.col = None
        if (self.currentPType == 5):
            self.cmap = self.appobj.cmap
        self.currentPType = i
		            
        #Allow changing of plot type before plotting
        if (self.currentVar != None or self.currentdVar != None):
            if self.currentPType == 1 and len(self.var.shape) < 3:
                self.error3DVar()
                self.ColorBar = None
            else:
                if self.vslicebox != None:
                    self.vslicebox.setParent(None)
                    self.vslicebox = None
                self.appobj.cs = None
                self.appobj.cs2 = None
                self.appobj.cs2label = None
                self.appobj.barbs = None
                self.appobj.vectors = None
                self.appobj.vectorkey = None
                self.appobj.domain_average = None
                self.appobj.recallProjection = True
                self.coasts = None
                self.countries = None
                self.states = None
                self.ColorBar = None
                if (self.currentVar != None):
                    if (self.dataSet.variableList[self.currentVar] != "REFL_10CM"):
                        self.colormax = None
                        self.colormin = None
                        self.ncontours = None
                if (self.currentdVar != None):
                    if (self.dataSet.dvarlist[self.currentdVar] != "refl"):
                        self.colormax = None
                        self.colormin = None
                        self.ncontours = None
            
                self.appobj.axes1[self.pNum-1] = None
                self.figure.clear()
 
                #Destroy vertical plot control widget if plot type is changed
                if self.vertControl != None:
                    self.vertControl.setParent(None)
                    self.vertControl = None

                #Destroy skew-T parcel control widget if plot type is changed
                if self.pControl != None:
                    self.pControl.setParent(None)
                    self.pControl = None

                #Destroy skew-T parcel control widget if plot type is changed
                if self.diffControl != None:
                    self.diffControl.setParent(None)
                    self.diffControl = None        

                #Skew-T
                if self.currentPType == 2:
                    #Change to mixed-layer CAPE plot
                    #i = np.where(np.array(self.dataSet.dvarlist) == 'CAPE_ML')
                    #if self.derivedVar == False:
                    #    self.selectionChangedVar(i[0][0])
                    #    self.replot2d = True
                    #else:
                    #    if self.currentdVar != i[0][0]:
                    #       self.selectionChangedVar(i[0][0])
                    #       self.replot2d = True
                    #    else:
                    #       self.replot2d = False

                    #Create tab to control parcel type
                    self.pControl = QGroupBox()
                    parTitle = 'Parcel Selection'
                    self.pControlLayout = QVBoxLayout(self.pControl)

                    selectVarWidget = QWidget()
                    selectVarWidgetLayout = QHBoxLayout()
                    selectVarWidget.setLayout(selectVarWidgetLayout)

                    selectVarLabel = QLabel()
                    selectVarLabel.setText('Parcel:')

                    selectPVar = QComboBox()
                    selectPVar.setStyleSheet(Layout.QComboBox())
                    self.Plist = ['Surface Based', 'Mixed Layer', 'Most Unstable']
                    selectPVar.addItems(self.Plist)
                    selectPVar.setSizeAdjustPolicy(QComboBox.AdjustToContents)
                    selectPVar.currentIndexChanged.connect(self.selectionChangeParcel)

                    selectVarWidgetLayout.addWidget(selectVarLabel)
                    selectVarWidgetLayout.addWidget(selectPVar)
                    self.pControlLayout.addWidget(selectVarWidget)
                    self.optionTabs.addTab(self.pControl,parTitle)
                    self.optionTabs.setCurrentIndex(self.optionTabs.count()-1)
                    self.tabbingLayout.addWidget(self.optionTabs)
                    self.pltFxn(self.pNum)

                #Vertical Profiles
                if self.currentPType == 3:
                    #Create tab to control vertical plot variable
                    self.vertControl = QGroupBox()
                    vertTitle = 'Vertical Profile Control'
                    self.vertControlLayout = QVBoxLayout(self.vertControl)
            
                    selectVarWidget = QWidget()
                    selectVarWidgetLayout = QHBoxLayout()
                    selectVarWidget.setLayout(selectVarWidgetLayout)

                    selectVarLabel = QLabel()
                    selectVarLabel.setText('Vertical Variable:')

                    selectVVar = QComboBox()
                    selectVVar.setStyleSheet(Layout.QComboBox())
                    self.Vvarlist = self.dataSet.threeDVars
                    selectVVar.addItems(self.Vvarlist)
                    selectVVar.setSizeAdjustPolicy(QComboBox.AdjustToContents)
                    selectVVar.currentIndexChanged.connect(self.selectionChangeVerticalVar)       
                    selectVarWidgetLayout.addWidget(selectVarLabel)
                    selectVarWidgetLayout.addWidget(selectVVar)
                    self.vertControlLayout.addWidget(selectVarWidget)
                    self.optionTabs.addTab(self.vertControl,vertTitle)
                    self.optionTabs.setCurrentIndex(self.optionTabs.count()-1)
                    self.tabbingLayout.addWidget(self.optionTabs)                  

                #Time Series
                if self.currentPType == 4:
                    pass     

                #Difference Plot
                if self.currentPType == 5:
                    #Reset colorbar settings
                    self.colormax = None
                    self.colormin = None
                    self.ncontours = None
                    self.appobj.cmap = self.cmap
                    #Create tab to control the difference plot options
                    self.diffControl = QGroupBox()
                    diffTitle = 'Difference Plot Control'
                    self.diffControlLayout = QVBoxLayout(self.diffControl)

                    selectdataWidget = QWidget()
                    selectdataWidgetLayout = QHBoxLayout(selectdataWidget)

                    selectdataLabel = QLabel()
                    selectdataLabel.setText('Comparison dataset:')

                    count = 0
                    self.dsetlist = []
                    for i in self.dSet:
                        if count == 0:
                            count += 1
                        else:
                            dsetname = os.path.basename(str(self.dSet[count].path))
                            self.dsetlist.append(dsetname+' [Dataset '+str(count)+']')
                            count += 1

                    self.selectData = QComboBox()
                    self.selectData.setStyleSheet(Layout.QComboBox())
                    self.selectData.addItems(self.dsetlist)
                    self.selectData.setCurrentIndex(len(self.dsetlist)-1)
                    self.selectData.setSizeAdjustPolicy(QComboBox.AdjustToContents)
                    self.selectData.activated.connect(self.selectionChangeDiffData)

                    selectdataWidgetLayout.addWidget(selectdataLabel)
                    selectdataWidgetLayout.addWidget(self.selectData)
                    self.diffControlLayout.addWidget(selectdataWidget)
                    self.optionTabs.addTab(self.diffControl,diffTitle)
                    self.optionTabs.setCurrentIndex(self.optionTabs.count()-1)
                    self.tabbingLayout.addWidget(self.optionTabs)

                    #Popup to select dataset for difference
                    self.dControl = QDialog(self.appobj)
                    self.dControl.resize(self.dControl.minimumSizeHint())
                    self.dControl.setWindowTitle("Difference Dataset Selection")
                    dControlLayout = QVBoxLayout(self.dControl)
        
                    self.buttonGroup = QButtonGroup()
                    
                    #Create a button for each dataset
                    count = 0
                    self.buttonname = []
                    for i in self.dsetlist:
                        self.button_name = QRadioButton("{}".format(i))
                        self.button_name.setObjectName("radiobtn_{}".format(i))
                        self.buttonname.append(self.button_name)
                        dControlLayout.addWidget(self.button_name)
                        self.buttonGroup.addButton(self.button_name,count)
                        self.button_name.clicked.connect(lambda:self.determineRadioSelection())
                        count += 1    
                    self.dControl.show() 
		    
                self.pltFxn(self.pNum) 
        else:
            self.errorPlotChange()
            self.currentPType = 0
            self.dataSet.selectPlotType.setCurrentIndex(self.currentPType)

        
    def determineRadioSelection(self):
        for i in range(len(self.dsetlist)):
            if self.buttonname[i].isChecked(): 
                self.dControl.close()
                self.selectData.setCurrentIndex(i)
                self.selectionChangeDiffData(i)                

    def selectionChangeDiffData(self,i):
        self.vartitle_ext = '('+self.dataSet.selectDset.itemText(self.dataSet.selectDset.currentIndex())+' - '+self.dsetlist[i]+')'
        #Add 1 to i because the dSet variable starts with a None
        i += 1
        self.diffdata = self.dSet[i]
        self.diffdata.setGrid(self.currentGrid, update = False)
        self.diffdata.setTimeIndex(self.currentTime)
        self.readField()
        self.readDiffField()
        self.pltFxn(self.pNum)

    def selectionChangeVerticalVar(self,i):
        self.selectedVvar = i
        self.vplot = self.dataSet.readNCVariable(self.Vvarlist[i],
                     barbs=self.appobj.plotbarbs, contour2=self.appobj.plotcontour2)
        self.VertvarTitle = self.dataSet.description +' ('+self.dataSet.units
        self.VertvarTitle = self.VertvarTitle +') \n'+self.dataSet.getTime()

    def selectionChangeParcel(self,i):
        self.selectedParcel = i
        if self.selectedParcel == 0:
            self.skewParcel = 'SB'
        elif self.selectedParcel == 1:
            self.skewParcel = 'ML'
        else:
            self.skewParcel = 'MU'

    def selectionChangeVar(self,i):
        self.currentVar = i
        self.derivedVar = False
        self.readField()
        if self.currentPType == 1 and len(self.var.shape) < 3:
            self.error3DVar()
        else:
            if self.colorlock == False:
                #Reset colorbar settings
                self.colormax = None
                self.colormin = None
                self.ncontours = None
            if self.colorbox is not None:
                self.colorbox.setParent(None)
                self.colorbox = None
            if self.vslicebox is not None:
                self.vslicebox.setParent(None)
                self.vslicebox = None
            #self.currentVar = i
            #self.derivedVar = False
            #self.readField()
            if (self.dataSet.variableList[self.currentVar] == 'REFL_10CM'):
                self.extend = 'max'
                self.colormin = 0
                self.colormax = 80
                if self.cmap != 'pyart_NWSRef':
                    self.appobj.cmap = self.cmap
                self.cmap = 'pyart_NWSRef'
            else:
                self.extend = 'both'
                self.cmap = self.appobj.cmap
            ndim   = len(self.var.shape)
            if self.currentPType == 3 or self.currentPType == 2:
                self.replot2d = True
                self.col = None
                self.row = None
            if self.currentPType == 4:
                self.replot2d = True
                self.col = None
                self.row = None
                self.newvar = True
            #Difference plot - read in difference data
            if self.currentPType == 5:
                self.readDiffField()
            if(ndim > 1) :
                if(len(self.var.shape) == 2):
                    self.nz=1
                else:
                    self.nz=self.var.shape[0]
                self.selectLevel.clear()
                self.updateListView = False
                for j in range(0,self.nz):
                    self.selectLevel.addItem(str(j+1))
                self.updateListView = True
                self.pltFxn(self.pNum) 

    def selectionChangedVar(self,i):
        self.currentdVar = i
        self.derivedVar = True
        self.readField()
        if self.currentPType == 1 and len(self.var.shape) < 3:
            self.error3DVar()
        else:
            if self.colorlock == False:
                #Reset colorbar settings
                self.colormax = None
                self.colormin = None
                self.ncontours = None
            if self.colorbox is not None:
                self.colorbox.setParent(None)
                self.colorbox = None
            if (self.dataSet.dvarlist[self.currentdVar] == 'refl'):
                self.extend = 'max'
                self.colormin = 0
                self.colormax = 80
                if self.cmap != 'pyart_NWSRef':
                    self.appobj.cmap = self.cmap
                self.cmap = 'pyart_NWSRef'
            else:
                self.extend = 'both'
                self.cmap = self.appobj.cmap

            ndim = len(self.var.shape)
            if self.currentPType == 3 or self.currentPType == 2:
                self.replot2d = True
                self.row = None
                self.col = None
            if self.currentPType == 4:
                self.replot2d = True
                self.col = None
                self.row = None
                self.newvar = True
            #Difference plot - read in difference data
            if self.currentPType == 5:
                self.readDiffField()
            if(ndim > 1) :
                if(len(self.var.shape) == 2):
                    self.nz=1
                else:
                    self.nz=self.var.shape[0]
                self.selectLevel.clear()
                self.updateListView = False
                for j in range(0,self.nz):
                    self.selectLevel.addItem(str(j+1))
                self.updateListView = True
                self.pltFxn(self.pNum)

    def selectionChangeGrid(self,i):
        self.ColorBar = None
        self.appobj.cs = None
        self.appobj.cs2 = None
        self.appobj.barbs = None
        self.appobj.vectors = None
        self.appobj.vectorkey = None
        self.appobj.cs2label = None
        self.appobj.domain_average=None
        self.coasts = None
        self.countries = None
        self.states = None
        self.appobj.recallProjection = True
        self.currentGrid = i+1
        self.currentTime = 0
        if self.colorlock == False:
            #Reset colorbar settings
            self.colormax = None
            self.colormin = None
            self.ncontours = None       
        varname = None
        if self.derivedVar == True:
            if self.currentdVar != None:
                varname = self.dataSet.dvarlist[self.currentdVar]
        else:
            if self.currentVar != None:
                varname = self.dataSet.variableList[self.currentVar]
        self.dataSet.setTimeIndex(self.currentTime)
        self.dataSet.setGrid(self.currentGrid)
        self.selectTime.clear()
        self.selectTime.addItems(self.dataSet.timeList[self.dataSet.currentGrid])
        self.selectVar.clear()
        self.selectVar.addItems(self.dataSet.variableList)
        if varname != None and self.derivedVar == True:
            self.currentdVar = np.where(np.array(self.dataSet.dvarlist) == varname)[0]
            if (len(self.currentdVar) > 0):
                self.currentdVar = self.currentdVar[0]
                self.dataSet.selectdVar.setCurrentIndex(self.currentdVar)
            else:
                self.varError(varname=varname)
                self.currentVar = None
                self.currentdVar = None
        if varname != None and self.derivedVar == False:
            self.currentVar = np.where(np.array(self.dataSet.variableList) == varname)[0]
            if (len(self.currentVar) > 0):
                self.currentVar = self.currentVar[0]
                self.selectVar.setCurrentIndex(self.currentVar)
            else:
                self.varError(varname=varname) 
                self.currentVar = None
                self.currentdVar = None
        self.readField()
        #Difference plot - read in difference data
        if self.currentPType == 5:
            self.diffdata.setTimeIndex(self.currentTime)
            self.diffdata.setGrid(self.currentGrid)
            self.readDiffField()
        self.appobj.axes1[self.pNum-1] = None
        self.figure.clear()
        if self.currentPType == 1 and len(self.var.shape) < 3:
            self.error3DVar()
        else:
            self.pltFxn(self.pNum)
        
    def selectionChangeLevel(self,i):
        if(i != -1 and self.updateListView == True):
            if self.colorlock == False:
                #Reset colorbar settings
                self.colormax = None
                self.colormin = None
                self.ncontours = None

            self.currentLevel = i
            self.readField()
            #Difference plot - read in difference data
            if self.currentPType == 5:
                self.readDiffField()
            self.pltFxn(self.pNum)

    def selectionChangeTime(self,i):
        self.currentTime = i
        self.dataSet.setTimeIndex(self.currentTime)
        if self.colorlock == False:
            #Reset colorbar settings
            self.colormax = None
            self.colormin = None
            self.ncontours = None
        if self.currentPType == 2:
            self.replot2d = True
            self.col = None
            self.row = None
        if self.currentPType == 3:
            self.replot2d = True
            self.col = None
            self.row = None
            self.selectionChangeVerticalVar(self.selectedVvar)
        self.readField()
        #Difference plot - read in difference data
        
        if self.currentPType == 5:
            self.diffdata.setTimeIndex(self.currentTime)
            self.readDiffField()
        self.pltFxn(self.pNum)

    def selectionChangeDset(self,i):
        self.currentDset = i+1
        if self.colorbox is not None:
            self.colorbox.setParent(None)
            self.colorbox = None
        if self.vslicebox is not None:
            self.vslicebox.setParent(None)
            self.vslicebox = None
        self.appobj.cbar.deleteLater()
        self.appobj.cbar = None
        self.appobj.cs = None
        self.appobj.cs2 = None
        self.appobj.barbs = None
        self.appobj.vectors = None
        self.appobj.vectorkey = None
        self.appobj.cs2label = None
        self.appobj.domain_average = None
        self.ColorBar = None
        self.derivedVar = False
        self.appobj.recallProjection = True
        #Plus 1 is needed because None is technically the first index
        self.dataSet = self.dSet[self.currentDset]
        self.getControlBar()
        self.currentGrid = 1
        self.currentTime = 0
        self.dataSet.setTimeIndex(self.currentTime)
        self.dataSet.setGrid(self.currentGrid, update=None)
        self.currentVar = None
        self.currentdVar = None
        self.selectVar.clear()
        self.selectVar.addItems(self.dataSet.variableList)
        #self.appobj.axes1.remove(self.appobj.axes1[self.pNum-1])
        self.appobj.axes1[self.pNum-1] = None
        self.figure.clear()
        self.pltFxn(self.pNum)

    def nxtButtonAction(self):
        self.currentTime+=1
        if self.currentTime == self.dataSet.ntimes*self.dataSet.getNumFiles():
            self.currentTime = 0
        self.selectTime.setCurrentIndex(self.currentTime)
        self.dataSet.setTimeIndex(self.currentTime)
        if self.currentPType == 2:
            self.replot2d = True
            self.col = None
            self.row = None
        if self.currentPType == 3:
            self.replot2d = True
            self.col = None
            self.row = None
            self.selectionChangeVerticalVar(self.selectedVvar)        
        self.readField()
        if self.colorlock == False:
            #Reset colorbar settings
            self.colormax = None
            self.colormin = None
            self.ncontours = None
        #Difference plot - read in difference data
        if self.currentPType == 5:
            self.diffdata.setTimeIndex(self.currentTime)
            self.readDiffField()
        self.pltFxn(self.pNum)
	 
    def prevButtonAction(self):
        self.currentTime-=1
        if self.currentTime == -1:
            self.currentTime = self.dataSet.getNumFiles()*self.dataSet.ntimes-1
        self.selectTime.setCurrentIndex(self.currentTime)
        self.dataSet.setTimeIndex(self.currentTime)
        if self.currentPType == 2:
            self.replot2d = True
            self.col = None
            self.row = None
        if self.currentPType == 3:
            self.replot2d = True
            self.col = None
            self.row = None
            self.selectionChangeVerticalVar(self.selectedVvar)          
        self.readField()
        if self.colorlock == False:
            #Reset colorbar settings
            self.colormax = None
            self.colormin = None
            self.ncontours = None
        #Difference plot - read in difference data
        if self.currentPType == 5:
            self.diffdata.setTimeIndex(self.currentTime)
            self.readDiffField()
        self.pltFxn(self.pNum)

    def selectionChangeOrient(self,i):
        self.currentOrient = i
        if self.orientList[self.currentOrient] == 'xz':
            self.refboxLabel.setText('Lat:')
        if self.orientList[self.currentOrient] == 'yz':
            self.refboxLabel.setText('Lon:')

    def enterPress(self):
        self.ref_pt = np.float(self.refbox.text())
        self.minpress = np.float(self.pbox.text())
        self.figure.clear()
        self.pltFxn(self.pNum)

    def enterCValue(self):
        self.colormin = np.float(self.selectMin.text())
        self.colormax = np.float(self.selectMax.text())
        self.ncontours = np.float(self.selectcontours.text())
        self.pltFxn(self.pNum)
    
    def error3DVar(self):
        if self.derivedVar == True:
            varname = self.dataSet.dvarlist[self.currentdVar]
        else:
            varname = self.dataSet.variableList[self.currentVar]
        msg = QMessageBox(self.appobj)
        msg.setIcon(QMessageBox.Information)
        msg.setText(varname+' is not a 3D variable. The selected plot type requires a 3D variable.')
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()  

    def errorPlotChange(self):
        msg = QMessageBox(self.appobj)
        msg.setIcon(QMessageBox.Information)
        msg.setText('Please select a variable before changing the plot type.')
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()      

    def varError(self,varname):
        msg = QMessageBox(self.appobj)
        msg.setIcon(QMessageBox.Information)
        msg.setText(varname+' is not available in the new grid\n Select a new variable.')
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()
    
    def errorSelectVVar(self):
        msg = QMessageBox(self.appobj)
        msg.setIcon(QMessageBox.Information)
        msg.setText('Please select a variable from the pop-up Vertical Profile Control')
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def shiftedColorMap(self, cmap=None, start=0.0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
        #From Paul H on stackoverflow
        '''
        Function to offset the "center" of a colormap. Useful for
        data with a negative min and positive max and you want the
        middle of the colormap's dynamic range to be at zero

        Input
        -----
        cmap : The matplotlib colormap to be altered
        start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
        midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
        stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
        '''
        cdict = {
            'red': [],
            'green': [],
            'blue': [],
            'alpha': []
        }

        # regular index to compute the colors
        reg_index = np.linspace(start, stop, 257)

        # shifted index to match the data
        shift_index = np.hstack([
            np.linspace(0.0, midpoint, 128, endpoint=False), 
            np.linspace(midpoint, 1.0, 129, endpoint=True)
        ])

        for ri, si in zip(reg_index, shift_index):
            r, g, b, a = cmap(ri)

            cdict['red'].append((si, r, r))
            cdict['green'].append((si, g, g))
            cdict['blue'].append((si, b, b))
            cdict['alpha'].append((si, a, a))

        self.newmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
        plt.register_cmap(cmap=self.newmap)
        #self.newmap = cdict        

    #return newcmap

###############################################################################
####                           End PlotSlab() Object                       ####
###############################################################################


###############################################################################
####                         Begin AppForm() Object                        ####
###############################################################################
        
class AppForm(QMainWindow):
    
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        #self.data = self.get_data2()
        self.plotCount = 0
        self.plotTab = []
        self.create_main_frame()
        self.setWindowTitle('PyGEOMET')
        
        #Main Menu
        extractAction = QAction("&Exit",self)
        extractAction.setShortcut("Ctrl+Q")
        extractAction.setStatusTip("Leave the App")
        extractAction.triggered.connect(self.close_program)        
        
        #openDir = QAction("&Select Directory", self)
        #openDir.setShortcut("Ctrl+D")
        #openDir.setStatusTip('Open Directory')
        #openDir.triggered.connect(self.wrfOpen)
        
        openEditor = QAction("&Editor", self)
        openEditor.setShortcut("Ctrl+E")
        openEditor.setStatusTip('Open Editor')
        
        #Create Dataset menu bar###############################
        #WRF
        openWRF = QAction("&WRF-ARW", self)
        openWRF.setShortcut("Ctrl+W")
        openWRF.setStatusTip('Open WRF Directory')
        openWRF.triggered.connect(self.wrfOpen)
        
        #GOES Class
        openGOESClass = QAction("&GOES (Class)", self)
        openGOESClass.setShortcut("Ctrl+G")
        openGOESClass.setStatusTip('Open GOES (Class) Directory')
        openGOESClass.triggered.connect(self.goesClassOpen)        

        #GOES UAH
        openGOESUAH = QAction("&GOES (UAH)", self)
        openGOESUAH.setShortcut("Ctrl+U")
        openGOESUAH.setStatusTip('Open GOES (UAH) Directory')
        openGOESUAH.triggered.connect(self.goesOpen)        

        #NCEP/NCAR Reanalysis OpenDAP
        openNCEPNCAR = QAction("&NCEP/NCAR Reanalysis", self)
        openNCEPNCAR.setShortcut("Ctrl+N")
        openNCEPNCAR.setStatusTip('Open NCEPNCAR')
        openNCEPNCAR.triggered.connect(self.NcepNcar) 

        #MERRA OpenDAP
        openMERRA = QAction("&MERRA", self)
        openMERRA.setShortcut("Ctrl+M")
        openMERRA.setStatusTip('Open MERRA')
        openMERRA.triggered.connect(self.merraOpen)

        #NEXRAD Radar Amazon Web Services (Make sure you have an account set up)
        openNEXRAD = QAction("&NEXRAD",self)
        openNEXRAD.setShortcut("Ctrl+R")
        openNEXRAD.setStatusTip('Open NEXRAD')
        openNEXRAD.triggered.connect(self.nexradOpen)

        #CMAQ 
        openCMAQ = QAction("&CMAQ",self)
        openCMAQ.setShortcut("Ctrl+Q")
        openCMAQ.setStatusTip('Open CMAQ')
        openCMAQ.triggered.connect(self.cmaqOpen)

        #WRF met_em files
        openMET = QAction("&MET",self)
        openMET.setShortcut("Ctrl+I")
        openMET.setStatusTip('Open MET')
        openMET.triggered.connect(self.metOpen)

        ####End Dataset menu bar ##############################

        # Get EOM
        getEOM = QAction('&Eq. of Motion', self)
        getEOM.setStatusTip("Calculate EOM Components")
        getEOM.triggered.connect(self.EOMget)

        #Color Palettes
        self.colorlist = ['jet','brg','rainbow','bwr','YlOrRd','viridis','magma','gray','pyart_NWSRef']
        jet = QAction("&Default (jet)",self)
        jet.triggered.connect(lambda: self.selectionChangeColorPallete(0))
        brg = QAction(QIcon('brg.png'),"&BlueRedGreen",self)
        brg.triggered.connect(lambda: self.selectionChangeColorPallete(1))
        rainbow = QAction("&Rainbow",self)
        rainbow.triggered.connect(lambda: self.selectionChangeColorPallete(2))
        bwr = QAction("&BlueWhiteRed",self)
        bwr.triggered.connect(lambda: self.selectionChangeColorPallete(3))
        ylorrd = QAction("&YellowOrangeRed",self)
        ylorrd.triggered.connect(lambda: self.selectionChangeColorPallete(4))
        viridis = QAction("&Viridis",self)
        viridis.triggered.connect(lambda: self.selectionChangeColorPallete(5))
        magma = QAction("&Magma",self)
        magma.triggered.connect(lambda: self.selectionChangeColorPallete(6))
        gray = QAction("&Gray",self)
        gray.triggered.connect(lambda: self.selectionChangeColorPallete(7))
        ref = QAction("&NWS Reflectivity",self)
        ref.triggered.connect(lambda: self.selectionChangeColorPallete(8))

        #Create map background menu bar
        defaultClear = QAction("&Default (Clear)", self)
        defaultClear.triggered.connect(lambda: self.changeBackground(0))

        blueMarble = QAction("&Blue Marble", self)
        blueMarble.triggered.connect(lambda: self.changeBackground(1))

        shadedRelief = QAction("&Shaded Relief", self)
        shadedRelief.triggered.connect(lambda: self.changeBackground(2))

        topo = QAction("&Topo Map", self)
        topo.triggered.connect(lambda: self.changeBackground(3))

        coarse = QAction("&Coarse", self)
        coarse.triggered.connect(lambda: self.changeResolution(0))

        low = QAction("&Low", self)
        low.triggered.connect(lambda: self.changeResolution(1))

        intermediate = QAction("&Intermediate", self)
        intermediate.triggered.connect(lambda: self.changeResolution(2))

        high = QAction("&High", self)
        high.triggered.connect(lambda: self.changeResolution(3))

        full = QAction("&Full", self)
        full.triggered.connect(lambda: self.changeResolution(4))

        self.statusBar()

        mainMenu = self.menuBar()
        mainMenu.setStyleSheet(Layout.QMenuBar())
        fileMenu = mainMenu.addMenu('&File')
        fileMenu.setStyleSheet(Layout.QMenu())
        #fileMenu.addAction(openDir)
        fileMenu.addAction(extractAction)

        editorMenu = mainMenu.addMenu('&Editor')
        editorMenu.setStyleSheet(Layout.QMenu())
        editorMenu.addAction(openEditor)

        dataSetMenu = mainMenu.addMenu('&Add Dataset')
        dataSetMenu.setStyleSheet(Layout.QMenu())

        #WRF
        menuWRF = dataSetMenu.addMenu('&WRF')
        menuWRF.setStyleSheet(Layout.QMenu())
        menuWRF.addAction(openWRF)
        menuWRF.addAction(openMET)

        #GOES
        menuGOES = dataSetMenu.addMenu('&GOES')
        menuGOES.setStyleSheet(Layout.QMenu())
        menuGOES.addAction(openGOESClass)
        menuGOES.addAction(openGOESUAH)

        #Remaining Datasets
        dataSetMenu.addAction(openCMAQ)
        dataSetMenu.addAction(openNCEPNCAR)
        dataSetMenu.addAction(openMERRA)
        dataSetMenu.addAction(openNEXRAD)

        #set up the plot settings menu
        plotMenu = mainMenu.addMenu('&Plot Settings')
        plotMenu.setStyleSheet(Layout.QMenu())
        calcMenu = mainMenu.addMenu('&Calculations')
        calcMenu.setStyleSheet(Layout.QMenu())
        menuEOM = calcMenu.addAction(getEOM)

        #Map background options
        mapMenu = plotMenu.addMenu('&Change Map Background')
        mapMenu.setStyleSheet(Layout.QMenu())
        mapMenu.addAction(defaultClear)
        mapMenu.addAction(blueMarble)
        mapMenu.addAction(shadedRelief)
        mapMenu.addAction(topo)

        #Map Resolution options
        resMenu = plotMenu.addMenu('&Change Map Resolution')
        resMenu.setStyleSheet(Layout.QMenu())
        resMenu.addAction(coarse)
        resMenu.addAction(low)
        resMenu.addAction(intermediate)
        resMenu.addAction(high)
        resMenu.addAction(full)

        #Colorbar options
        colorbarMenu = plotMenu.addMenu('&Change Colorbar')
        colorbarMenu.setStyleSheet(Layout.QMenu())
        colorbarMenu.addAction(jet)
        colorbarMenu.addAction(brg)
        colorbarMenu.addAction(rainbow)
        colorbarMenu.addAction(bwr)
        colorbarMenu.addAction(ylorrd)
        colorbarMenu.addAction(viridis)
        colorbarMenu.addAction(magma)
        colorbarMenu.addAction(gray)
        colorbarMenu.addAction(ref)
        
        #Contour Fill options
        contourMenu = plotMenu.addMenu('&Contour Fill Type')
        contourMenu.setStyleSheet(Layout.QMenu())
        cf = QAction("ContourF",self)
        cf.triggered.connect(lambda: self.ContourFPlot())
        pc = QAction("PColorMesh",self)
        pc.triggered.connect(lambda: self.PColorMeshPlot())
        contourMenu.addAction(cf)
        contourMenu.addAction(pc)

        #Overlay options
        overlayMenu = plotMenu.addMenu('&Overlays')
        overlayMenu.setStyleSheet(Layout.QMenu())
        wb = QAction("Wind Barbs",self)
        wb.triggered.connect(lambda: self.plotWindBarbs())
        wb.setShortcut("Ctrl+B")
        wb.setStatusTip('Overlay Wind Barbs')
        wv = QAction("Wind Vectors",self)
        wv.triggered.connect(lambda: self.plotWindVectors())
        wv.setShortcut("Ctrl+V")
        wv.setStatusTip('Overlay Wind Vectors')
        c2 = QAction("2nd Contour",self)
        c2.triggered.connect(lambda: self.plotContour2())
        c2.setShortcut("Ctrl+C")
        c2.setStatusTip('Overlay 2nd Contour')
        overlayMenu.addAction(c2)
        overlayMenu.addAction(wb)
        overlayMenu.addAction(wv)
        
        #### End Menu ####
        self.dataSet = [None]
        self.paths = [None]
        self.path = None 
        self.addplot = None
        self.cbar = None
        self.numDset = 0
        self.on_draw(None)
        self.axes1 = []   
        self.cs = None
        self.filltype = "contourf"
        self.cmap = 'jet'
        self.plotbarbs = False
        self.plotvectors = False
        self.plotcontour2 = False
        self.background = None
        self.recallProjection = True
        self.changeColor = False
        self.cs = None
        self.cs2 = None
        self.cs2label = None
        self.barbs = None
        self.vectors = None
        self.vectors2 = None
        self.vectorkey = None
        self.resolution = 'l'
        self.clear = False
 
    #Select the directory
    def wrfOpen(self):
        self.dname = 'WRF'
        #self.selectdataset.close()
        self.path = QFileDialog.getExistingDirectory(self, 'Select Directory')
        self.prefix = 'wrfout*'
        files = sorted(glob.glob(os.path.join(str(self.path),self.prefix)))
        #print(os.path.join(str(self.path),self.prefix))
        if self.path:
            if len(files) > 0:
                self.numDset += 1
                self.dataSet.append(WrfDataset.WrfDataset())
                self.paths.append(self.path)
                QApplication.setOverrideCursor(Qt.WaitCursor)
                self.dataSet[self.numDset].name(path=str(self.path) , prefix='wrfout')
                QApplication.restoreOverrideCursor()
            else:
                self.errorNoFilesFound()

    #Select the directory
    def metOpen(self):
        self.dname = 'MET'
        #self.selectdataset.close()
        self.path = QFileDialog.getExistingDirectory(self, 'Select Directory')
        self.prefix = 'met*em*'
        files = sorted(glob.glob(os.path.join(str(self.path),self.prefix)))
        #print(os.path.join(str(self.path),self.prefix))
        if self.path:
            if len(files) > 0:
                self.numDset += 1
                self.dataSet.append(METDataset.METDataset())
                self.paths.append(self.path)
                QApplication.setOverrideCursor(Qt.WaitCursor)
                self.dataSet[self.numDset].name(path=str(self.path) , prefix='met*em')
                QApplication.restoreOverrideCursor()
            else:
                self.errorNoFilesFound()

    #Select the directory
    def goesClassOpen(self):
        self.dname = 'GOES_Class'
        #self.selectdataset.close()
        self.path = QFileDialog.getExistingDirectory(self, 'Select Directory')
        self.prefix = 'goes*'
        files = sorted(glob.glob(os.path.join(str(self.path),self.prefix)))
        if self.path:
            if len(files) > 0:
                self.numDset += 1
                self.dataSet.append(GOESClass.GOESClassDataset())
                self.paths.append("GOES Class")
                QApplication.setOverrideCursor(Qt.WaitCursor)
                self.dataSet[self.numDset].name(path=str(self.path) , prefix='goes')
                QApplication.restoreOverrideCursor()
            else:
                self.errorNoFilesFound()

    #Select the directory
    def goesOpen(self):
        self.dname = 'GOES_UAH'
        #self.selectdataset.close()
        self.path = QFileDialog.getExistingDirectory(self, 'Select Directory')
        self.prefix = 'wrf??km*'
        files = sorted(glob.glob(os.path.join(str(self.path),self.prefix)))
        if self.path:
            if len(files) > 0:
                self.numDset += 1
                self.dataSet.append(GOES.GOESDataset())
                self.paths.append("GOES UAH")
                QApplication.setOverrideCursor(Qt.WaitCursor)
                self.dataSet[self.numDset].name(path=str(self.path) , prefix='wrf??km')
                QApplication.restoreOverrideCursor()
            else:
                self.errorNoFilesFound()

    #Select the directory
    def merraOpen(self):
        self.dname = 'MERRA'
        #self.selectdataset.close()
        self.numDset += 1
        QApplication.setOverrideCursor(Qt.WaitCursor)
        self.dataSet.append(MERRA.MERRAdataset())
        QApplication.restoreOverrideCursor()

    def NcepNcar(self):
        self.dname = 'NCEP'
        #self.selectdataset.close()
        self.numDset += 1
        QApplication.setOverrideCursor(Qt.WaitCursor)
        self.dataSet.append(NCAR.NCARdataset())
        QApplication.restoreOverrideCursor()

    def nexradOpen(self):
        self.dname = 'NEXRAD Radar'
        self.numDset += 1
        QApplication.setOverrideCursor(Qt.WaitCursor)
        self.dataSet.append(NEXRAD.Radardataset())
        QApplication.restoreOverrideCursor()

    #Select the directory
    def cmaqOpen(self):
        self.dname = 'CMAQ'
        #self.selectdataset.close()
        self.path = QFileDialog.getExistingDirectory(self, 'Select Directory')
        self.prefix = 'CCTM*'
        files = sorted(glob.glob(os.path.join(str(self.path),self.prefix)))
        if self.path:
            if len(files) > 0:
                self.numDset += 1
                self.dataSet.append(CmaqDataset.CmaqDataset())
                self.paths.append(self.path)
                QApplication.setOverrideCursor(Qt.WaitCursor)
                self.dataSet[self.numDset].name(path=str(self.path) , prefix='CCTM')
                QApplication.restoreOverrideCursor()
            else:
                self.errorNoFilesFound()

    def EOMget(self):
        if self.dataSet[self.numDset].dsetname == "MERRA":
            self.eom = True
            self.recallProjection = False
            self.on_draw(self.plotCount)
        else:
            self.EOMerror()
    def selectionChangeColorPallete(self,i):
        self.changeColor = True
        self.cmap = self.colorlist[i]
        self.on_draw(self.plotCount)

    #Function that changes the background color option when user selected
    def changeBackground(self,i):
        self.ColorBar = None
        self.recallProjection = True
        if i == 0:
            self.background = None
            self.clear = True
            #self.axes1[self.slbplt.pNum-1] = None
            #self.slbplt.figure.clear()
        elif i == 1:
            self.background = '.bluemarble'
        elif i == 2:
            self.background = '.shadedrelief'
        else:
            self.background = '.etopo'
        self.on_draw(self.plotCount)

    def changeResolution(self,i):
        self.recallProjection = True
        if i == 0:
            self.resolution = 'c'
        elif i == 1:
            self.resolution = 'l'
        elif i == 2:
            self.resolution = 'i'
        elif i == 3:
            self.resolution = 'h'
        else:
            self.resolution = 'f'
        self.on_draw(self.plotCount)

    #These two functions change the plot settings
    def PColorMeshPlot(self):
        if self.filltype == "contourf":
            if self.cs != None:
                for coll in self.cs.collections:
                    coll.remove()
        self.cs = None
        self.domain_average = None
        self.filltype = "pcolormesh"
        self.on_draw(self.plotCount)

    def ContourFPlot(self):
        self.cs = None
        self.domain_average = None
        self.filltype = "contourf"
        self.on_draw(self.plotCount)

    def plotWindBarbs(self):
        if self.plotbarbs == False:
            self.plotvectors = False
            if self.vectors != None:
                self.vectors.remove()
            if self.vectorkey != None:
                self.vectorkey.remove()
            self.vectors = None
            self.vectorkey = None
            self.plotbarbs = True
        else:
            self.plotvectors = False
            self.plotbarbs = False
            if self.barbs != None:
                self.barbs.remove()
                self.vectors2.remove()
            self.barbs = None
        self.on_draw(self.plotCount)

    def plotWindVectors(self):
        if self.plotvectors == False:
            self.plotbarbs = False
            if self.barbs != None:
                self.barbs.remove()
                self.vectors2.remove()
            self.plotvectors = True
            self.barbs = None
        else:
            self.plotvectors = False
            if self.vectors != None:
                self.vectors.remove()
            if self.vectorkey != None:
                self.vectorkey.remove()
            self.vectors = None
            self.vectorkey = None
        self.on_draw(self.plotCount)

    def plotContour2(self):
        if self.plotcontour2 == False:
            self.plotcontour2 = True
        else:
            self.plotcontour2 = False
            if self.cs2 != None:
                for coll in self.cs2.collections:
                    coll.remove()
                for labels in self.cs2label:
                    labels.remove()

            self.cs2 = None
        self.on_draw(self.plotCount)

    #Create exit pop up window    
    def close_program(self):
       
        choice = QMessageBox.question(self, 'Quit',"Exit the program?",
                                      QMessageBox.Yes | QMessageBox.No)
        
        if choice == QMessageBox.Yes:
            sys.exit()
        else:
            pass  
    
    #Create error message when files are not found in selected directory
    def errorNoFilesFound(self):
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Information)
        msg.setText('No '+self.dname+' files found in the selected directory.\n\n'+ 
                    self.dname+' files must be named '+self.prefix)
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def EOMerror(self):
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Information)
        msg.setText('Equations of Motion calculations are only available in the MERRA dataset at this time.')
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()


    def readfxn(self,wrfobj,varname):
        print('numGrids = ',wrfobj.numGrids)
        phyd = wrfobj.readNCVariable('P_HYD')
        self.dataSet.description = 'Hyrostatic Pressure'
        self.dataSet.units = 'hPa'
        print("dx1, dx2=",wrfobj.dx[0],
              wrfobj.dx[1],wrfobj.glats[0].shape,wrfobj.ny[0])
        return phyd/100.0
        
    def create_main_frame(self):
    
        self.main_frame = QWidget()
        qlayout = QHBoxLayout(self.main_frame)
        self.main_frame.setLayout(qlayout)
        #self.setGeometry(100,100,1200,750)

        #Get Desktop Information
        desktop = QDesktopWidget()
        screen_size = QRectF(desktop.screenGeometry(desktop.primaryScreen()))
        if (screen_size.width() > 3000.):
            self.screenx = screen_size.x() + screen_size.width()/2.
            self.screeny = screen_size.y() + screen_size.height() 
        else:
            self.screenx = screen_size.x() + screen_size.width()
            self.screeny = screen_size.y() + screen_size.height()
       
        #Scale GUI Window
        #self.resize(1200,750) 
        self.resize(self.screenx*.9,self.screeny*.85) 
        #self.setMinimumWidth(self.screenx*0.8) 

        self.controlPanel = QGroupBox()
        self.controlPanel.setStyleSheet(Layout.QGroupBox())
#        self.controlPanel.setMinimumWidth(.35*.8*self.screenx)
#        self.controlPanel.setMaximumWidth(.35*.8*self.screenx)
        self.controlPanel.setMinimumWidth(.40*.8*self.screenx)
        self.controlPanel.setMaximumWidth(.40*.8*self.screenx)
	
        self.cpLayout = QVBoxLayout(self.controlPanel)
        #self.DsetButton = QPushButton('Add Dataset', self)
        #self.DsetButton.resize(self.DsetButton.minimumSizeHint())
        #self.DsetButton.clicked.connect(self.addDataSetAction)
       
        selectPlotControlWidget = QWidget()
        selectPlotControlWidgetLayout = QHBoxLayout()
        selectPlotControlWidget.setLayout(selectPlotControlWidgetLayout)
        self.addButton = QPushButton('Add Plot', self)
        self.addButton.setStyleSheet(Layout.QPushButton())
        self.addButton.resize(self.addButton.minimumSizeHint())
        self.addButton.clicked.connect(self.addButtonAction)
        self.deleteButton = QPushButton('Delete Plot', self)
        self.deleteButton.setStyleSheet(Layout.QPushButton())
        self.deleteButton.resize(self.deleteButton.minimumSizeHint())
        self.deleteButton.clicked.connect(self.deleteButtonAction)
        selectPlotControlWidgetLayout.addWidget(self.addButton)
        selectPlotControlWidgetLayout.addWidget(self.deleteButton)
        
        #self.cpLayout.addWidget(self.DsetButton)
        self.cpLayout.addWidget(selectPlotControlWidget)
        self.plotControlTabs = QTabWidget(self.controlPanel)
        self.plotTab.append(QWidget())
        #self.plotControlTabs.addTab(self.plotTab[0],'GridView')
        self.cpLayout.addWidget(self.plotControlTabs)
 
        qscroll = QScrollArea(self.main_frame)
        qscroll.setStyleSheet(Layout.QScrollArea())
        qlayout.addWidget(qscroll)
        qlayout.addWidget(self.controlPanel)

        qscrollContents = QWidget()
        qscrollLayout = QVBoxLayout(qscrollContents)

        qscroll.setWidget(qscrollContents)
        qscroll.setWidgetResizable(True)
        
        qfigWidget = QWidget(qscrollContents)
        
        self.pltList = []
        self.plotLayout = QVBoxLayout()
        qfigWidget.setLayout(self.plotLayout)
        qscrollLayout.addWidget(qfigWidget)
        qscrollContents.setLayout(qscrollLayout)
        self.setCentralWidget(self.main_frame)

        self.move(self.screenx*.05,self.screeny*.025)
        self.setStyleSheet(Layout.QMain())

    def get_data2(self):
        return np.arange(20).reshape([4, 5]).copy()
    
    def addButtonAction(self):
        if len(self.dataSet) >= 2:
            self.plotCount = self.plotCount + 1
            self.cw = CanvasWidget(self,self.plotCount)
            slbplt = PlotSlab(dset=self.dataSet,AppWid=self)
            slbplt.setConnection(self.on_draw,self.plotCount)
            slbplt.plotCount = self.plotCount
            self.plotTab.append(QWidget())
            self.tabLayout = QVBoxLayout()
            self.plotTab[self.plotCount].setLayout(self.tabLayout)
            slbplt.getControlBar()
            self.cw.setPlot(slbplt)
            self.plotLayout.addWidget(self.cw)
            self.pltList.append(self.cw)
            tabName = 'Plot #'+ str(self.plotCount)
            self.plotControlTabs.addTab(self.plotTab[self.plotCount],tabName)
            self.plotControlTabs.setCurrentIndex(self.plotCount)
            self.on_draw(self.plotCount)
        else:
            self.errorDset()

    def deleteButtonAction(self):
        if (self.plotCount != 0):
            self.plotLayout.itemAt(self.plotControlTabs.currentIndex()).widget().deleteLater()
            #self.plotTab[self.plotControlTabs.currentIndex()] = None
            #self.pltList[self.plotControlTabs.currentIndex()] = None
            self.plotControlTabs.removeTab(self.plotControlTabs.currentIndex())
            self.on_draw(self.plotCount)

    def on_draw(self,plotnum):
        if len(self.pltList) > 0:
            if self.dname != 'NEXRAD Radar':
                self.pltList[plotnum-1].drawPlot()
            else:
                self.pltList[plotnum-1].drawPlot()

        self.main_frame.show()
        
    def on_key_press(self, event):
        print('you pressed', event.key)
        #implement the default mpl key press events described at
        # http://matplotlib.org/users/navigation_toolbar.html#navigation-keyboard-shortcuts
        key_press_handler(event, self.canvas, self.mpl_toolbar)

    def errorDset(self):
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Information)
        msg.setText("Please add a dataset before trying to create a plot.")
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

###############################################################################
####                           End AppForm() Object                        ####
###############################################################################

###############################################################################
####                       Begin DataSelector() Object                     ####
###############################################################################

class DataSelector:
    def testfxn(self,event):
        print ('leave_figure', event.canvas.figure)
        event.canvas.figure.patch.set_facecolor('blue')
    def __init__(self, line):
        self.points = 0
        self.startMark = None
        self.visible = False
        self.line = line
        self.userSelects = []
        self.numSelection = 0
        self.mode = 'None'
        self.modeDefs = {'s':'Select', 'S':'Select','d':'Delete','D':'Delete',
                         'v':'View','V': 'View', 'escape': 'None'}
        self.subMode = 'None'
        self.subModeDefs = {'l':'Line', 'L':'Line','p':'Polygon','P':'Polygon'}
        #self.xs = list(line.get_xdata())
        #self.ys = list(line.get_ydata())
        self.xs = []
        self.ys = []
        #print('init',self.xs)
        line.figure.canvas.setFocusPolicy( Qt.ClickFocus )
        line.figure.canvas.setFocus()
        self.cid = line.figure.canvas.mpl_connect('button_press_event',self)
        self.cid = line.figure.canvas.mpl_connect('key_press_event', self.keyPress)
        
    def keyPress(self, event):
        print('in keypress')
        if event.key in self.modeDefs and 1:
            self.mode = self.modeDefs[event.key]
            print('Mode = ',self.mode)
        if event.key in self.subModeDefs and self.mode != 'None':
            self.subMode = self.subModeDefs[event.key]
            print('sub Mode = ',self.subMode)
    
    def clearSelect(self):
        
        if(self.p1 != None):
            self.p1.remove()
            
        if(self.p2 != None):
            self.p2.remove()
        
        self.line.set_data(None, None)

        self.points = 0
        
        self.line.figure.canvas.draw()
        
    def addToLine(self,event):
        
        self.xs.append(event.xdata)
        
        self.ys.append(event.ydata)
        
        if(self.points == 0 and self.subMode == 'Line'):
            self.p1= event.inaxes.annotate("A", xy=(event.xdata,event.ydata), xycoords=("data"))
        else:
            if(self.points == 0 and self.subMode == 'Polygon'):
               self.startMark= event.inaxes.plot(event.xdata,event.ydata,'o')
               print(type(self.startMark),self.startMark)
            self.line.set_data(self.xs, self.ys)
              
        self.points = self.points+1

    def endLine(self,event):
    
        if (self.subMode == 'Line'):
            self.p2=event.inaxes.annotate("B", xy=(event.xdata,event.ydata),
                                          xycoords=("data"))   
        else:
            self.startMark[0].remove()
            self.xs.append(self.xs[0])
            self.ys.append(self.ys[0])
            self.line.set_data(self.xs, self.ys)
        self.userSelects.append((self.xs,self.ys))
        
        self.numSelection = self.numSelection+1

    def __call__(self, event):
        
        if event.inaxes!=self.line.axes: return
        if (self.mode == 'None'):
            if event.key in self.modeDefs:
                self.mode = self.modeDefs[event.key]
                print('Mode = ',self.mode)
        if(self.mode == 'Select'):
            if(event.button == 3 ):
                print("Erasing line")
                self.clearSelect()
            else:
                if(event.dblclick ):
                    self.endLine(event)
                else:
                    self.addToLine(event)
        
        self.line.figure.canvas.draw()
       
def main():
    app = QApplication(sys.argv)
    app.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    form = AppForm()
    form.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
