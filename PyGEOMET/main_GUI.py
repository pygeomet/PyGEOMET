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
import matplotlib.path as mplPath
import scipy.spatial as spatial
import netCDF4
import PyGEOMET.utils.wrf_functions as wrf
import PyGEOMET.utils.wrf_cython as wrf_cython
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
import PyGEOMET.datasets.SoundingDataset as SoundDataset
import PyGEOMET.datasets.NetCDFDataset as NetCDFDataset
import PyGEOMET.datasets.GOESRDataset as GOESRDataset
import PyGEOMET.utils.LayoutFormat as Layout
import PyGEOMET.utils.sensor_info as CRTMInfo
import PyGEOMET.utils.create_colormap as create_colormap
import PyGEOMET.utils.readObs as readObs
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

        #Horizontal Plot
        if (self.plotObj.currentPType == 'Horizontal Slice' or 
            self.plotObj.currentPType == '1-panel'):
            self.plotObj.plot()
            if self.plotObj.appobj.eom == True:
                if self.plotObj.cid is None:
                    self.plotObj.cid = self.canvas.mpl_connect('button_press_event', self.onclick)
                if self.plotObj.col is not None and self.plotObj.row is not None:
                    self.plotObj.getTable(self.plotObj.col,self.plotObj.row)
            self.canvas.draw()
        
        #Vertical Cross section
        if (self.plotObj.currentPType == 'Vertical Slice'):
            self.plotObj.plot()
            self.canvas.draw()

        #Skew T
        if (self.plotObj.currentPType == 'SkewT/Hodograph'):
            if self.plotObj.replot2d == True:
                self.plotObj.plot()
                self.canvas.draw()
            if self.plotObj.cid is None: 
                self.plotObj.cid = self.canvas.mpl_connect('button_press_event', self.onclick)
            if self.plotObj.col is not None and self.plotObj.row is not None:
                self.plotObj.plotSkewT(self.plotObj.col,self.plotObj.row)
                plt.ion()
                plt.show()  

        #Vertical Profile
        if (self.plotObj.currentPType == 'Vertical Profile'):
            if self.plotObj.replot2d == True:
                self.plotObj.plot()
                self.canvas.draw()
            if self.plotObj.cid is None:
                self.plotObj.cid = self.canvas.mpl_connect('button_press_event', self.onclick)
            if self.plotObj.col is not None and self.plotObj.row is not None:
                self.plotObj.plotVertPro(self.plotObj.col,self.plotObj.row)
                plt.ion()
                plt.show()

        #Time Series
        if (self.plotObj.currentPType == 'Time Series'):
            if self.plotObj.replot2d == True:
                self.plotObj.plot()
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
        if (self.plotObj.currentPType == 'Difference Plot'):
            self.plotObj.plot()
            self.canvas.draw()

        #Spatial average
        if (self.plotObj.currentPType == 'Spatial Stats'):
            self.plotObj.spatialStats()
            self.plotObj.plot()
            self.canvas.draw()

        #Anomaly
        if (self.plotObj.currentPType == 'Anomaly'):
            self.plotObj.calculateAnomaly()
            if (self.plotObj.cAnomaly == True):
                self.plotObj.plot()
                self.canvas.draw()

        #Surfae Statistics
        if (self.plotObj.currentPType == 'Surface Statistics'):
            self.plotObj.sfcStats()
            self.plotObj.plot()
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
        #Make sure the user clicks on the map
        #If they don't, ix and iy will be None. 
        #This will result in an error in the row/col calculations
        if (ix != None and iy != None):
            if (self.plotObj.dataSet.runType != 'IDEAL'):
                lon, lat  = self.plotObj.dataSet.map[self.plotObj.currentGrid-1](
                            self.plotObj.dataSet.glons[self.plotObj.currentGrid-1],
                            self.plotObj.dataSet.glats[self.plotObj.currentGrid-1])
                if (self.plotObj.dataSet.dsetname == 'Sounding'):
                    a = np.vstack((lat,lon)).T
                    pt = [iy, ix]
                    distance, index = spatial.KDTree(a).query(pt)
                    #print(distance)
                    #print(index)
                    #print(a[index])
                    #print(iy,ix)
                    col = index
                    row = index
                    print(self.plotObj.dataSet.stat_id[index])
                else:                    
                    clon,clat = self.plotObj.dataSet.map[self.plotObj.currentGrid-1](
                                self.plotObj.dataSet.lon0[self.plotObj.currentGrid-1],
                                self.plotObj.dataSet.lat0[self.plotObj.currentGrid-1])
                    tmp = np.argsort(np.abs(lat[:,int(self.plotObj.dataSet.nx[self.plotObj.currentGrid -1]/2)] - clat))[0]
                    tmp2 = np.argsort(np.abs(lon[tmp,:] - clon))[0]
                    row = np.argsort(np.abs(lat[:,tmp2] - iy))[0]
                    col  = np.argsort(np.abs(lon[row,:] - ix))[0]
            else:                
                
                row = int(iy*1000/self.plotObj.dataSet.dy[self.plotObj.currentGrid-1] +
                         (self.plotObj.dataSet.ny[self.plotObj.currentGrid-1]-1)/2)
                col = int(ix*1000/self.plotObj.dataSet.dx[self.plotObj.currentGrid-1] +
                         (self.plotObj.dataSet.nx[self.plotObj.currentGrid-1]-1)/2)

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
        else:
            self.errorClick()

    #Throw an error for a bad click
    def errorClick(self):
        msg = QMessageBox(self.plotObj.appobj)
        msg.setIcon(QMessageBox.Information)
        msg.setText('Make sure you clicked on a point within the map')
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()
 
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
        self.currentPType = self.dataSet.ptypes[0]
        self.currentOrient = 0
        self.nz = None
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
        self.calTimeRange = True
        self.colorbox = None
        self.vslicebox = None
        self.orientList = ['xz','yz','custom']
        self.appobj = AppWid
        self.filltype = "contourf"
        self.cmap2 = 'jet'
        self.plotbarbs = False
        self.plotvectors = False
        self.plotcontour2 = False
        self.selectedVvar = None
        self.selecteddVvar = None
        self.vplot = None
        self.vertControl = None
        self.pControl = None
        self.diffControl = None
        self.seriesControl = None
        self.selectedParcel = 0
        self.skewParcel = 'SB'
        self.background = None
        self.recallProjection = True
        self.clear = False
        self.changeColor = False
        self.ColorBar = None
        self.extend = 'both'
        self.plotcoasts = True
        self.plotcountries = True
        self.plotstates = True
        self.plotcounties = False
        self.plotcities = False
        self.plotlatlon = True
        self.cs = None
        self.cs2 = None
        self.cs2label = None
        self.barbs = None
        self.vectors = None
        self.vectorkey = None
        self.currentGrid = 1
        self.resolution = self.dataSet.resolution
        self.changeRes = False
        self.domain_average = None
        self.coasts = None
        self.countries = None
        self.states = None
        self.counties = None
        self.cities = False
        self.parallels = None
        self.meridians = None
        self.cmap = self.cmap2
        self.appobj.eom = None
        self.colorlock = False
        self.crtmbox = None
        self.crtmSensor = None
        self.crtmChannel = None
        self.crtmVar = None
        self.userGrid = False
        self.east = None
        self.north = None
        self.west = None
        self.south = None
        #Previous plot type (used to create vertical difference plot)
        self.prevPType = None
        self.yinterval = 1
        self.xinterval = 1
        self.max_val = None
        self.min_val = None
        #Set spatial statistics types
        self.spatialOptions = ['Mean','Max','Min']
        self.spatialControl = None
        self.statsControl = None
        self.statsChangeTime = True
        self.selectedStatType = 'Bias'

        #self.cPIndex = self.appobj.plotControlTabs.currentIndex()
        if 'runType' not in dir(self.dataSet):
            self.dataSet.runType = None

    def setFigure(self, fig):
        if(fig != None):
            self.figure = fig

    def setConnection(self,fxn=None,plotnum=None):
        if(fxn != None):
            self.pltFxn = fxn
        if(plotnum != None):
            self.pNum = plotnum
           
    def readField(self):
        #Updata data here to accomidate multiple grids with multiple plots
        self.dataSet.setGrid(self.currentGrid)
        if not self.derivedVar:
            if self.currentVar != None:
                ##Updata data here to accomidate multiple grids with multiple plots
                #self.dataSet.setGrid(self.currentGrid)
                self.var = self.dataSet.readNCVariable(self.dataSet.variableList[self.currentVar],
                    barbs=self.plotbarbs, vectors = self.plotvectors,
                    contour2=self.plotcontour2)
                self.varTitle = self.dataSet.description +' ('+self.dataSet.units
                self.varTitle = self.varTitle +') \n'+self.dataSet.getTime()
                #account for unstaggering of the grids in 3-dimensional variables
                if self.dataSet.dsetname == 'WRF' or self.dataSet.dsetname == 'MET':
                    if len(self.var.shape) == 3:
                        if self.var.shape[0] == self.dataSet.nz[self.currentGrid-1]:
                            self.var = wrf.unstaggerZ(self.var)
                        if self.var.shape[1] == self.dataSet.ny[self.currentGrid-1]:
                            self.var = wrf.unstaggerY(self.var)
                        if self.var.shape[2] == self. dataSet.nx[self.currentGrid-1]:
                            self.var = wrf.unstaggerX(self.var) 
                        self.varTitle = self.varTitle + ', Level=' + str(self.dataSet.levelList[self.currentLevel])
                    else:
                        if self.var.shape[0] == self.dataSet.ny[self.currentGrid-1]:
                            plotself.var = wrf.unstaggerY(self.var)
                        if self.var.shape[1] == self.dataSet.nx[self.currentGrid-1]:
                            self.var = wrf.unstaggerX(self.var)
                self.diffvar = self.var
            else:
                pass
        else:
            if self.currentdVar != None:
                ##Updata data here to accomidate multiple grids with multiple plots
                #self.dataSet.setGrid(self.currentGrid)
                t0 = time.clock()
                t1 = time.time()
                QApplication.setOverrideCursor(Qt.WaitCursor)
                dvar = wrf_dvar.WRFDerivedVar(dset = self.dataSet, 
                                              var = self.dataSet.dvarlist[self.currentdVar],
                                              ptype = self.currentPType, sensor = self.crtmSensor,
                                              channel = self.crtmChannel, 
                                              path = self.appobj.main_path,
                                              req_var = self.crtmVar)
                QApplication.restoreOverrideCursor()
                print( time.clock() - t0,"seconds process time plots")
                print( time.time() - t1,"seconds wall time plots")

                self.var = dvar.var
                self.var2 = dvar.var2
                self.varTitle = dvar.varTitle
                #Short title used in time series plot
                self.sTitle = dvar.sTitle
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
            #Updata data here to accomidate multiple grids with multiple plots
            self.dataSet.setGrid(self.currentGrid)
            ind = np.where(np.array(self.diffdata.variableList) == self.dataSet.variableList[self.currentVar])[0]
            self.diffvar = self.diffdata.readNCVariable(self.diffdata.variableList[ind],
                barbs=self.plotbarbs, vectors = self.plotvectors,
                contour2=self.plotcontour2)
            #account for unstaggering of the grids in 3-dimensional variables
            if self.dataSet.dsetname == 'WRF' or self.dataSet.dsetname == 'MET':
                if len(self.diffvar.shape) == 3:
                    if self.diffvar.shape[0] == self.diffdata.nz[self.currentGrid-1]:
                        self.diffvar = wrf.unstaggerZ(self.diffvar)
                    if self.diffvar.shape[1] == self.diffdata.ny[self.currentGrid-1]:
                        self.diffvar = wrf.unstaggerY(self.diffvar)
                    if self.diffvar.shape[2] == self. diffdata.nx[self.currentGrid-1]:
                        self.diffvar = wrf.unstaggerX(self.diffvar)
                else:
                    if self.diffvar.shape[0] == self.diffdata.ny[self.currentGrid-1]:
                        self.diffvar = wrf.unstaggerY(self.diffvar)
                    if self.diffvar.shape[1] == self.diffdata.nx[self.currentGrid-1]:
                        self.diffvar = wrf.unstaggerX(self.diffvar)
        else:
            #Updata data here to accomidate multiple grids with multiple plots
            self.dataSet.setGrid(self.currentGrid)
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
        #Initialize the options for the dataset's
        # first plot type
        self.initializePlotOptions()
 
    def lockColorBar(self):
        self.colorlock = True
  
    def unlockColorBar(self):
        self.colorlock = False

    def controlColorBar(self):

        boxTitleString = 'Colorbar Control'
        self.colorbox = QGroupBox()
        self.colorbox.setStyleSheet(Layout.QGroupBox())

        colorboxLayout = QVBoxLayout()
        self.colorbox.setLayout(colorboxLayout)

        radio = QWidget()
        radiolayout = QHBoxLayout()
        radio.setLayout(radiolayout)
        self.lock = QRadioButton("Lock Color Range")
        self.lock.setStyleSheet(Layout.QCheckBox())
        self.lock.clicked.connect(lambda:self.lockColorBar())
        self.unlock = QRadioButton("Unlock Color Range")
        self.unlock.setStyleSheet(Layout.QCheckBox())
        self.unlock.clicked.connect(lambda:self.unlockColorBar())
        self.unlock.setChecked(True)
        radiolayout.addWidget(self.lock)
        radiolayout.addWidget(self.unlock)        

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
        if self.dataSet.runType == 'IDEAL':
            self.refboxLabel.setText('y-displacement:')
        else:
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
        import time
        t0 = time.clock()
        t1 = time.time()        
        #Check if colormap was changed
        if self.changeColor == True:
            #Need to have the values set = satellite color tables
            #Need to be brightness temperature or radiance
            #Check if derived var is set first because not every data set has
            # a derived variable list
            if (self.max_val != None and self.min_val != None
                and (self.currentdVar != None or self.dataSet.dsetname == 'GOES R' 
                or self.dataSet.dsetname == 'GOES Class' or self.dataSet.dsetname == 'GOES UAH')):
                #Next check if brightness temp/radiance
                if (self.currentdVar != None):
                    if (self.dataSet.dvarlist[self.currentdVar] == 'BrightTemp/Radiance'):
                        self.colormin = self.min_val
                        self.colormax = self.max_val
                        self.ncontours = 41
                        self.colorlock = True
                        self.lock.setChecked(True)
                if (self.dataSet.dsetname == 'GOES R' or self.dataSet.dsetname == 'GOES Class' 
                    or self.dataSet.dsetname == 'GOES UAH'):
                        self.colormin = self.min_val
                        self.colormax = self.max_val
                        self.ncontours = 41
                        self.colorlock = True
                        self.lock.setChecked(True)
            self.cmap = self.cmap2      

        #Set background map
        alpha = 1
        #Used for switching back to clear background
        if self.clear:
            self.figure.clear()
            self.clear = False
            #self.appobj.cs = None
            self.cs = None
            self.coasts = None
            self.states = None
            self.countries = None
            self.counties = None
            self.meridians = None
            self.parallels = None
            self.ColorBar = None
            self.domain_average = None

        #Set map background
        if self.background is not None:
            pre = 'self.dataSet.map['
            gridnum = str(self.currentGrid-1)+']'
            end = self.background+'(ax=self.axes1['+str(self.pNum-1)+'])'
            exec(pre+gridnum+end,None,None)
            #Make variable transparent
            alpha=0.75

        #Recall projection when grid/dataset is changed.
        if self.recallProjection == True: 
            if (self.currentPType == 'Vertical Slice' or (self.currentPType == 'Difference Plot' and
                self.prevPType == 'Vertical Slice')):
                self.figure.clear()
                #Set things to none so we don't try to remove them later
                self.cs = None
                self.cs2 = None
                self.barbs = None
                self.vectors = None
                self.vectorkey = None
                self.domain_average = None
                self.ColorBar = None
            self.axes1 = self.figure.add_subplot(111)
            self.dataSet.resolution = self.resolution
            if (self.dataSet.map[self.currentGrid-1] != None and self.changeRes == False and self.userGrid == False):
                self.dataSet.map[self.currentGrid-1].ax = self.axes1
            else:
                #Reset projection to change map resolution and define map boundaries
                #Note self.currentGrid-1 is not used because setProjection subtracts 1 from input
                self.dataSet.setProjection(self.currentGrid,axs=self.axes1,west=self.west,north=self.north,
                                           east=self.east,south=self.south)
                #Reset to false so set projection isnt called again
                self.userGrid = False
                self.changeRes = False

        #Set/Force NEXRAD settings
        if self.dataSet.dsetname == 'NEXRAD Radar':
            self.nz = 1
            self.filltype = 'pcolormesh' 
            #plot the parallels and meridians each time
            self.parallels = None
            self.meridians = None       
            #Based on called variable color and range from pyART
            self.cmap = self.dataSet.cmap
            self.colormin = self.dataSet.range[0]
            self.colormax = self.dataSet.range[1]
            self.ncontours = 41.

        #Remove Geography - have to do it before contouring otherwise it
        #                   won't be visible
        # Also, have to do it outside of var plotting loop to allow removal
        #       of geography before a variable is selected
        if self.coasts != None :
            self.coasts.remove()
        if self.countries != None:
            self.countries.remove()
        if self.states != None:
            self.states.remove()
        if self.counties != None:
            self.counties.remove()

        #Determine if a variable was passed - plots grid if not (technically you dont have to select a variable
        # to set currentVar of currentdVar in order to still create a surface stat plot)
        if (self.currentVar != None or self.currentdVar != None 
            or self.currentPType == 'Surface Statistics'):
            #Determine if the passed variable is 3D
            if (self.nz == 1 and self.currentPType != 'Surface Statistics'):
                #Make sure the variable is 3D for vertical cross section
                if (self.currentPType == 'Vertical Slice' or (self.currentPType == 'Difference Plot' and
                    self.prevPType == 'Vertical Slice')):
                    self.error3DVar()
                elif (self.currentPType == 'Difference Plot'):
                    pltfld = self.var - self.diffvar
                elif (self.currentPType == 'Spatial Stats'):
                    #Determine the stat type
                    statType = self.spatialOptions[self.spatialSType.currentIndex()]
                    if (statType == 'Mean'):
                        pltfld = np.nanmean(self.statVar,axis=0)
                    elif (statType == 'Max'):
                        pltfld = np.nanmax(self.statVar,axis=0)
                    else:
                        pltfld = np.nanmin(self.statVar,axis=0)
                elif (self.currentPType == 'Anomaly'):
                    pltfld = self.var - self.meanAnomaly
                else:
                    pltfld = self.var
            elif (self.nz > 1 and self.currentPType != 'Surface Statistics'):
                #If vertical cross section take the whole array
                if (self.currentPType == 'Vertical Slice'):
                    var = self.var
                #Create vertical slice difference plot
                elif (self.currentPType == 'Difference Plot' and self.prevPType == 'Vertical Slice'):
                    var = self.var - self.diffvar
                elif (self.currentPType == 'Difference Plot'):
                    pltfld = self.var[self.currentLevel] - self.diffvar[self.currentLevel]
                elif (self.currentPType == 'Spatial Stats'):
                    #Determine the stat type
                    statType = self.spatialOptions[self.spatialSType.currentIndex()]
                    #Setup a dummy variable to get the 3d (t,y,x) variable from 4d (t,z,y,x)
                    dum = np.squeeze(self.statVar[:,self.currentLevel,:,:])
                    if (statType == 'Mean'):
                        pltfld = np.nanmean(dum,axis=0)
                    elif (statType == 'Max'):
                        pltfld = np.nanmax(dum,axis=0)                    
                    else:
                        pltfld = np.nanmin(dum,axis=0)
                elif (self.currentPType == 'Anomaly'):
                    pltfld = self.var[self.currentLevel] - self.meanAnomaly
                else:
                    pltfld = self.var[self.currentLevel]
            elif (self.currentPType == 'Surface Statistics'):
                pltfld = self.stat
 
            #Set time string for saving the figure
            # timeObj is None in the NEXRAD dataset
            if (self.dataSet.dsetname != 'NEXRAD Radar'):
                time_string = self.dataSet.timeObj.strftime("%Y%m%d_%H%M%S")
                #Get variable name
                if (self.currentPType == 'Surface Statistics'):
                    if (self.selectedStatType == 'Bias'):
                        varname = 'MeanBias'
                    else:
                        varname = 'RMSE' 
                else:
                    if not self.derivedVar:
                        varname = self.dataSet.variableList[self.currentVar]
                    else:
                        varname = self.dataSet.dvarlist[self.currentdVar]

                #change the default plot name - for saving figure
                self.figure.canvas.get_default_filename = lambda: (varname + '_' + time_string + '.png')

            #Set up vslice variables 
            if (self.currentPType == 'Vertical Slice' or (self.currentPType == 'Difference Plot' and 
                self.prevPType == 'Vertical Slice')):
                #Create user control box 
                if self.vslicebox is None:
                    self.verticalSliceControl()
                #Get variables
                #CMAQ doesn't these variables so pass 
                if (self.dataSet.dsetname != 'CMAQ'):
                    self.u10, self.v10, w, press, height = self.dataSet.getVertVars() 
                    #Set second contour variable to w
                    var2 = w
                #Get variable dimensions
                dims = var.shape

                #Get xz orientation values
                if self.orientList[self.currentOrient] == 'xz':
                    diff = abs(self.dataSet.glats[self.currentGrid-1][:,0] - self.ref_pt)
                    ind = np.argsort(diff)
                    horiz = np.tile(np.sum(self.dataSet.glons[self.currentGrid-1][ind[0:4],:],axis=0)/4.,(dims[0],1))
                    horiz2 = np.tile(np.sum(self.dataSet.glats[self.currentGrid-1][ind[0:4],:],axis=0)/4.,(dims[0],1))
                    pltfld = np.squeeze(np.sum(var[:,ind[0:4],:],axis=1)/4.)
                    if (self.dataSet.dsetname != 'CMAQ'):
                        plevs = np.squeeze(np.sum(press[:,ind[0:4],:],axis=1)/4.)
                        if self.dataSet.dsetname != 'MET':
                            pvar2 = np.squeeze(np.sum(var2[:,ind[0:4],:],axis=1)/4.)
                            wwind = np.squeeze(np.sum(w[:,ind[0:4],:],axis=1)/4.)
                        hgt = np.squeeze(np.sum(height[:,ind[0:4],:],axis=1)/4.)
                        hwind = np.squeeze(np.sum(self.u10[:,ind[0:4],:],axis=1)/4.)
                    else:
                        plevs = np.swapaxes(np.tile(np.linspace(1,dims[0],dims[0]),(dims[2],1)),0,1)
                    
                    xdim = self.dataSet.nx[self.currentGrid-1]-1

                #Get yz orientation values
                if self.orientList[self.currentOrient] == 'yz':
                    diff = abs(self.dataSet.glons[self.currentGrid-1][0,:]-self.ref_pt)
                    ind = np.argsort(diff)
                    horiz = np.tile(np.sum(self.dataSet.glats[self.currentGrid-1][:,ind[0:4]],axis=1)/4.,(dims[0],1))
                    horiz2 =np.tile(np.sum(self.dataSet.glons[self.currentGrid-1][:,ind[0:4]],axis=1)/4.,(dims[0],1))
                    pltfld = np.squeeze(np.sum(var[:,:,ind[0:4]],axis=2)/4.)
                    if (self.dataSet.dsetname != 'CMAQ'):
                        plevs = np.squeeze(np.sum(press[:,:,ind[0:4]],axis=2)/4.)
                        if self.dataSet.dsetname != 'MET':
                            pvar2 = np.squeeze(np.sum(var2[:,:,ind[0:4]],axis=2)/4.)
                            wwind = np.squeeze(np.sum(w[:,:,ind[0:4]],axis=2)/4.)
                        hgt = np.squeeze(np.sum(height[:,:,ind[0:4]],axis=2)/4.)
                        hwind = np.squeeze(np.sum(self.v10[:,:,ind[0:4]],axis=2)/4.)
                    else:
                        plevs = np.swapaxes(np.tile(np.linspace(1,dims[0],dims[0]),(dims[1],1)),0,1)
                    
                    xdim = self.dataSet.ny[self.currentGrid-1]-1

            #Initialize colorbar range and increment
            if (self.colormin is None or self.colormax is None or self.ncontours is None):
                self.colormin = np.nanmin(pltfld)
                self.colormax = np.nanmax(pltfld)
                self.ncontours = 41.
                #Check for NAN
                if self.colormin != self.colormin:
                    if (self.currentPType == 'Difference Plot' or self.currentPType == 'Surface Statistics'):
                        self.colormin = -1.
                    else:
                        self.colormin = 0.
                if self.colormax != self.colormax:
                    self.colormax = 1.         
                #Make sure max and min are not equal       
                if self.colormax == self.colormin:
                    if self.colormax == 0:
                        self.colormax = 1.
                    else:
                        self.colormax = np.abs(self.colormax) * 2
            
            #For difference plot to normalize colorscale
            if (self.currentPType == 'Difference Plot' or (self.currentPType == 'Surface Statistics' 
                and self.selectedStatType == 'Bias')):
                if (self.colormin >= 0):
                    self.colormin = -1.*self.colormax
                if (self.colormax <= 0):
                    self.colormax = -1.*self.colormin                            
                #Colormap normalization
                midpoint = 1 - self.colormax/(self.colormax+abs(self.colormin))
                self.shiftedColorMap(cmap=matplotlib.cm.seismic, midpoint=midpoint, name='shifted')
                #self.cmap = self.newmap
                self.cmap = 'shifted'
            #Sets up user colorbar control   
            if self.colorbox is None:
                self.controlColorBar()
            else:
                #Change colorbox control values dynamically
                self.selectMin.setText(str(self.colormin))
                self.selectMax.setText(str(self.colormax))
                self.selectcontours.setText(str(self.ncontours))

            #Remove old items before creating new contour
            #Remove second control
            #if self.appobj.cs2 != None:
            #    for coll in self.appobj.cs2.collections:
            #        coll.remove()
            #    for labels in self.appobj.cs2label:
            #        labels.remove()
            if self.cs2 != None:
                for coll in self.cs2.collections:
                    coll.remove()
                for labels in self.cs2label:
                    labels.remove()
            #Remove wind barbs
            if self.barbs != None:
                self.barbs.remove()
                self.vectors2.remove()
            #Remove vectors
            if self.vectors != None:
                self.vectors.remove()
            #Remove vector magnitude legend
            if self.vectorkey != None:
                self.vectorkey.remove()
                
            #Define plotting levels
            lvls = np.linspace(self.colormin,self.colormax,self.ncontours)
            #Create plots - either contourf or pcolormesh
            if (self.filltype == "contourf" and self.currentPType != "Surface Statistics"):
                #Remove old contour before plotting - don't have to recall projection
                if (self.cs is not None):
                    for coll in self.cs.collections:
                        coll.remove()
                #Create contour on a map
                #Check if map exists to plot on
                #Second condition allows the creation of a vertical slice difference plot
                if (self.currentPType == 'Vertical Slice' or (self.currentPType == 'Difference Plot' and
                    self.prevPType == 'Vertical Slice')):
                    self.cs = self.axes1.contourf(
                                      horiz,
                                      plevs, pltfld,
                                      levels=lvls,
                                      cmap=plt.cm.get_cmap(str(self.cmap)),
                                      extend=self.extend)
                #if the dataset doesn't have a map, use the grid
                # created in the Dataset object. Note:
                # the x,y grid needs to be glons, glats here
                elif (self.dataSet.map[self.currentGrid-1] == None):
                    self.cs = self.axes1.contourf(
                                      self.dataSet.glons[self.currentGrid-1][::self.yinterval,::self.xinterval],
                                      self.dataSet.glats[self.currentGrid-1][::self.yinterval,::self.xinterval],
                                      pltfld[::self.yinterval,::self.xinterval],
                                      levels=lvls,
                                      cmap=plt.cm.get_cmap(str(self.cmap)),
                                      alpha=alpha,extend=self.extend)
                    #restrict the axes for visual appeal
                    #self.appobj.axes1[self.currentGrid-1].set_xlim(
                    #self.dataSet.glons[self.currentGrid-1].min(),
                    #self.dataSet.glons[self.currentGrid-1].max())
                    #restrict the axes for visual appeal
                    #self.appobj.axes1[self.currentGrid-1].set_ylim(
                    #self.dataSet.glats[self.currentGrid-1].min(),
                    #self.dataSet.glats[self.currentGrid-1].max())
                else:
                    self.cs = self.dataSet.map[self.currentGrid-1].contourf(
                                      self.dataSet.glons[self.currentGrid-1][::self.yinterval,::self.xinterval],
                                      self.dataSet.glats[self.currentGrid-1][::self.yinterval,::self.xinterval],
                                      pltfld[::self.yinterval,::self.xinterval],
                                      levels=lvls,
                                      latlon=True,extend=self.extend,
                                      cmap=plt.cm.get_cmap(str(self.cmap)),
                                      alpha=alpha, ax=self.axes1)

            elif (self.filltype == "pcolormesh" and self.currentPType != "Surface Statistics"):
                #Mask values less than the minimum for pcolormesh
                pltfld = np.ma.masked_less_equal(pltfld,self.colormin)
                #Normalize colors between min and max
                norm = matplotlib.colors.Normalize(vmin=np.amin(lvls),vmax=np.amax(lvls))
                #Clear figure if there are masked values - create new axis and set the projection
                if (np.ma.count_masked(pltfld).any()):
                    self.axes1 = None
                    self.domain_average = None
                    self.ColorBar = None
                    self.meridians = None
                    self.parallels = None
                    self.figure.clear()
                    self.axes1 = self.figure.add_subplot(111)
                    self.dataSet.resolution = self.resolution
                    if (self.dataSet.map[self.currentGrid-1] != None):
                        self.dataSet.map[self.currentGrid-1].ax = self.axes1
                
                #Create pcolormesh on a map
                #Check if map exists to plot on
                #Second condition allows the creation of a vertical slice difference plot
                if (self.currentPType == 'Vertical Slice' or (self.currentPType == 'Difference Plot' and
                    self.prevPType == 'Vertical Slice')):
                    self.cs = self.axes1.pcolormesh(
                                      horiz, 
                                      plevs, pltfld, 
                                      norm=norm,
                                      cmap=plt.cm.get_cmap(str(self.cmap)))
                #if the dataset doesn't have a map, use the grid
                # created in the Dataset object. Note:
                # the x,y grid needs to be glons, glats here
                elif (self.dataSet.map[self.currentGrid-1] == None):
                    self.cs = self.axes1.pcolormesh(
                                      self.dataSet.glons[self.currentGrid-1][::self.yinterval,::self.xinterval],
                                      self.dataSet.glats[self.currentGrid-1][::self.yinterval,::self.xinterval],
                                      pltfld[::self.yinterval,::self.xinterval],
                                      norm=norm,
                                      cmap=plt.cm.get_cmap(str(self.cmap)),
                                      alpha=alpha)
                    #restrict the axes for visual appeal
                    #self.appobj.axes1[self.currentGrid-1].set_xlim(
                    #self.dataSet.glons[self.currentGrid-1].min(),
                    #self.dataSet.glons[self.currentGrid-1].max())
                    #restrict the axes for visual appeal
                    #self.appobj.axes1[self.currentGrid-1].set_ylim(
                    #self.dataSet.glats[self.currentGrid-1].min(),
                    #self.dataSet.glats[self.currentGrid-1].max())
                else:
                    self.cs = self.dataSet.map[self.currentGrid-1].pcolormesh(
                                      self.dataSet.glons[self.currentGrid-1][::self.yinterval,::self.xinterval],
                                      self.dataSet.glats[self.currentGrid-1][::self.yinterval,::self.xinterval],
                                      pltfld[::self.yinterval,::self.xinterval],
                                      norm=norm,
                                      latlon=True,cmap=plt.cm.get_cmap(str(self.cmap)),
                                      alpha=alpha, ax=self.axes1)

            elif (self.currentPType == "Surface Statistics"):
                #Normalize colors between min and max
                norm = matplotlib.colors.Normalize(vmin=np.amin(lvls),vmax=np.amax(lvls))
                self.cs = self.dataSet.map[self.currentGrid-1].scatter(
                                      self.obslon,self.obslat,4*abs(pltfld),c=pltfld,
                                      norm=norm,cmap=plt.cm.get_cmap(str(self.cmap)),
                                      latlon=True)
                  
            #Clear old colorbar
            if self.ColorBar != None:
                self.ColorBar.remove()
            
            #Create colorbar
            #Second condition allows the creation of a vertical slice difference plot colorbar
            if (self.currentPType ==  'Vertical Slice' or (self.currentPType == 'Difference Plot' and
                self.prevPType == 'Vertical Slice')):
                self.ColorBar = self.figure.colorbar(self.cs, ax=self.axes1, 
                                                     orientation='horizontal', pad=0.1)  
            else:
                #Create colorbar axis    
                divider = make_axes_locatable(self.axes1)
                cax = divider.append_axes('right', size="5%", pad=0.1)
                self.ColorBar = self.figure.colorbar(self.cs,cax=cax, orientation='vertical')         

            #Extend colorbar if contourf - pcolormesh doesn't allow extend at this point
            self.ColorBar.extend = self.extend
           
            #Create colorbar ticks 
            self.ColorBar.ax.tick_params(labelsize=9)                
            self.axes1.set_title(self.varTitle,fontsize = 10)

            #Call function to determine values on map while you hoover your mouse
            if (self.currentPType == "Surface Statistics"):
                pass #Gives an error for this plot, haven't looked into the problem
            else:
                self.displayValuesXYZ(pltfld)
            
            #Display domain average on the map
            if self.domain_average != None:
                self.domain_average.remove()
            #Calculate domain average
            if (self.currentPType == "Surface Statistics"):
                davg = self.dom_stat
                if (self.selectedStatType == 'Bias'):
                    name = "Mean Bias: "
                elif (self.selectedStatType == 'Error'):
                    name = "RMSE: "
            else:
                davg = np.nanmean(pltfld)
                name = "Domain Average: "
            self.domain_average = self.axes1.text(0.95, -0.12,
                 (name + str(davg)),
                 verticalalignment='bottom',horizontalalignment='right',
                 transform = self.axes1.transAxes,
                 color='k',fontsize=10)

            #Plot wind barbs or vectors
            if ((self.plotbarbs == True or self.plotvectors == True) and self.dataSet.dsetname != 'CMAQ'):
                #Create map coordinates for wind barbs and vectors
                if (self.currentPType == 'Vertical Slice' and self.dataSet.dsetname != 'MET'):
                    xb = horiz
                    yb = plevs
                    self.u10 = hwind
                    self.v10 = wwind
                    ydim = self.dataSet.nz[self.currentGrid-1]-1
                    #Set wind/vector intervals
                    if (self.dataSet.nz[self.currentGrid-1] < 60.):
                        interval = 2
                    else:
                        interval = 4
                    if (self.dataSet.ny[self.currentGrid-1] < 250. and self.dataSet.nx[self.currentGrid-1] < 250.):
                        interval2 = 5
                    else:
                        interval2 = 10
                elif (self.currentPType != 'Vertical Slice'):                    
                    xb, yb = self.dataSet.map[self.currentGrid-1](
                             self.dataSet.glons[self.currentGrid-1],
                             self.dataSet.glats[self.currentGrid-1])
                                
                    #Get the winds for non derived variables
                    #Winds are already set for derived variables in DerivedVar.py
                    if not self.derivedVar:
                        #Read field to get winds
                        self.readField()
                        #Set winds to current level
                        if self.nz != 1:
                           self.u10 = self.dataSet.u10[self.currentLevel]
                           self.v10 = self.dataSet.v10[self.currentLevel]
                        else:
                           self.u10 = self.dataSet.u10
                           self.v10 = self.dataSet.v10

                    #Handle winds for MERRA and NCAR Reanalysis datasets 
                    if (self.dataSet.dsetname == "MERRA" or self.dataSet.dsetname == 'NCEP/NCAR Reanalysis II'):
                        if len(self.u10.shape) == 3:
                            self.u10 = self.dataSet.u10[self.currentLevel]
                        if len(self.v10.shape) == 3:
                            self.v10 = self.dataSet.v10[self.currentLevel]
                        #Set grid interval
                        if self.dataSet.dsetname == 'MERRA':
                            interval = 20
                        if self.dataSet.dsetname == 'NCEP/NCAR Reanalysis II':
                            interval = 8

                    #Set WRF wind barb/vector interval based on grid spacing
                    #These are just set based on what looked acceptable
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
                    interval2 = interval
                    xdim = self.dataSet.nx[self.currentGrid-1]-1
                    ydim = self.dataSet.ny[self.currentGrid-1]-1
                else:
                    #don't allow metfiles to try and plot vertical velocity
                    self.plotbarbs = False ; self.plotvectors = False
                #Plot wind barbs
                if (self.plotbarbs == True):
                    self.barbs = self.axes1.barbs(
                        xb[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-int(interval2/2.):interval2], 
                        yb[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-int(interval2/2.):interval2], 
                        self.u10[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-int(interval2/2.):interval2], 
                        self.v10[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-int(interval2/2.):interval2], 
                        length = 5, sizes={"emptybarb":0.0}, barbcolor='k', flagcolor='k',linewidth=0.50, pivot='middle')

                    #Create less than 5 m/s wind direction vectors 
                    unorm = self.u10/np.sqrt(np.power(self.u10,2) + np.power(self.v10,2))
                    vnorm = self.v10/np.sqrt(np.power(self.u10,2) + np.power(self.v10,2))
                    self.vectors2 = self.axes1.quiver(
                        xb[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-int(interval2/2.):interval2],
                        yb[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-int(interval2/2.):interval2],
                        unorm[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-1-int(interval2/2.):interval2],
                        vnorm[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-int(interval2/2.):interval2],
                        pivot='middle', headwidth=0, headlength=0, headaxislength=0,scale=60,width=0.0015)

                #Plot vectors
                if self.plotvectors == True:
                        #Calculate the max wind speed based on the selected grid points - used for scaling and legend
                        maxspeed = np.nanmax((self.u10[::interval,::interval2]**2+self.v10[::interval,::interval2]**2)**(1./2.))
                        self.appobj.vectors = self.axes1.quiver(
                        xb[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-int(interval2/2.):interval2],
                        yb[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-int(interval2/2.):interval2],
                        self.u10[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-1-int(interval2/2.):interval2],
                        self.v10[int(interval/2.):ydim-int(interval/2.):int(interval/2.),interval2:xdim-int(interval2/2.):interval2],                        
                        color='k',units='inches',pivot='mid',scale=maxspeed*4)
                        #Plot vector legend
                        self.vectorkey = self.axes1.quiverkey(self.appobj.vectors, 0.05, -0.08, 
                                                int(maxspeed*.75), str(int(maxspeed*.75))+' $m s^{-1}$',labelpos='E')

            #Second contour
            if (self.plotcontour2 == True and self.dataSet.dsetname != 'CMAQ'):
                if (self.currentPType == 'Vertical Slice'):
                    if (self.dataSet.dsetname != 'MET'):
                        self.cs2 = self.axes1.contour(horiz, plevs,
                                          pvar2, colors='k', linewidths=1.5)
                    else:
                        self.plotcontour2 = False
                else:
		              #Get second contour field
                    self.readField()
                    #Derived variables have their own 2nd contour options at this point
                    if not self.derivedVar:
                        self.var2 = pltfld
                    #Create contour
                    self.cs2 = self.dataSet.map[self.currentGrid-1].contour(
                                      self.dataSet.glons[self.currentGrid-1],
                                      self.dataSet.glats[self.currentGrid-1],
                                      self.var2,
                                      latlon=True,
                                      colors='k',
                                      linewidths=1.5,
                                      ax = self.axes1)
                     
                #Sets up displayed values
                if self.plotcontour2 == True:
                    if np.abs(self.var2).max() <= 1:
                        fmt = '%8.6f'
                    else:
                        fmt = '%d'
                    self.cs2label = self.cs2.clabel(inline=1, 
                                           font=10, fmt=fmt)

        #Create geography
        #Second condition on vertical slice allows changing of grid when Vertical Slice
        # selected but no variable selected
        if (self.currentPType == 'Difference Plot' and self.prevPType == 'Vertical Slice'):
           flag = True
        else:
           flag = False

        if ((self.currentPType != 'Vertical Slice' 
            or (self.currentVar == None and self.currentdVar == None))
            and self.dataSet.map[self.currentGrid-1] != None and not flag): 
            #print ("here", self.currentPType, self.prevPType)
            #print ("Map resolution:", self.dataSet.map[self.currentGrid-1].resolution)
            #print ("Current Grid:", self.currentGrid-1, self.currentGrid)
            #for i in range(len(self.dataSet.map)):
            #   print ("All Map resolution:", i, self.dataSet.map[i].resolution)
            #print ("\n Coasts:", self.coasts, "\n")
            #print ("\n Countries:", self.countries, "\n")
            #print ("\n Counties:", self.counties, "\n")
            #print ("\n States:", self.states, "\n")
            if self.plotcoasts:
                #print("Plot coasts:",self.dataSet.map[self.currentGrid-1].resolution)
                self.coasts = self.dataSet.map[self.currentGrid-1].drawcoastlines(ax=self.axes1)
            if self.plotcountries:
                self.countries = self.dataSet.map[self.currentGrid-1].drawcountries(ax=self.axes1)
            if self.plotstates:
                #print("plotting states")
                self.states = self.dataSet.map[self.currentGrid-1].drawstates(ax=self.axes1)
            if self.plotcounties:
                self.counties = self.dataSet.map[self.currentGrid-1].drawcounties(ax=self.axes1,linewidth=0.5)
            if self.plotcities and self.cities is False:
                for info,shape in zip(self.dataSet.map[self.currentGrid-1].cities_info,self.dataSet.map[self.currentGrid-1].cities):
                    if info['UATYP10'] is 'U':
                        poly = mplPath.Path(np.array(shape))
                        x,y = zip(*shape)
                        self.dataSet.map[self.currentGrid-1].plot(x,y, color='k', linewidth=0.5, ax=self.axes1)
                self.cities = True
            # draw parallels
            if self.dataSet.projectionType == "robin" or self.dataSet.projectionType == "geos" or \
               self.dataSet.projectionType == "ortho" or self.dataSet.projectionType == "aeqd":
                parallels = np.arange(-90.,90.,15)
                meridians = np.arange(0.,360., 30)
            else:
                #Create intervals than look resonable
                plim = max(int((np.nanmax(np.array(np.abs(self.dataSet.glats[self.currentGrid-1])))
                               -np.nanmin(np.array(np.abs(self.dataSet.glats[self.currentGrid-1]))))/10.),1) 
                parallels = np.arange(-90.,90.,plim)
                mlim = max(int((np.nanmax(np.array(np.abs(self.dataSet.glons[self.currentGrid-1])))
                               -np.nanmin(np.array(np.abs(self.dataSet.glons[self.currentGrid-1]))))/10.),1)
                meridians = np.arange(0.,360.,mlim)

            #Can't put labels on Geostationary, Orthographic or Azimuthal Equidistant basemaps 	
            if self.dataSet.projectionType == "ortho" and self.dataSet.projectionType == "aeqd" and self.dataSet.projectionType == "geos":
                if self.plotlatlon:
                    linewidth=0.5
                else:
                    linewidth=0
                # draw parallels
                self.dataSet.map[self.currentGrid-1].drawparallels(parallels,ax=self.axes1)
                # draw meridians
                self.dataSet.map[self.currentGrid-1].drawmeridians(meridians,ax=self.axes1)
            else:	
                if self.plotlatlon:
                    linewidth=0.5
                else:
                    linewidth=0
                # draw parallels
                if self.parallels == None:
                    self.parallels = self.dataSet.map[self.currentGrid-1].drawparallels(parallels,labels=[1,0,0,0],fontsize=8,
                                 ax=self.axes1,linewidth=linewidth)
                # draw meridians
                if self.meridians == None:
                    self.meridians = self.dataSet.map[self.currentGrid-1].drawmeridians(meridians,labels=[0,0,0,1],fontsize=8,
                                 ax=self.axes1,linewidth=linewidth)
            #Switch to not call projection again until grid or dataset is changed - for speed
            self.recallProjection = False
        elif ((self.currentPType == 'Vertical Slice' or (self.currentPType == 'Difference Plot' and 
               self.prevPType == 'Vertical Slice')) and (self.currentVar != None or
              self.currentdVar != None)):
            ax1 = self.axes1
            #Have to handle vertical slice of CMAQ differently because Pressure/height
            #  isn't in the standard output
            if (self.dataSet.dsetname == 'CMAQ'):
                ax1.set_ylabel("Model Level")
                #Set the x and y limits
                ax1.set_ylim(1,dims[0])
                ax1.set_xlim(np.amin(horiz),np.amax(horiz))
                ax1.set_title(self.varTitle,fontsize = 10)                
            else:
                #Set up map insert
                ax1.set_ylabel("Pressure [hPa]")
                ax1.set_yscale('log')
                ax1.set_xlim(horiz[0,:].min(),horiz[0,:].max())
                ax1.set_ylim(plevs.max(), self.minpress)
                subs = [2,3,5,7,8.5]
                loc2 = matplotlib.ticker.LogLocator(base=10., subs=subs)
                ax1.yaxis.set_minor_locator(loc2)
                ax1.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
                ax1.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
                ax1.set_title(self.varTitle,fontsize = 10)
                if self.dataSet.dsetname == "MERRA":
                    if self.dataSet.grid[5] == 'P':
                        min_ind = 0
                    else:
                        min_ind = -1
                else:
                    min_ind = 0
                ax1.fill_between(horiz[min_ind,:],1100., plevs[0,:], facecolor='peru')
                ax1.plot(horiz[min_ind,:],plevs[0,:],color='black')  
                ax2 = self.axes1.twinx()
                ax2.set_ylabel("Altitude [km]")
                diff = abs(plevs[:,0] - self.minpress)
                ind = np.argsort(diff)
                ax2.set_ylim(hgt.min(), hgt[ind[0],0])
                fmt = matplotlib.ticker.FormatStrFormatter("%g")
                ax2.yaxis.set_major_formatter(fmt)
            #Create map inset
            if (self.dataSet.map[self.currentGrid-1] != None):
                axins = inset_axes(ax1,width="30%",height="30%",
                                   loc=1,borderpad=0)
                map2 = self.dataSet.map[self.currentGrid-1]
                map2.ax = axins
                map2.etopo()
                map2.drawcoastlines()
                map2.drawcountries()
                map2.drawstates()
                # draw parallels.
                plim = max(int((abs(self.dataSet.glats[self.currentGrid-1]).max()-abs(self.dataSet.glats[self.currentGrid-1]).min())/2.5),0.5)
                parallels = np.arange(-90.,90.,plim)
                map2.drawparallels(parallels,labels=[0,1,0,0],fontsize=8)
                # draw meridians
                mlim = max(int((abs(self.dataSet.glons[self.currentGrid-1]).max()-abs(self.dataSet.glons[self.currentGrid-1]).min())/2.5),0.5)
                meridians = np.arange(0.,360.,mlim)
                map2.drawmeridians(meridians,labels=[0,0,1,0],fontsize=8)
                if self.orientList[self.currentOrient] == 'xz':
                    xx, yy = map2(horiz[0,:], horiz2[0,:])
                if self.orientList[self.currentOrient] == 'yz':
                    xx, yy = map2(horiz2[0,:],horiz[0,:])
                map2.plot(xx,yy,color='red',linewidth=2)
            #Have to delete the axis everytime
            self.recallProjection = True

        #Add radar location as a black dot - cone of silence
        if self.dataSet.dsetname == 'NEXRAD Radar':
            self.dataSet.map[self.currentGrid-1].plot(self.dataSet.lon0[0],
                    self.dataSet.lat0[0], marker='o',markersize=4,color='black',latlon=True, ax = self.axes1)

        #Plot sounding locations on the map and add a title
        if (self.dataSet.dsetname == 'Sounding'):
            #Plot station markers
            self.dataSet.map[self.currentGrid-1].plot(self.dataSet.glons[self.currentGrid-1],
                                                      self.dataSet.glats[self.currentGrid-1],
                                                      'bo',markersize=6,latlon=True)
            #Add Title
            soundingTitle = 'Click on a Sounding Location to Create a Plot'
            self.axes1.set_title(soundingTitle,fontsize = 12)
        #Keep the same colormap
        self.changeColor = False
        
        # Nair
        line, = self.axes1.plot([0], [0])  # empty line
        self.dataSelector = DataSelector(self,line)
        
        import time
        print( time.clock() - t0,"seconds process time plots")
        print( time.time() - t1,"seconds wall time plots")

    def plotSkewT(self,i,j):

        #Determine the Dataset Type     
        if (self.dataSet.dsetname == 'Sounding'):
            self.dataSet.getObsFile(ind = i,year = self.dataSet.year, month = self.dataSet.month,
                                    day = self.dataSet.day, hour = self.dataSet.hour)
            #Get above ground level height
            z = self.dataSet.h - self.dataSet.h[0]
            self.skewt = SkewT.SkewTobj(temp = self.dataSet.t, pressure = self.dataSet.p,
                                        Z = z, qv = self.dataSet.q,
                                        td = self.dataSet.dew, parcel = self.skewParcel,
                                        u = self.dataSet.u, v = self.dataSet.v )

        #The WRF MET files have different pressure, temp and humidity variables
        elif (self.dataSet.dsetname == 'MET'):
            #Get Variables
            p = np.squeeze(self.dataSet.readNCVariable('PRES'))/100. #Convert from pascal
            t = np.squeeze(self.dataSet.readNCVariable('TT'))-273.15 #Convert from K
            rh = np.squeeze(self.dataSet.readNCVariable('RH'))/100. #Convert from % 
            u = self.dataSet.readNCVariable('UU')
            v = self.dataSet.readNCVariable('VV')
           
            #Estimate height from pressure
            height = np.log(p[:,j,i]/p[0,j,i])*(-8000.)
             
 
            #Make sure RH values are valid (0-100%)
            rh[rh > 1.] = 1.
            rh[rh < 0.0] = 0.0
            #Calculate vapor pressure from temperature and RH
            qv = 6.11*np.exp((17.67*t)/(243.5+t))*rh
            self.skewt = SkewT.SkewTobj(temp = t[:,j,i], pressure = p[:,j,i], qv = qv[:,j,i], parcel = self.skewParcel, 
                                        Z = height,u = u[:,j,i], v = v[:,j,i])
                       

        else:
            #WRF
            #Get Variables 
            p = np.squeeze(self.dataSet.readNCVariable('P'))
            pb = np.squeeze(self.dataSet.readNCVariable('PB'))
            theta_prime = np.squeeze(self.dataSet.readNCVariable('T'))
            qv = np.squeeze(self.dataSet.readNCVariable('QVAPOR'))
            ph = self.dataSet.readNCVariable('PH')
            phb = self.dataSet.readNCVariable('PHB')
            u_in = self.dataSet.readNCVariable('U')
            v_in = self.dataSet.readNCVariable('V')
       
            #Construct full fields        
            pressure = (p + pb) * 0.01
            temp = (theta_prime + 300.) * (pressure/1000.) ** (287.04/1004.)
            height = (ph + phb)/9.81
            
            #Unstagger
            #get the height array in reference to ground level
            z = wrf.unstaggerZ(height) - wrf.unstaggerZ(height)[0,:,:]
            u = wrf.unstaggerX(u_in)
            v = wrf.unstaggerY(v_in)  

            self.skewt = SkewT.SkewTobj(temp = temp[:,j,i], pressure = pressure[:,j,i], 
                                        Z = z[:,j,i], qv = qv[:,j,i], parcel = self.skewParcel,
                                        u = u[:,j,i], v = v[:,j,i])

        #Add Sounding Station Name to Title
        if (self.dataSet.dsetname == 'Sounding'):
            self.skewt.fig.suptitle(self.dataSet.stat_id[i]+ ' : ' +self.dataSet.getTime())
        else: 
            self.skewt.fig.suptitle(self.dataSet.getTime())

        #Create Figure
        self.skewt.fig.canvas.set_window_title('Skew-T at i='+str(i)+' j='+str(j))
        self.skewt.ax1.set_title(self.Plist[self.selectedParcel]+' Profile')   
        self.skewt.ax2.set_title('Hodograph [kt]')     
        self.skewt.ax1.set_xlabel('Temperature [$^o$C]')# \n'+self.dataSet.getTime())
        

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
            subs = [0.5,2,3,4,5,6,7,8,9]
            loc = matplotlib.ticker.LogLocator(base=10., subs=subs)
            ax.yaxis.set_minor_locator(loc)
            ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
            ax.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
            #fmt = matplotlib.ticker.FormatStrFormatter("%g")
            #ax.yaxis.set_major_formatter(fmt)
            ax.margins(0.02)
            #Turn major and minor grid lines on
            ax.yaxis.grid(True, which='both')
        else:
           self.errorSelectVVar()

    def plotTimeSeries(self):

        #Check if the time indices were changed
        sind = self.startTime.currentIndex()
        eind = self.endTime.currentIndex()
        #If the time indices changed, get the variable over the new
        # time range
        if (sind != self.timeSeriesStart or eind != self.timeSeriesEnd):
            self.createTimeSeries(self.timeI,self.timeJ)
        
        #Get start and end time
        start = self.time_series[0].strftime("%b %d %Y, %H:%M:%S")
        end = self.time_series[-1].strftime("%b %d %Y, %H:%M:%S UTC")
        #If plotting a derived variable, use the set short title
        if (not self.derivedVar):
            varTitle = self.dataSet.description +' ('+self.dataSet.units+')'
        else:
            varTitle = self.sTitle
        varTitle = varTitle + '\n' + start + '  -  ' + end
        varTitle = varTitle + '\n i = '+ str(self.col) + ', j = ' + str(self.row) 
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
            varTitle = varTitle + '\n Level=' + str(self.dataSet.levelList[self.currentLevel])
        else:
            #2D Variable
            ax.plot_date(ind,self.var_series,'bo-')
        
        ax.set_title(varTitle)        

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
        QApplication.setOverrideCursor(Qt.WaitCursor) 
        #Save the indices so function can be called within the plotobj
        self.timeI = i
        self.timeJ = j
        #Hold the current time for resetting the data
        input_time_index = self.dataSet.currentTimeIndex
        #Start time
        self.timeSeriesStart = self.startTime.currentIndex()
        #End Time
        self.timeSeriesEnd = self.endTime.currentIndex()
        #Get variable over time range
        varSeries, self.time_series = self.getVarOverRange(self.timeSeriesStart,self.timeSeriesEnd)
        varSeries = np.array(varSeries)
        #If the number of dimensions is less than 3 then the selected variable is 2D
        if (len(self.var.shape) < 3):
            self.var_series = varSeries[:,j,i]
        else:
            self.var_series = varSeries[:,:,j,i]
        self.dataSet.setTimeIndex(input_time_index)
        QApplication.restoreOverrideCursor()

    #This function calculates the spatial statistics (mean,max,min) over a specified time interval
    def spatialStats(self):
        QApplication.setOverrideCursor(Qt.WaitCursor)
        #Save the dataset current time index
        input_time_index = self.dataSet.currentTimeIndex
        #Start time
        timeStatStart = self.spatialSTime.currentIndex()
        #self.timeSeriesEnd = self.endTime.getCurrentIndex()
        #End Time
        timeStatEnd = self.spatialETime.currentIndex()
        #Get variable over time range
        if (self.calTimeRange == True):
            self.statVar, self.timeRange = self.getVarOverRange(timeStatStart,timeStatEnd)        
            self.statVar = np.array(self.statVar)
            #Reset the data
            self.dataSet.setTimeIndex(input_time_index)
        else:
            #Reset
            self.calTimeRange = True
        #Create the plot title - different than the one created in readField
        #If plotting a derived variable, used the short title created in readField
        start = self.timeRange[0].strftime("%b %d %Y, %H:%M:%S")
        end = self.timeRange[-1].strftime("%b %d %Y, %H:%M:%S UTC")
        if (self.derivedVar == True):
            self.varTitle = self.spatialOptions[self.spatialSType.currentIndex()] +' '
            self.varTitle = self.varTitle + self.sTitle
        else:
            self.varTitle = self.spatialOptions[self.spatialSType.currentIndex()] +' '
            self.varTitle = self.varTitle + self.dataSet.description +' ('+self.dataSet.units+')'
        self.varTitle = self.varTitle + '\n' + start + '  -  ' + end     
        #Add the level to the variable title if 3D
        if (self.nz > 1):
            self.varTitle = self.varTitle + '\n Level = ' + str(self.currentLevel+1)
        QApplication.restoreOverrideCursor()

    #This function returns an array containing a variable over a given time range
    def getVarOverRange(self,start,end):
        #Create list to hold values at selected point and time
        var_sum = []
        timeRange = []
        append_var = var_sum.append
        append_time = timeRange.append
        #If the number of dimensions is less than 3 then the selected variable is 2D
        if (len(self.var.shape) < 3):
             for i in range(start,end+1):
                timeIndex = self.dataSet.setTimeIndex(i)
                self.dataSet.getTime()
                append_time(self.dataSet.timeObj)#.strftime("%Y-%m-%d %H:%M:%S"))
                self.readField()
                append_var(self.var)
        else:
             for i in range(start,end+1):
                timeIndex = self.dataSet.setTimeIndex(i)
                self.dataSet.getTime()
                append_time(self.dataSet.timeObj)
                self.readField()
                append_var(self.var)

        return var_sum, timeRange


    #This function is used to calculate anomalies
    def calculateAnomaly(self):
        #Check if the mean needs to be calculated again
        if (self.calTimeRange == True):
            #Save the dataset current time index
            input_time_index = self.dataSet.currentTimeIndex
            #Get the current selected hour,minute,second - used to get the simulation mean at that time
            cHour = self.dataSet.timeObj.strftime("%H:%M")
            #Determine number of times in loaded files 
            numfiles = self.dataSet.getNumFiles()
            ntimes = self.dataSet.ntimes
            #Loop through the files times and save times matching the current time of day
            meanVar = []
            for ii in range(numfiles*ntimes):
                self.dataSet.setTimeIndex(ii)
                #Update the time object
                self.dataSet.getTime()
                nHour = self.dataSet.timeObj.strftime("%H:%M")
                if (nHour == cHour):
                    #Get the current variable
                    self.readField()
                    #Check if the field has a z-component
                    if (self.nz == 1):
                        meanVar.append(self.var)
                    else:
                        meanVar.append(self.var[self.currentLevel])
            #Get the original data back
            self.dataSet.setTimeIndex(input_time_index)
            self.readField()
            #Check if there is at least more than 1 file occurring at the same time of day
            if (len(meanVar) >= 2):
                self.meanAnomaly = np.nanmean(meanVar,axis=0)
                self.cAnomaly = True
                self.varTitle = 'Anomalous ' + self.varTitle
            else:
                self.errorAnomaly()
                #Switch back to horizontal plot
                self.currentPType = 'Horizontal Slice'
                ind = np.where(np.array(self.dataSet.ptypes) == self.currentPType)[0]
                self.selectPlotType.setCurrentIndex(ind)
                self.cAnomaly = False
                self.pltFxn(self.pNum)       
        else:
            self.cAnomaly = True
            self.calTimeRange = True 

    #This function extracts the loaded observation data at the current time
    # and calculates the statitsics (Bias and RMSE)
    def sfcStats(self):
        #Check if time was changed
        if (self.statsChangeTime):
            #Get the index where the obs times match the current time
            ctime = self.dataSet.getTime(obj=True)
            ind = np.where(self.obsdate == ctime)[0]
        
            #Read in the model variables
            u10 = self.dataSet.readNCVariable('U10')
            v10 = self.dataSet.readNCVariable('V10')
            sinalpha = self.dataSet.readNCVariable('SINALPHA')
            cosalpha = self.dataSet.readNCVariable('COSALPHA')
            t2 = self.dataSet.readNCVariable('T2')
            q2 = self.dataSet.readNCVariable('Q2')*1000. #Convert to g/kg

            #Calculate Earth relative model winds and convert to wind speed
            # and meteorological wind direction in degrees
            uearth = u10*cosalpha - v10*sinalpha
            vearth = v10*cosalpha + u10*sinalpha
            wrf_ws = np.sqrt(uearth*uearth + vearth*vearth)
            wrf_wd = 270 - np.arctan2(vearth,uearth)*(180./np.pi) 

            #Get WRF lower left corner
            x0wrf, y0wrf = self.dataSet.map[self.currentGrid-1](
                           self.dataSet.glons[self.currentGrid-1][0,0],
                           self.dataSet.glats[self.currentGrid-1][0,0])
            x0wrf = x0wrf - self.dataSet.dx[self.currentGrid-1]/2.
            y0wrf = y0wrf - self.dataSet.dx[self.currentGrid-1]/2.

            #Determine the i and j indexs on the WRF grid cooresponding to the
            # obs locations and perform bilinear interpolation of the model to 
            # the obs point
            self.obslon = []
            self.obslat = []
            ws_bias = []
            wd_bias = []
            t_bias = []
            q_bias = []
            ws_err = []
            wd_err = []
            t_err = []
            q_err = []
            for i in range(len(ind)):
                x, y = self.dataSet.map[self.currentGrid-1](self.olon[ind][i],self.olat[ind][i])
                ii = (x-x0wrf)/self.dataSet.dx[self.currentGrid-1]
                jj = (y-y0wrf)/self.dataSet.dx[self.currentGrid-1]
                #Check if the observation falls within the grid
                if (np.floor(ii) > 0 and np.floor(ii) < self.dataSet.nx[self.currentGrid-1]-2 and
                    np.floor(jj) > 0 and np.floor(jj) < self.dataSet.ny[self.currentGrid-1]-2):

                    #Observations
                    ws = self.ws_obs[ind][i]
                    wd = self.wd_obs[ind][i]
                    t = self.t_obs[ind][i]
                    q = self.q_obs[ind][i]
                    self.obslon.append(self.olon[ind][i])
                    self.obslat.append(self.olat[ind][i])

                    #Bilinear interpolation to observation point
                    #Wind Speed
                    term1 = (1./((np.floor(ii)+1 - np.floor(ii))*(np.floor(jj)+1 - np.floor(jj))))
                    term2 = wrf_ws[np.int(jj),np.int(ii)]*(np.floor(ii)+1-ii)*(np.floor(jj)+1-jj)
                    term3 = wrf_ws[np.int(jj),np.int(ii)+1]*(ii-np.floor(ii))*(np.floor(jj)+1-jj)
                    term4 = wrf_ws[np.int(jj)+1,np.int(ii)]*(np.floor(ii)+1-ii)*(jj-np.floor(jj))
                    term5 = wrf_ws[np.int(jj)+1,np.int(ii)+1]*(ii-np.floor(ii))*(jj-np.floor(jj))
                    ws_mod = term1*(term2+term3+term4+term5)

                    #Wind Direction
                    term2 = wrf_wd[np.int(jj),np.int(ii)]*(np.floor(ii)+1-ii)*(np.floor(jj)+1-jj)
                    term3 = wrf_wd[np.int(jj),np.int(ii)+1]*(ii-np.floor(ii))*(np.floor(jj)+1-jj)
                    term4 = wrf_wd[np.int(jj)+1,np.int(ii)]*(np.floor(ii)+1-ii)*(jj-np.floor(jj))
                    term5 = wrf_wd[np.int(jj)+1,np.int(ii)+1]*(ii-np.floor(ii))*(jj-np.floor(jj))
                    wd_mod = term1*(term2+term3+term4+term5)

                    #Temperature
                    term2 = t2[np.int(jj),np.int(ii)]*(np.floor(ii)+1-ii)*(np.floor(jj)+1-jj)
                    term3 = t2[np.int(jj),np.int(ii)+1]*(ii-np.floor(ii))*(np.floor(jj)+1-jj)
                    term4 = t2[np.int(jj)+1,np.int(ii)]*(np.floor(ii)+1-ii)*(jj-np.floor(jj))
                    term5 = t2[np.int(jj)+1,np.int(ii)+1]*(ii-np.floor(ii))*(jj-np.floor(jj))
                    t_mod = term1*(term2+term3+term4+term5)

                    #Water Vapor Mixing Ratio
                    term2 = q2[np.int(jj),np.int(ii)]*(np.floor(ii)+1-ii)*(np.floor(jj)+1-jj)
                    term3 = q2[np.int(jj),np.int(ii)+1]*(ii-np.floor(ii))*(np.floor(jj)+1-jj)
                    term4 = q2[np.int(jj)+1,np.int(ii)]*(np.floor(ii)+1-ii)*(jj-np.floor(jj))
                    term5 = q2[np.int(jj)+1,np.int(ii)+1]*(ii-np.floor(ii))*(jj-np.floor(jj))
                    q_mod = term1*(term2+term3+term4+term5)

                    #Calculate bias
                    ws_bias.append(ws_mod - ws)
                    wd_bias.append(wd_mod - wd)
                    t_bias.append(t_mod - t)
                    q_bias.append(q_mod - q)

                    #Calculate Error
                    ws_err.append((ws_mod - ws)**2)
                    wd_err.append((wd_mod - wd)**2)
                    t_err.append((t_mod - t)**2)
                    q_err.append((q_mod - q)**2)

            self.ws_bias = np.array(ws_bias)
            self.ws_err = np.array(ws_err)
            self.wd_bias = np.array(wd_bias)
            self.wd_err = np.array(wd_err)
            self.t_bias = np.array(t_bias)
            self.t_err = np.array(t_err)
            self.q_bias = np.array(q_bias)
            self.q_err = np.array(q_err)        
        
        #Determine selected variable
        if (self.selectedStatsVar == 'Wind Speed'):
            bias_var = self.ws_bias
            mean_sq = self.ws_err
            err_var = np.sqrt(self.ws_err)
            units = 'm s$^{-1}$'
        elif (self.selectedStatsVar == 'Wind Direction'):
            bias_var = self.wd_bias
            mean_sq = self.wd_err
            err_var = np.sqrt(self.wd_err)
            units = 'degrees'
        elif (self.selectedStatsVar == 'Temperature'):
            bias_var = self.t_bias
            mean_sq = self.t_err
            err_var = np.sqrt(self.t_err)
            units = 'K'
        elif (self.selectedStatsVar == 'Water Vapor Mixing Ratio'):
            bias_var = self.q_bias
            mean_sq = self.q_err
            err_var = np.sqrt(self.q_err)
            units = 'g kg$^{-1}$'

        #Calculate mean bias
        mb = np.nanmean(bias_var)
        #Calculate RMSE
        rmse = np.sqrt(np.nanmean(mean_sq))
   
        #Determine selected stat type
        if (self.selectedStatType == 'Bias'):
           self.stat = bias_var
           self.dom_stat = mb
           self.varTitle = 'Model ' + self.selectedStatsVar + ' Bias: '
           self.varTitle = self.varTitle + '('+ units + ') \n'+self.dataSet.getTime()
        elif (self.selectedStatType == 'Error'):
           self.stat = err_var
           self.dom_stat = rmse
           self.varTitle = 'Model ' + self.selectedStatsVar + ' Error: '
           self.varTitle = self.varTitle + '('+ units + ') \n'+self.dataSet.getTime()
        
        #Reset 
        self.statsChangeTime = True

        #Recall projection
        self.recallProjection = True
        self.clear = True

    def displayValuesXYZ(self,var):
        self.axes1.imshow(var, interpolation='nearest')
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
                if self.dataSet.runType == 'IDEAL':
                    row = int(y*1000/self.dataSet.dy[self.currentGrid-1] + (self.dataSet.ny[self.currentGrid-1]-1)/2)
                    col = int(x*1000/self.dataSet.dx[self.currentGrid-1] + (self.dataSet.nx[self.currentGrid-1]-1)/2)
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

            #if (self.currentPType == 'SkewT/Hodograph' or 
            #    self.currentPType == 'Vertical Profile' or 
            #    self.currentPType == 'Time Series'):
            #   return 'i=%6.2f, j=%6.2f, value=%10.6f'%(col, row, z)
            #else:
            return 'lon=%6.2f, lat=%6.2f, i=%5i, j=%5i, value=%10.6f'%(x1, y1, col, row, z)
        if self.dataSet.dsetname != 'NEXRAD Radar':
            self.axes1.format_coord = format_coord

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

    #Initialization function for the plots
    # Needed for datasets that don't have specific variables
    # Also for changing the plot type before selecting a variable
    def initializePlotOptions(self):

        #Skew-T
        if (self.currentPType == 'SkewT/Hodograph'):

            self.replot2d = True

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

        #Vertical Profiles
        if (self.currentPType == 'Vertical Profile'):
            self.replot2d = True

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

            #Derived vertical variable control
            selectdVarWidget = QWidget()
            selectdVarWidgetLayout = QHBoxLayout()
            selectdVarWidget.setLayout(selectdVarWidgetLayout)

            selectdVarLabel = QLabel()
            selectdVarLabel.setText('Derived Vertical Variable:')

            selectdVVar = QComboBox()
            selectdVVar.setStyleSheet(Layout.QComboBox())
            self.dVvarlist = ['wind-3d']
            selectdVVar.addItems(self.dVvarlist)
            selectdVVar.setSizeAdjustPolicy(QComboBox.AdjustToContents)
            selectdVVar.activated.connect(self.selectionChangedVerticalVar)
            selectdVarWidgetLayout.addWidget(selectdVarLabel)
            selectdVarWidgetLayout.addWidget(selectdVVar)
            self.vertControlLayout.addWidget(selectdVarWidget)


            self.optionTabs.addTab(self.vertControl,vertTitle)
            self.optionTabs.setCurrentIndex(self.optionTabs.count()-1)
            self.tabbingLayout.addWidget(self.optionTabs)

        #Surface Statistics
        if (self.currentPType == 'Surface Statistics'):
            #Define intial variable and statistic type
            self.selectedStatsVar = 'Wind Speed'
            self.selectedStatType = 'Bias'
            #Get start and enddate of WRF files
            begdate = self.dataSet.timeList[self.dataSet.currentGrid][0]
            enddate = self.dataSet.timeList[self.dataSet.currentGrid][-1]
            #Read in the observation file
            obsFile = QFileDialog.getOpenFileName(self.appobj, 'Select File')
            if (len(obsFile) > 0):
                QApplication.setOverrideCursor(Qt.WaitCursor)
                self.obsdate, self.olat, self.olon, \
                              self.t_obs, self.q_obs, \
                              self.ws_obs, self.wd_obs = readObs.readObs(obsFile[0],begdate,enddate)
                QApplication.restoreOverrideCursor()
            else:
                print("File was not selected")

            #Create tab to control sfc stats variable
            self.statsControl = QGroupBox()
            statsTitle = 'Surface Statistic Control'
            self.statsControlLayout = QVBoxLayout(self.statsControl)

            selectVarWidget = QWidget()
            selectVarWidgetLayout = QVBoxLayout()
            selectVarWidget.setLayout(selectVarWidgetLayout)

            #Variable control
            selectVarLabel = QLabel()
            selectVarLabel.setText('Variable:')
            selectSfcVar = QComboBox()
            selectSfcVar.setStyleSheet(Layout.QComboBox())
            self.sfcList = ['Wind Speed', 'Temperature', 'Water Vapor Mixing Ratio']
            selectSfcVar.addItems(self.sfcList)
            selectSfcVar.setSizeAdjustPolicy(QComboBox.AdjustToContents)
            selectSfcVar.currentIndexChanged.connect(self.selectionChangedStatVar)
        
            #Stat type control
            selectStatTLabel = QLabel()
            selectStatTLabel.setText('Statistic:')
            selectStatVar = QComboBox()
            selectStatVar.setStyleSheet(Layout.QComboBox())
            self.statList = ['Bias', 'Error']
            selectStatVar.addItems(self.statList)
            selectStatVar.setSizeAdjustPolicy(QComboBox.AdjustToContents)
            selectStatVar.currentIndexChanged.connect(self.selectionChangedStatType)        
        
            #Add widgets
            selectVarWidgetLayout.addWidget(selectVarLabel)
            selectVarWidgetLayout.addWidget(selectSfcVar)
            selectVarWidgetLayout.addWidget(selectStatTLabel)
            selectVarWidgetLayout.addWidget(selectStatVar)
            self.statsControlLayout.addWidget(selectVarWidget)

            self.optionTabs.addTab(self.statsControl,statsTitle)
            self.optionTabs.setCurrentIndex(self.optionTabs.count()-1)
            self.tabbingLayout.addWidget(self.optionTabs)            

        #Time Series
        if (self.currentPType == 'Time Series'):
            self.replot2d = True

            #Create tab to control the time series time options
            self.seriesControl = QGroupBox()
            sTitle = 'Time Series Plot Control'
            self.seriesControlLayout = QVBoxLayout(self.seriesControl)

            #Start time control
            startWidget = QWidget()
            startWidgetLayout = QHBoxLayout(startWidget)

            startLabel = QLabel()
            startLabel.setText('Start Time:')

            self.startTime = QComboBox()
            self.startTime.setStyleSheet(Layout.QComboBox())
            self.startTime.addItems(self.dataSet.timeList[self.currentGrid])
            self.startTime.setCurrentIndex(0)
            self.startTime.setSizeAdjustPolicy(QComboBox.AdjustToContents)

            #End time control
            endWidget = QWidget()
            endWidgetLayout = QHBoxLayout(endWidget)

            endLabel = QLabel()
            endLabel.setText('End Time:')

            self.endTime = QComboBox()
            self.endTime.setStyleSheet(Layout.QComboBox())
            self.endTime.addItems(self.dataSet.timeList[self.currentGrid])
            self.endTime.setCurrentIndex(len(self.dataSet.timeList[self.currentGrid])-1)
            self.endTime.setSizeAdjustPolicy(QComboBox.AdjustToContents)

            #Plot button
            SeriesPlotButton = QPushButton('Replot')
            SeriesPlotButton.setStyleSheet(Layout.QPushButton3())
            SeriesPlotButton.resize(SeriesPlotButton.minimumSizeHint())
            SeriesPlotButton.clicked.connect(self.getTimeSeriesStartEnd)

            #Add widgets and create Tab
            startWidgetLayout.addWidget(startLabel)
            startWidgetLayout.addWidget(self.startTime)
            endWidgetLayout.addWidget(endLabel)
            endWidgetLayout.addWidget(self.endTime)
            self.seriesControlLayout.addWidget(startWidget)
            self.seriesControlLayout.addWidget(endWidget)
            self.seriesControlLayout.addWidget(SeriesPlotButton)
            self.optionTabs.addTab(self.seriesControl,sTitle)
            self.optionTabs.setCurrentIndex(self.optionTabs.count()-1)
            self.tabbingLayout.addWidget(self.optionTabs)

            #Popup to select the initial start and end time
            #self.tsControl = QDialog(self.appobj)
            #self.tsControl.resize(self.tsControl.minimumSizeHint())
            #self.tsControl.setWindowTitle("Time Series Start and End Time")
            #tsControlLayout = QVBoxLayout(self.tsControl)

            #Add widgets to the popup
            #tsControlLayout.addWidget(self.startWidget)
            #tsControlLayout.addWidget(self.endWidget)

            #self.tsControl.show()
        
        #Difference plot
        if (self.currentPType == 'Difference Plot' or (self.currentPType == 'Difference Plot' and
            self.prevPType == 'Vertical Slice')):
            #Reset colorbar settings
            self.colormax = None
            self.colormin = None
            self.ncontours = None
            self.cmap2 = self.cmap
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

        #Spatial Statistics
        if (self.currentPType == 'Spatial Stats'):
            #Initialize to True when selected
            self.calTimeRange = True

            #Create tab to control the spatial stats time options and condition
            self.spatialControl = QGroupBox()
            Title = 'Spatial Statistics Plot Control'
            self.spatialControlLayout = QVBoxLayout(self.spatialControl)

            #Statistics control
            statsWidget = QWidget()
            statsWidgetLayout = QHBoxLayout(statsWidget)
            statsLabel = QLabel()
            statsLabel.setText('Statistic Type:')

            self.spatialSType = QComboBox()
            self.spatialSType.setStyleSheet(Layout.QComboBox())
            self.spatialSType.addItems(self.spatialOptions)
            self.spatialSType.setCurrentIndex(0)
            self.spatialSType.setSizeAdjustPolicy(QComboBox.AdjustToContents)

            #Start time control
            startWidget = QWidget()
            startWidgetLayout = QHBoxLayout(startWidget)

            startLabel = QLabel()
            startLabel.setText('Start Time:')

            self.spatialSTime = QComboBox()
            self.spatialSTime.setStyleSheet(Layout.QComboBox())
            self.spatialSTime.addItems(self.dataSet.timeList[self.currentGrid])
            self.spatialSTime.setCurrentIndex(0)
            self.spatialSTime.setSizeAdjustPolicy(QComboBox.AdjustToContents)

            #End time control
            endWidget = QWidget()
            endWidgetLayout = QHBoxLayout(endWidget)

            endLabel = QLabel()
            endLabel.setText('End Time:')

            self.spatialETime = QComboBox()
            self.spatialETime.setStyleSheet(Layout.QComboBox())
            self.spatialETime.addItems(self.dataSet.timeList[self.currentGrid])
            self.spatialETime.setCurrentIndex(len(self.dataSet.timeList[self.currentGrid])-1)
            self.spatialETime.setSizeAdjustPolicy(QComboBox.AdjustToContents)

            #Plot button
            SpatialPlotButton = QPushButton('Create Plot')
            SpatialPlotButton.setStyleSheet(Layout.QPushButton3())
            SpatialPlotButton.resize(SpatialPlotButton.minimumSizeHint())
            SpatialPlotButton.clicked.connect(self.getTimeStatsStartEnd)

            #Add widgets and create Tab
            statsWidgetLayout.addWidget(statsLabel)
            statsWidgetLayout.addWidget(self.spatialSType)            
            startWidgetLayout.addWidget(startLabel)
            startWidgetLayout.addWidget(self.spatialSTime)
            endWidgetLayout.addWidget(endLabel)
            endWidgetLayout.addWidget(self.spatialETime)
            self.spatialControlLayout.addWidget(statsWidget)
            self.spatialControlLayout.addWidget(startWidget)
            self.spatialControlLayout.addWidget(endWidget)
            self.spatialControlLayout.addWidget(SpatialPlotButton)
            self.optionTabs.addTab(self.spatialControl,Title)
            self.optionTabs.setCurrentIndex(self.optionTabs.count()-1)
            self.tabbingLayout.addWidget(self.optionTabs)

        #Anomaly
        if (self.currentPType == 'Anomaly'):
            #Initialize to True when selected
            self.calTimeRange = True

    #Function that handles the plot type changes
    def selectionChangePlot(self,i):
        
        #Save previous plot type
        self.prevPType = self.currentPType

        #Reset clicked point values
        if self.cid is not None:
            self.cid = None
            self.row = None
            self.col = None

        #If changed from difference plot - get the old colormap back
        if (self.currentPType == 'Difference Plot' or self.currentPType == 'Surface Statistics'):
            self.cmap = self.cmap2
            
        #Set the new plot type
        self.currentPType = self.dataSet.ptypes[i]
		           
        #Remove all added items
        #Necessary to prevent error in plot after
        # clearing the figure below
        self.cs = None
        self.cs2 = None
        self.cs2label = None
        self.barbs = None
        self.vectors = None
        self.vectorkey = None
        self.domain_average = None
        #Get the new projection
        self.recallProjection = True
        self.coasts = None
        self.countries = None
        self.states = None
        self.counties = None
        self.parallels = None
        self.meridians = None
        #Removes the colorbar
        self.ColorBar = None
        #Reset the color scale unless plotting WRF reflectivity
        #Not all datasets have a dvarlist so only check if WRF at the moment
        #Revist this as additional datasets are added (if they have derived vars)
        if (self.dataSet.dsetname == 'WRF'):
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
        
        #Remove the plotting axes and clear the figure
        self.axes1 = None
        self.figure.clear()

        #Make sure to replot cities if necessary
        if (self.currentPType != 'Vertical Slice'):
            if (self.plotcities == True):
                self.cities = False

        #Remove vertical slice controls if necessary
        if self.vslicebox != None:
            self.vslicebox.setParent(None)
            self.vslicebox = None
 
        #Destroy vertical plot control widget if necessary 
        if self.vertControl != None:
            self.vertControl.setParent(None)
            self.vertControl = None
                
        #Destroy skew-T parcel control widget if necessary
        if self.pControl != None:
            self.pControl.setParent(None)
            self.pControl = None

        #Clear difference plot control widget if necessary
        if self.diffControl != None:
            self.diffControl.setParent(None)
            self.diffControl = None     

        #Clear time series plot control widget if necessary
        if self.seriesControl != None:
            self.seriesControl.setParent(None)
            self.seriesControl = None

        #Clear Spatial stats plot control widget if necessary
        if self.spatialControl != None:
            self.spatialControl.setParent(None)
            self.spatialControl = None

        #Clear Sfc Stats plot control widget if necessary
        if self.statsControl != None:
            self.statsControl.setParent(None)
            self.statsControl = None        
                 
        #Initialize the plot options    
        self.initializePlotOptions()    

        #Don't allow plotting if no variable is selected
        if (self.currentVar != None or self.currentdVar != None 
            or self.currentPType == 'Surface Statistics'):
            #Updata data here to accomidate multiple grids with multiple plots
            self.dataSet.setGrid(self.currentGrid)
            #Make sure there is a 3D variable selected before allowing 
            # a vertical slice plot to be made
            if (self.currentPType == 'Vertical Slice' and len(self.var.shape) < 3):
                #Switch back to the previous plot
                self.currentPType = self.prevPType
                #Determine previous plot type index
                ind = np.where(np.array(self.dataSet.ptypes) == self.currentPType)[0]
                self.selectPlotType.setCurrentIndex(ind)
                #Throw an error
                self.error3DVar()   
            else:            
                self.pltFxn(self.pNum) 

    #Changes the variable
    def selectionChangeVar(self,i):
        self.currentVar = i
        self.derivedVar = False
        self.readField()
        #initialize an index for derived variable
        # not used but needed in setPlotVars
        #self.currentdVar = 0
        self.setPlotVars()

    #Changes the derived variable
    def selectionChangedVar(self,i):
        self.currentdVar = i
        self.derivedVar = True
        #initialize an index for derived variable
        # not used but needed in setPlotVars
        self.currentVar = 0
        #Check if selection is simulated brightness temp
        #Wait until selections are made
        if (self.dataSet.dvarlist[self.currentdVar] == 'BrightTemp/Radiance'):
            #Create an initial sensor and channel selection box
            self.crtmControl = QDialog(self.appobj)
            self.crtmControl.resize(self.crtmControl.minimumSizeHint())
            self.crtmControl.setMinimumHeight(self.appobj.screeny*.25)
            self.crtmControl.setMinimumWidth(self.appobj.screenx*.25)
            self.crtmControl.setWindowTitle("Select sensor and channel")
            crtmControlLayout = QVBoxLayout(self.crtmControl)
            
            #Create sensor combo box
            sensorBoxLabel = QLabel()
            sensorBoxLabel.setText('Sensor Type: \n (Note: Visible channel calculations \n take longer than IR)')
            self.sensorBox = QComboBox()
            self.sensorBox.setStyleSheet(Layout.QComboBox())
            #Define sensor types
            self.sensor_names = ['GOES-16 ABI (IR)', 'GOES-16 ABI (Vis)',
                                 'GOES-15 Imager (IR)', 'GOES-15 Imager (Vis)',
                                 'GOES-14 Imager (IR)', 'GOES-14 Imager (Vis)',
                                 'GOES-13 Imager (IR)', 'GOES-13 Imager (Vis)',
                                 'GOES-12 Imager (IR)', 'GOES-12 Imager (Vis)'] 
            self.sensorBox.addItems(self.sensor_names)
            self.sensorBox.setSizeAdjustPolicy(QComboBox.AdjustToContents)
            #Connect to function to change the channel list
            self.sensorBox.currentIndexChanged.connect(self.getNewSensorInfo)

            #Create Variable Selection combo box
            # Brightness Temperature or Radiance
            varBoxLabel = QLabel()
            varBoxLabel.setText('Output Variable:')
            self.varBox = QComboBox()
            self.varBox.setStyleSheet(Layout.QComboBox())
            #Define Variable Options
            self.crtmVarList = ['Brightness Temperature', 'Radiance',
                                'Up Radiance', 'Down Radiance',
                                'Solar Downwelling Radiance']                       
            self.varBox.addItems(self.crtmVarList)
            self.varBox.setSizeAdjustPolicy(QComboBox.AdjustToContents)

            #Create channel combo box  
            channelBoxLabel = QLabel()
            channelBoxLabel.setText('Sensor Channel:')
            self.channelBox = QComboBox()
            self.channelBox.setStyleSheet(Layout.QComboBox())
            #Define sensor types
            self.sensor, self.channels, self.names = CRTMInfo.getSensorInfo(self.sensor_names[0])
            self.channelBox.addItems(self.names)
            self.channelBox.setSizeAdjustPolicy(QComboBox.AdjustToContents)
             
            #Create button to submit sensor and channel selections
            subButton = QPushButton()
            subButton.setStyleSheet(Layout.QPushButton3())
            subButton.resize(subButton.minimumSizeHint())
            subButton.setText('Submit')
            subButton.clicked.connect(self.crtmOptionsControl)

            #Add the Widgets to the layout
            crtmControlLayout.addWidget(sensorBoxLabel)
            crtmControlLayout.addWidget(self.sensorBox)
            crtmControlLayout.addWidget(varBoxLabel)
            crtmControlLayout.addWidget(self.varBox)
            crtmControlLayout.addWidget(channelBoxLabel)
            crtmControlLayout.addWidget(self.channelBox)
            crtmControlLayout.addWidget(subButton)
             
            #Display the selection window
            self.crtmControl.show()

        else:        
            self.readField()    
            self.setPlotVars()

    #This function creates the CRTM options tab
    #Needed because pop window only shows when first
    # switching to brightness temperature variable
    def crtmOptionsControl(self):

        #Close the pop window
        self.crtmControl.close()
        #Save selections
        sens = self.sensorBox.currentIndex()  
        chan = self.channelBox.currentIndex()
        var = self.varBox.currentIndex()

        #Define group box and give it a name
        boxTitleString = 'CRTM Control'
        self.crtmbox = QGroupBox()
        self.crtmbox.setStyleSheet(Layout.QGroupBox())

        crtmLayout = QVBoxLayout()
        self.crtmbox.setLayout(crtmLayout)

        #Sensor Box 
        sensorBoxLabel = QLabel()
        sensorBoxLabel.setText('Sensor Type: \n (Note: Visible channel calculations \n take longer than IR)')
        self.sensorBox = QComboBox()
        self.sensorBox.setStyleSheet(Layout.QComboBox())
        self.sensorBox.addItems(self.sensor_names)
        self.sensorBox.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        #Set index based on first user selection
        self.sensorBox.setCurrentIndex(sens)
        #Connect to function to change the channel list
        self.sensorBox.currentIndexChanged.connect(self.getNewSensorInfo)

        #Create Variable Selection combo box
        # Brightness Temperature or Radiance
        varBoxLabel = QLabel()
        varBoxLabel.setText('Output Variable:')
        self.varBox = QComboBox()
        self.varBox.setStyleSheet(Layout.QComboBox())
        #Define Variable Options
        self.varBox.addItems(self.crtmVarList)
        #Set index based on first user selection
        self.varBox.setCurrentIndex(var)
        self.varBox.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        
        #Create channel combo box
        channelBoxLabel = QLabel()
        channelBoxLabel.setText('Sensor Channel:')
        self.channelBox = QComboBox()
        self.channelBox.setStyleSheet(Layout.QComboBox())
        self.channelBox.addItems(self.names)
        self.channelBox.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        #Set index based on first user selection
        self.channelBox.setCurrentIndex(chan)

        #Replot button
        self.crtmPlotButton = QPushButton('Replot')
        self.crtmPlotButton.setStyleSheet(Layout.QPushButton3())
        self.crtmPlotButton.resize(self.crtmPlotButton.minimumSizeHint())
        self.crtmPlotButton.clicked.connect(self.selectionCRTM)

        #Connect widgets
        crtmLayout.addWidget(sensorBoxLabel)
        crtmLayout.addWidget(self.sensorBox)
        crtmLayout.addWidget(varBoxLabel)
        crtmLayout.addWidget(self.varBox)
        crtmLayout.addWidget(channelBoxLabel)
        crtmLayout.addWidget(self.channelBox)
        crtmLayout.addWidget(self.crtmPlotButton)

        #Add to options tab
        self.optionTabs.addTab(self.crtmbox,boxTitleString)
        self.optionTabs.setCurrentIndex(self.optionTabs.count()-1)
        self.tabbingLayout.addWidget(self.optionTabs)

        #Pass the initial selections on to complete plotting
        self.selectionCRTM()

    #This function sets the channels corresponding to the selected
    # sensor
    def getNewSensorInfo(self):
        
        #Determine the sensor index
        sindex = self.sensorBox.currentIndex()

        #Clear the current channel list
        self.channelBox.clear()
        
        #Add the new channel list
        #Get the new channels and names
        self.sensor, self.channels, self.names = CRTMInfo.getSensorInfo(self.sensor_names[sindex])
        self.channelBox.addItems(self.names)

    #This function controls the CRTM sensor and channel selections
    def selectionCRTM(self):
        
        #Get the sensor selection
        #self.crtmSensor = self.sensor[self.sensorBox.currentIndex()]
        self.crtmSensor = self.sensor        

        #Get the channel selection
        self.crtmChannel = self.channels[self.channelBox.currentIndex()]

        #Get the variable selection
        self.crtmVar = self.crtmVarList[self.varBox.currentIndex()]

        self.readField()
        self.setPlotVars()

    #This function is called to set plotting variables
    #  after a varaible (var or dvar) is changed.
    #  Also sets selections after the plot type is changed
    def setPlotVars(self):

        #Must be a 3D variable if on vertical slice
        if (self.currentPType == 'Vertical Slice' and len(self.var.shape) < 3):
            self.error3DVar()
        else:
            #If color scale is locked - unlock it
            if self.colorlock == True:
                self.colorlock = False
                self.unlock.setChecked(True) 
            #Setting these to None will force a re-scale in plot
            self.colormax = None
            self.colormin = None
            self.ncontours = None           
            
            #Force plotting of WRF radar reflectivity to pyART colors/scale
            #Not all datasets have a dvarlist so only check if WRF at the moment
            #Revist this as additional datasets are added (if they have derived vars)
            if (self.dataSet.dsetname == 'WRF'):
                #Need to check if a derived var has been selected before calling
                # something based on its index
                if (self.currentdVar != None):
                    #Clear CRTM controls if necessary
                    if (self.crtmbox is not None and 
                        self.dataSet.dvarlist[self.currentdVar] != 'BrightTemp/Radiance'):
                        self.crtmbox.setParent(None)
                        self.crtmbox = None
                    if (self.dataSet.dvarlist[self.currentdVar] == 'refl' and self.derivedVar == True):
                        self.extend = 'max'
                        self.colormin = 0
                        self.colormax = 80
                        self.ncontours = 41
                        #Lock the colorbar
                        self.colorlock = True
                        #Sets up user colorbar control if this is first variable
                        # Need it to initialize self.lock
                        if self.colorbox is None:
                            self.controlColorBar()
                        self.lock.setChecked(True)
                        if self.cmap != 'pyart_NWSRef':
                            self.cmap2 = self.cmap
                        self.cmap = 'pyart_NWSRef'
                    #Else is needed to change color back after being in reflectivity mode
                    else:
                        self.extend = 'both'
                        self.cmap = self.cmap2
                #Need an additional check for the non-derived reflectivity in WRF
                if (self.dataSet.variableList[self.currentVar] == 'REFL_10CM' and self.derivedVar == False):
                    self.extend = 'max'
                    self.colormin = 0
                    self.colormax = 80
                    self.ncontours = 41
                    #Lock the colorbar
                    self.colorlock = True
                    #Sets up user colorbar control if this is first variable
                    # Need it to initialize self.lock
                    if self.colorbox is None:
                        self.controlColorBar()
                    self.lock.setChecked(True)
                    if self.cmap != 'pyart_NWSRef':
                        self.cmap2 = self.cmap
                    self.cmap = 'pyart_NWSRef'
                #Else is needed to change color back after being in reflectivity mode
                else:
                    self.extend = 'both'
                    self.cmap = self.cmap2

            #Get the dimensions of the new variable
            ndim   = len(self.var.shape)
            #Replot horizontal slice for skewT, vertical profile and 
            #   time series plots 
            if (self.currentPType == 'SkewT/Hodograph' or 
                self.currentPType == 'Vertical Profile' or
                self.currentPType == 'Time Series'):
                self.replot2d = True
                self.col = None
                self.row = None
                #If time series also signal a new var
                if (self.currentPType == 'Time Series'):
                    self.newvar = True
            #Difference plot - read in difference data
            if (self.currentPType == 'Difference Plot'):
                self.readDiffField()
            #Make sure at least a 2D variables was selected
            if(ndim > 1) :
                #Set levels to 1 is 2D only
                # otherwise get 3D levels
                if(len(self.var.shape) == 2):
                    self.nz=1
                else:
                    self.nz=self.var.shape[0]
                #Clear level list and get new list
                self.selectLevel.clear()
                self.updateListView = False
                for j in range(0,self.nz):
                    self.selectLevel.addItem(str(j+1))
                self.updateListView = True
                #Call the plot function
                self.pltFxn(self.pNum)        
        
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
        self.selecteddVvar = None
        self.selectedVvar = i
        self.vplot = self.dataSet.readNCVariable(self.Vvarlist[i],
                     barbs=self.plotbarbs, contour2=self.plotcontour2)
        self.VertvarTitle = self.dataSet.description +' ('+self.dataSet.units
        self.VertvarTitle = self.VertvarTitle +') \n'+self.dataSet.getTime()

    def selectionChangedVerticalVar(self,i):
        self.selectedVvar = None
        self.selecteddVvar = i
        dVvar = wrf_dvar.WRFDerivedVar(dset = self.dataSet,
                                       var = self.dVvarlist[self.selecteddVvar],
                                       ptype = self.currentPType, sensor = self.crtmSensor,
                                       channel = self.crtmChannel,
                                       path = self.appobj.main_path,
                                       req_var = self.crtmVar)
        self.vplot = dVvar.var
        self.VertvarTitle = dVvar.varTitle

    #This function is called when the sfc statistics variable is changed
    def selectionChangedStatVar(self,i):
        self.selectedStatsVar = self.sfcList[i]
        #Time didnt change so set to False
        self.statsChangeTime = False
        #Reset color scale
        self.colormin = None
        self.colormax = None
        self.ncontours = None
        self.pltFxn(self.pNum)

    #This function is called when the sfc statistics type is changed
    def selectionChangedStatType(self,i):
        self.selectedStatType = self.statList[i]
        #Time didnt change so set to False
        self.statsChangeTime = False
        if (self.selectedStatType == 'Error'):
            self.cmap = self.cmap2
        #Reset color scale
        self.colormin = None 
        self.colormax = None 
        self.ncontours = None
        self.pltFxn(self.pNum)

    def selectionChangeParcel(self,i):
        self.selectedParcel = i
        if self.selectedParcel == 0:
            self.skewParcel = 'SB'
        elif self.selectedParcel == 1:
            self.skewParcel = 'ML'
        else:
            self.skewParcel = 'MU'

    def selectionChangeGrid(self,i):
        self.ColorBar = None
        self.cs = None
        self.cs2 = None
        self.barbs = None
        self.vectors = None
        self.vectorkey = None
        self.cs2label = None
        self.domain_average = None
        self.coasts = None
        self.countries = None
        self.states = None
        self.counties = None
        self.meridians = None
        self.parallels = None
        self.recallProjection = True
        self.currentGrid = i+1
        #Check if we need to change the map resolution
        if (self.dataSet.map[self.currentGrid-1].resolution != self.resolution):
            self.changeRes = True
        if (self.plotcities == True):    
            self.cities = False
            self.plotcities = False
            self.plotCities(plot=False)
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
        #Need this before adding the new variables
        #However, this gets called again later in readField so this
        # reads the data twice for now. Fix it later!
        self.dataSet.setGrid(self.currentGrid)
        self.selectTime.clear()
        self.selectTime.addItems(self.dataSet.timeList[self.currentGrid])
        self.selectVar.clear()
        self.dataSet.updateData()
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
        #If in time series - need to replot 2D field
        if (self.currentPType == 'Time Series'):
            self.replot2d = True  
            #Setting these to None keeps the previous plot from being
            # created again
            self.col = None
            self.row = None
        #Difference plot - read in difference data
        if (self.currentPType == 'Difference Plot'):
            self.diffdata.setTimeIndex(self.currentTime)
            self.diffdata.setGrid(self.currentGrid)
            self.readDiffField()
        self.axes1 = None
        self.figure.clear()
        #Need double conditional since var will be None if a variable
        # has not been selected
        if (self.currentPType == 'Vertical Slice' and
            (self.currentVar != None or self.currentdVar != None)):
            if (len(self.var.shape) < 3):
                self.error3DVar()
            else:
                self.pltFxn(self.pNum)
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
            if (self.currentPType == 'Difference Plot'):
                self.readDiffField()
            self.calTimeRange = False
            self.pltFxn(self.pNum)

    def selectionChangeDset(self,i):
        self.currentDset = i+1
        if self.colorbox is not None:
            self.colorbox.setParent(None)
            self.colorbox = None
            #If a colorbox exist, these will also exist
            #If color scale is locked - unlock it
            if self.colorlock == True:
                self.colorlock = False
            #Setting these to None will force a re-scale in plot
            self.colormax = None
            self.colormin = None
            self.ncontours = None
        if self.vslicebox is not None:
            self.vslicebox.setParent(None)
            self.vslicebox = None
        self.appobj.cbar.deleteLater()
        self.appobj.cbar = None
        self.cs = None
        self.cs2 = None
        self.barbs = None
        self.vectors = None
        self.vectorkey = None
        self.cs2label = None
        self.domain_average = None
        self.ColorBar = None
        self.coasts = None
        self.countries = None
        self.counties = None
        self.states = None
        self.derivedVar = False
        self.parallels = None
        self.meridians = None
        self.recallProjection = True
        self.currentGrid = 1
        self.currentTime = 0
        self.currentVar = None
        self.currentdVar = None
        #Check if switching from the Soundings dataset
        # Will have to change the plot type as Soundings only plot a Skew-T
        if (self.dataSet.dsetname == "Sounding"):
            #Switch to horizontal plot
            self.currentPType = 'Horizontal Slice' 
        #Plus 1 is needed because None is technically the first index
        self.dataSet = self.dSet[self.currentDset]
        #Check if switching to Soundings dataset
        # Will have to change the plot type as Soundings only plot a Skew-T
        if (self.dataSet.dsetname == "Sounding"):
            #Also replot map with station locations
            self.replot2d = True
            #Switch to horizontal plot
            self.currentPType = 'SkewT/Hodograph'
            #Need to call this after setting the plot
            # in order to initialize correctly
            self.getControlBar()        
        else:
            self.getControlBar()
            self.dataSet.setTimeIndex(self.currentTime)
            self.dataSet.setGrid(self.currentGrid, update=None)
            self.selectVar.clear()
            self.selectVar.addItems(self.dataSet.variableList)
        self.axes1 = None
        self.figure.clear()
        if 'runType' not in dir(self.dataSet):
            self.dataSet.runType = None
        self.pltFxn(self.pNum)

    #Function connected to the time list combobox
    def selectionChangeTime(self,i):
        self.currentTime = i
        self.dataSet.setTimeIndex(self.currentTime,update=False)
        self.timeChangeSettings()

    #Function connected to the next button
    def nxtButtonAction(self):
        self.currentTime+=1
        if self.currentTime == self.dataSet.ntimes*self.dataSet.getNumFiles():
            self.currentTime = 0
        self.selectTime.setCurrentIndex(self.currentTime)
        self.dataSet.setTimeIndex(self.currentTime,update=False)
        self.timeChangeSettings()	
 
    #Function connected to the previous button
    def prevButtonAction(self):
        self.currentTime-=1
        if self.currentTime == -1:
            self.currentTime = self.dataSet.getNumFiles()*self.dataSet.ntimes-1
        self.selectTime.setCurrentIndex(self.currentTime)
        self.dataSet.setTimeIndex(self.currentTime,update=False)
        self.timeChangeSettings()

    #All time change functions call this
    # function to set all of necessary plot
    # specific settings. Reduces redundancy.
    def timeChangeSettings(self):
        #Handle color scale settings
        if self.colorlock == False:
            #Reset colorbar settings
            self.colormax = None
            self.colormin = None
            self.ncontours = None
        #Replot 2d plot, set col and row to None so
        # a new window will not pop-up with changing time.
        #If we want to change the pop with time remove setting
        # self.col and self.row to None
        if (self.currentPType == 'SkewT/Hodograph' or
            self.currentPType == 'Vertical Profile' or 
            self.currentPType == 'Time Series'):
            self.replot2d = True
            self.col = None
            self.row = None
            if (self.currentPType == 'Vertical Profile' and self.selectedVvar != None):
                self.selectionChangeVerticalVar(self.selectedVvar)
            #Derived vertical profile
            if (self.currentPType == 'Vertical Profile' and self.selecteddVvar != None):
                self.selectionChangedVerticalVar(self.selecteddVvar)
        #GOES-R dataset handles getting the variable for plotting
        if (self.dataSet.dsetname == "GOES R" or self.dataSet.dsetname == "GOES Class"):
            #Have to update data here because readField isn't called
            self.dataSet.updateData()
            self.dataSet.advanceTime(self)
        else:
            self.readField()

        #Difference plot - read in difference data
        if (self.currentPType == 'Difference Plot'):
            print (self.currentTime)
            self.diffdata.setTimeIndex(self.currentTime)
            self.readDiffField()
                        
        #Make sure a variable is selected before plotting
        #Needed to allow changing of plot type 
        if (self.currentVar != None or self.currentdVar != None 
            or self.currentPType == 'Surface Statistics'):
            self.pltFxn(self.pNum)

    def selectionChangeOrient(self,i):
        self.currentOrient = i
        if self.orientList[self.currentOrient] == 'xz':
            if self.dataSet.runType == 'IDEAL':
                self.refboxLabel.setText('y-displacement:')
            else:
                self.refboxLabel.setText('Lat:')
        if self.orientList[self.currentOrient] == 'yz':
            if self.dataSet.runType == 'IDEAL':
                self.refboxLabel.setText('x-displacement:')
            else:
                self.refboxLabel.setText('Lon:')

    #Function to add wind barbs to the plot
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
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Function to add wind vectors to the plot
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
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Function to add second contour to the plot
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
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Function to toggle coastlines on/off
    def plotCoastlines(self):
        if self.plotcoasts == False:
            self.plotcoasts = True
        else:
            self.plotcoasts = False
            if self.coasts != None:
                self.coasts.remove()
                self.coasts = None
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Function to toggle counties on/off
    def plotCountries(self):
        if self.plotcountries == False:
            self.plotcountries = True
        else:
            self.plotcountries = False
            if self.countries != None:
                self.countries.remove()
                self.countries = None
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Function to toggle states on/off
    def plotStates(self):
        if self.plotstates == False:
            self.plotstates = True
        else:
            self.plotstates = False
            if self.states != None:
                self.states.remove()
                self.states = None
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Function to toggle counties on/off
    def plotCounties(self):
        if self.plotcounties == False:
            self.plotcounties = True
        else:
            self.plotcounties = False
            if self.counties != None:
                self.counties.remove()
                self.counties = None
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Function to toggle cities on/off (US only - read from utils/city_shapefiles/ directory)
    def plotCities(self,plot=None):
        if self.plotcities == False:
            self.plotcities = True
            self.dataSet.map[self.currentGrid-1].readshapefile(self.appobj.main_path + "utils/city_shapefile/cb_2016_us_ua10_500k",'cities',drawbounds=False)
        else:
            self.plotcities = False
            self.cities = False
            #if self.cities != None:
            #    self.cities.remove()
            #    self.cities = None
        self.calTimeRange = False
        if (plot is None):
            self.pltFxn(self.pNum)

    #Function that changes the background color option when user selected
    def changeBackground(self,i):
        #self.ColorBar = None
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
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    def changeResolution(self,i):
        self.recallProjection = True
        self.changeRes = True
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
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #This function calls plotting routine after changing the time index
    # Needed so the user doesn't have to click the map again just to 
    # change the time
    def getTimeSeriesStartEnd(self):
        self.pltFxn(self.pNum)

    #This function calls plotting routine after changing the stats options
    def getTimeStatsStartEnd(self):
        self.pltFxn(self.pNum)
        
    def enterPress(self):
        self.ref_pt = np.float(self.refbox.text())
        self.minpress = np.float(self.pbox.text())
        self.figure.clear()
        self.pltFxn(self.pNum)

    def enterCValue(self):
        self.colormin = np.float(self.selectMin.text())
        self.colormax = np.float(self.selectMax.text())
        self.ncontours = np.float(self.selectcontours.text())
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Function to define subsample interval
    def setSubsampleValue(self):

        #Create a box for the user to input the subsample interval
        self.subBox = QDialog(self.appobj)
        self.subBox.setStyleSheet(Layout.QGroupBox())
        subBoxControlLayout = QHBoxLayout(self.subBox)

        subLabel = QLabel()
        subLabel.setText('Subsample Interval:')
        self.selectSubsample = QLineEdit()
        self.selectSubsample.setStyleSheet(Layout.QLineEdit())
        self.selectSubsample.setText(str(self.xinterval))
        self.selectSubsample.returnPressed.connect(self.enterSubValue)

        #Add widgets
        subBoxControlLayout.addWidget(subLabel)
        subBoxControlLayout.addWidget(self.selectSubsample)
        self.subBox.show()

    #Function to set the subsample interval
    def enterSubValue(self):
        self.subBox.close()
        self.xinterval = np.int(self.selectSubsample.text())
        self.yinterval = np.int(self.selectSubsample.text())
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Change from contourf to pcolormesh
    def PColorMeshPlot(self):
        if self.filltype== "contourf":
            if self.cs != None:
                for coll in self.cs.collections:
                    coll.remove()
        self.parallels = None
        self.meridians = None
        self.cs = None
        self.domain_average = None
        self.filltype = "pcolormesh"
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Change plot from pcolormesh to contour
    def ContourFPlot(self):
        self.cs = None
        if self.domain_average != None:
           self.domain_average.remove()
        self.domain_average = None
        self.parallels = None
        self.meridians = None
        self.filltype = "contourf"
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #This function controls the color palette selection
    def selectionChangeColorPalette(self,i):
        self.changeColor = True
        self.cmap2 = self.appobj.colorlist[i]
        #Create cmap for Satellite Water Vapor
        if (self.cmap2 == 'sat_WV'):
            self.max_val = 273.
            self.min_val = 163.
            position = [0,(195.-self.min_val)/(self.max_val-self.min_val),
                        (223.-self.min_val)/(self.max_val-self.min_val),
                        (243.-self.min_val)/(self.max_val-self.min_val),
                        (261.-self.min_val)/(self.max_val-self.min_val),
                        (263.-self.min_val)/(self.max_val-self.min_val),1]
            ctable_path = os.path.join(self.appobj.main_path,'utils','colortables',self.cmap2)
            new_cmap = create_colormap.make_cmap(ctable_path,position=position,bit=False)
            plt.register_cmap(cmap=new_cmap)
        #Create cmap for Satellite IR
        elif (self.cmap2 == 'sat_IR'):
            self.max_val = 300.
            self.min_val = 170.
            position = [0,(200.-self.min_val)/(self.max_val-self.min_val),
                        (208.-self.min_val)/(self.max_val-self.min_val),
                          (218.-self.min_val)/(self.max_val-self.min_val),
                          (228.-self.min_val)/(self.max_val-self.min_val),
                          (245.-self.min_val)/(self.max_val-self.min_val),
                          (253.-self.min_val)/(self.max_val-self.min_val),
                          (258.-self.min_val)/(self.max_val-self.min_val),1]
            ctable_path = os.path.join(self.appobj.main_path,'utils','colortables',self.cmap2)
            new_cmap = create_colormap.make_cmap(ctable_path,position=position,bit=False)
            plt.register_cmap(cmap=new_cmap)
        #Create cmap for Cloud Albedo
        elif (self.cmap2 == 'cloud_albedo'):
            self.max_val = 1.0
            self.min_val = 0.0
            position = [0,(0.2-self.min_val)/(self.max_val-self.min_val),
                        (0.4-self.min_val)/(self.max_val-self.min_val),
                          (0.6-self.min_val)/(self.max_val-self.min_val),
                          (0.8-self.min_val)/(self.max_val-self.min_val),
                          1]
            ctable_path = os.path.join(self.appobj.main_path,'utils','colortables',self.cmap2)
            new_cmap = create_colormap.make_cmap(ctable_path,position=position,bit=False)
            plt.register_cmap(cmap=new_cmap)
        else:
            self.max_val = None
            self.min_val = None
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #Function to allow the user to turn on and off the lat/lon grid
    def plotLatLonGrid(self):
        if self.plotlatlon == False:
            self.plotlatlon = True
            for par in self.parallels:
                self.parallels[par][0][0].remove()
                if (len(self.parallels[par][1]) != 0):
                    self.parallels[par][1][0].remove()
            for mer in self.meridians:
                self.meridians[mer][0][0].remove()
                if (len(self.meridians[mer][1]) != 0):
                    self.meridians[mer][1][0].remove()
        else:
            self.plotlatlon = False
            for par in self.parallels:
                self.parallels[par][0][0].remove()
                if (len(self.parallels[par][1]) != 0):
                    self.parallels[par][1][0].remove()
            for mer in self.meridians:
                self.meridians[mer][0][0].remove()
                if (len(self.meridians[mer][1]) != 0):
                    self.meridians[mer][1][0].remove()
        self.parallels = None
        self.meridians = None
        self.calTimeRange = False
        self.pltFxn(self.pNum)

    #This function creates the popup so users can input
    # a map boundaries
    def createMapBoundaries(self):
        #Create a selection box
        self.mapBoundary = QDialog(self.appobj)
        self.mapBoundary.resize(self.mapBoundary.minimumSizeHint())
        self.mapBoundary.setMinimumHeight(self.appobj.screeny*.25)
        self.mapBoundary.setMinimumWidth(self.appobj.screenx*.25)
        self.mapBoundary.setWindowTitle("Set the map boundary")
        mapControlLayout = QGridLayout()
        self.mapBoundary.setLayout(mapControlLayout)

        #Create North/East/South/West Boxes
        northControl = QWidget()
        northLayout = QVBoxLayout()
        northControl.setLayout(northLayout)        
        northBox = QLabel()
        northBox.setText('North (latitude):')
        self.northLat = QLineEdit()
        self.northLat.setStyleSheet(Layout.QLineEdit())
        self.northLat.setText(str(self.north))
        northLayout.addWidget(northBox)
        northLayout.addWidget(self.northLat)

        eastControl = QWidget()
        eastLayout = QVBoxLayout()
        eastControl.setLayout(eastLayout)
        eastBox = QLabel()
        eastBox.setText('East (longitude):')
        self.eastLon = QLineEdit()
        self.eastLon.setStyleSheet(Layout.QLineEdit())
        self.eastLon.setText(str(self.east))
        eastLayout.addWidget(eastBox)
        eastLayout.addWidget(self.eastLon)

        southControl = QWidget()
        southLayout = QVBoxLayout()
        southControl.setLayout(southLayout)
        southBox = QLabel()
        southBox.setText('South (latitude):')
        self.southLat = QLineEdit()
        self.southLat.setStyleSheet(Layout.QLineEdit())
        self.southLat.setText(str(self.south))
        southLayout.addWidget(southBox)
        southLayout.addWidget(self.southLat)

        westControl = QWidget()
        westLayout = QVBoxLayout()
        westControl.setLayout(westLayout)
        westBox = QLabel()
        westBox.setText('West (longitude):')
        self.westLon = QLineEdit()
        self.westLon.setStyleSheet(Layout.QLineEdit())
        self.westLon.setText(str(self.west))
        westLayout.addWidget(westBox)
        westLayout.addWidget(self.westLon)        

        #Create an enter button
        enterButton = QPushButton('Enter')
        enterButton.setStyleSheet(Layout.QPushButton3())
        enterButton.resize(enterButton.minimumSizeHint())
        enterButton.clicked.connect(self.setMapBoundaries)

        #Create a button to set default values
        resetButton = QPushButton('Reset Domain')
        resetButton.setStyleSheet(Layout.QPushButton3())
        resetButton.resize(resetButton.minimumSizeHint())
        resetButton.clicked.connect(self.resetMapBoundaries)

        #Connect widgets
        mapControlLayout.addWidget(northControl,0,1)
        mapControlLayout.addWidget(westControl,1,0)
        mapControlLayout.addWidget(eastControl,1,2)
        mapControlLayout.addWidget(southControl,2,1)
        mapControlLayout.addWidget(enterButton,3,0)
        mapControlLayout.addWidget(resetButton,3,2)

        #Display the selection window
        self.mapBoundary.show()

    #Function to allow the user to set the map lat/lon boundaries
    def setMapBoundaries(self):
        #Check for input errors
        if (self.southLat.text() == 'None' or self.northLat.text() == 'None' or 
            self.westLon.text() == 'None' or self.eastLon.text() == 'None'): 
            #Display error message
            msg = QMessageBox(self.appobj)
            msg.setIcon(QMessageBox.Information)
            msg.setText('Please set all values before hitting enter')
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
        elif (np.float(self.southLat.text()) > 90. or np.float(self.southLat.text()) < -90. or
              np.float(self.northLat.text()) > 90. or np.float(self.northLat.text()) < -90.):
            #Display error message
            msg = QMessageBox(self.appobj)
            msg.setIcon(QMessageBox.Information)
            msg.setText('North and South must be between -90 and 90')
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
        elif (np.float(self.northLat.text()) <  np.float(self.southLat.text())):
            #Display error message
            msg = QMessageBox(self.appobj)
            msg.setIcon(QMessageBox.Information)
            msg.setText('North must be greater than South')
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
        elif (np.float(self.eastLon.text()) > 180. or np.float(self.eastLon.text()) < -180. or
              np.float(self.westLon.text()) > 180. or np.float(self.westLon.text()) < -180.):
            #Display error message
            msg = QMessageBox(self.appobj)
            msg.setIcon(QMessageBox.Information)
            msg.setText('East and West must be between -180 and 180')
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
        elif (np.float(self.eastLon.text()) <  np.float(self.westLon.text())):
            #Display error message
            msg = QMessageBox(self.appobj)
            msg.setIcon(QMessageBox.Information)
            msg.setText('East must be greater than West')
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
        else:
            self.userGrid = True
            self.south = np.float(self.southLat.text())
            self.north = np.float(self.northLat.text())
            self.west = np.float(self.westLon.text())
            self.east = np.float(self.eastLon.text())
            self.mapBoundary.close()
            self.recallProjection = True
            self.figure.clear()
            self.ColorBar = None
            self.cs = None
            self.cs2 = None
            self.barbs = None
            self.vectors = None
            self.vectorkey = None
            self.cs2label = None
            self.domain_average = None
            self.coasts = None
            self.countries = None
            self.states = None
            self.counties = None
            self.meridians = None
            self.parallels = None
            self.calTimeRange = False
            self.pltFxn(self.pNum)

    #Function to set values back to default of None
    def resetMapBoundaries(self):
        self.mapBoundary.close()
        self.userGrid = True
        self.south = None
        self.west = None
        self.north = None
        self.east = None
        self.recallProjection = True
        self.figure.clear()
        self.ColorBar = None
        self.cs = None
        self.cs2 = None
        self.barbs = None
        self.vectors = None
        self.vectorkey = None
        self.cs2label = None
        self.domain_average = None
        self.coasts = None
        self.countries = None
        self.states = None
        self.counties = None
        self.meridians = None
        self.parallels = None
        self.calTimeRange = False
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

    def errorAnomaly(self):
        msg = QMessageBox(self.appobj)
        msg.setIcon(QMessageBox.Information)
        msg.setText('The current dataset does not contain more than one value at the current time of day')
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
        self.currentPlot = 0
        self.plotCount = 0
        #self.plotTab = []
        self.create_main_frame()
        self.setWindowTitle('PyGEOMET')
        
        #Main Menu
        extractAction = QAction("&Exit",self)
        extractAction.setShortcut("Ctrl+Q")
        extractAction.setStatusTip("Leave the App")
        extractAction.triggered.connect(self.close_program)        
        
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

        #GOES R
        openGOESR = QAction("&GOES R", self)
        openGOESR.setShortcut("Ctrl+S")
        openGOESR.setStatusTip('Open GOES R Directory')
        openGOESR.triggered.connect(self.goesROpen)

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

        #Soundings
        openSound = QAction("&Soundings",self)
        openSound.setShortcut("Ctrl+S")
        openSound.setStatusTip('Open Sounding')
        openSound.triggered.connect(self.soundingOpen)

        #Generic netCDF files
        openNetCDF = QAction("&Generic netCDF",self)
        openNetCDF.setShortcut("Ctrl+C")
        openNetCDF.setStatusTip('Open netCDF')
        openNetCDF.triggered.connect(self.netcdfOpen)

        ####End Dataset menu bar ##############################

        # Get EOM
        self.getEOM = QAction('&Eq. of Motion', self)
        #Disable until a plot is added
        self.getEOM.setEnabled(False)
        self.getEOM.setStatusTip("Calculate EOM Components")
        self.getEOM.triggered.connect(self.EOMget)

        #Color Palettes
        self.colorlist = ['jet','brg','rainbow','bwr','RdBu','YlOrRd',
                          'viridis','magma','gray','pyart_NWSRef','sat_WV',
                          'sat_IR','cloud_albedo']
        jet = QAction("&Default (jet)",self)
        jet.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(0))
        brg = QAction(QIcon('brg.png'),"&BlueRedGreen",self)
        brg.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(1))
        rainbow = QAction("&Rainbow",self)
        rainbow.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(2))
        bwr = QAction("&BlueWhiteRed",self)
        bwr.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(3))
        rdbu = QAction("&RedBlue",self)
        rdbu.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(4))
        ylorrd = QAction("&YellowOrangeRed",self)
        ylorrd.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(5))
        viridis = QAction("&Viridis",self)
        viridis.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(6))
        magma = QAction("&Magma",self)
        magma.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(7))
        gray = QAction("&Gray",self)
        gray.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(8))
        ref = QAction("&NWS Reflectivity",self)
        ref.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(9))
        satWV = QAction("&Satellite Water Vapor",self)
        satWV.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(10))
        satIR = QAction("&Satellite IR",self)
        satIR.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(11))
        cldA = QAction("&Cloud Albedo",self)
        cldA.triggered.connect(lambda: self.slbplt[self.currentPlot].selectionChangeColorPalette(12))

        #Create map background menu bar
        defaultClear = QAction("&Default (Clear)", self)
        defaultClear.triggered.connect(lambda: self.slbplt[self.currentPlot].changeBackground(0))

        blueMarble = QAction("&Blue Marble", self)
        blueMarble.triggered.connect(lambda: self.slbplt[self.currentPlot].changeBackground(1))

        shadedRelief = QAction("&Shaded Relief", self)
        shadedRelief.triggered.connect(lambda: self.slbplt[self.currentPlot].changeBackground(2))

        topo = QAction("&Topo Map", self)
        topo.triggered.connect(lambda: self.slbplt[self.currentPlot].changeBackground(3))

        coarse = QAction("&Coarse", self)
        coarse.triggered.connect(lambda: self.slbplt[self.currentPlot].changeResolution(0))

        low = QAction("&Low", self)
        low.triggered.connect(lambda: self.slbplt[self.currentPlot].changeResolution(1))

        intermediate = QAction("&Intermediate", self)
        intermediate.triggered.connect(lambda: self.slbplt[self.currentPlot].changeResolution(2))

        high = QAction("&High", self)
        high.triggered.connect(lambda: self.slbplt[self.currentPlot].changeResolution(3))

        full = QAction("&Full", self)
        full.triggered.connect(lambda: self.slbplt[self.currentPlot].changeResolution(4))

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
        menuGOES.addAction(openGOESR)
        menuGOES.addAction(openGOESUAH)

        #Remaining Datasets
        dataSetMenu.addAction(openCMAQ)
        dataSetMenu.addAction(openNCEPNCAR)
        dataSetMenu.addAction(openMERRA)
        dataSetMenu.addAction(openNEXRAD)
        dataSetMenu.addAction(openSound)
        dataSetMenu.addAction(openNetCDF)

        #set up the plot settings menu
        self.plotMenu = mainMenu.addMenu('&Plotting Options')
        #Use the disable style sheet to make the text look "Grayed out"
        self.plotMenu.setStyleSheet(Layout.DisQMenu())

        #Calculations Menu
        self.calcMenu = mainMenu.addMenu('&Calculations')
        #Use the disable style sheet to make the text look "Grayed out"
        self.calcMenu.setStyleSheet(Layout.DisQMenu())
        menuEOM = self.calcMenu.addAction(self.getEOM)

        #Map settings
        self.mapMenu = self.plotMenu.addMenu('&Map Settings')
        #Disable until a plot is added
        self.mapMenu.setEnabled(False)

        #Map background options
        bckMenu = self.mapMenu.addMenu('&Change Map Background')
        bckMenu.setStyleSheet(Layout.QMenu())
        bckMenu.addAction(defaultClear)
        bckMenu.addAction(blueMarble)
        bckMenu.addAction(shadedRelief)
        bckMenu.addAction(topo)

        #Map Resolution options
        resMenu = self.mapMenu.addMenu('&Change Map Resolution')
        resMenu.setStyleSheet(Layout.QMenu())
        resMenu.addAction(coarse)
        resMenu.addAction(low)
        resMenu.addAction(intermediate)
        resMenu.addAction(high)
        resMenu.addAction(full)

        #User-defined Map Boundaries
        mc = QAction("Map Boundaries",self)
        mc.triggered.connect(lambda: self.slbplt[self.currentPlot].createMapBoundaries())
        mc.setStatusTip("Set the Lat/Lon boundaries for the Basemap")
        self.mapMenu.addAction(mc)

        self.pMenu  = self.plotMenu.addMenu('&Plot Settings')
        #Disable until a plot is added
        self.pMenu.setEnabled(False)
        #Colorbar options
        colorbarMenu = self.pMenu.addMenu('&Change Colorbar')
        colorbarMenu.setStyleSheet(Layout.QMenu())
        colorbarMenu.addAction(jet)
        colorbarMenu.addAction(brg)
        colorbarMenu.addAction(rainbow)
        colorbarMenu.addAction(bwr)
        colorbarMenu.addAction(rdbu)
        colorbarMenu.addAction(ylorrd)
        colorbarMenu.addAction(viridis)
        colorbarMenu.addAction(magma)
        colorbarMenu.addAction(gray)
        colorbarMenu.addAction(ref)
        satMenu = colorbarMenu.addMenu('&Satellite')
        satMenu.addAction(satWV)
        satMenu.addAction(satIR)
        satMenu.addAction(cldA)
 
        #Contour Fill options
        contourMenu = self.pMenu.addMenu('&Contour Fill Type')
        contourMenu.setStyleSheet(Layout.QMenu())
        cf = QAction("ContourF",self)
        cf.triggered.connect(lambda: self.slbplt[self.currentPlot].ContourFPlot())
        pc = QAction("PColorMesh",self)
        pc.triggered.connect(lambda: self.slbplt[self.currentPlot].PColorMeshPlot())
        contourMenu.addAction(cf)
        contourMenu.addAction(pc)
         
        #Overlay options
        self.overlayMenu = self.plotMenu.addMenu('&Overlays')
        #Disable until a plot is added
        self.overlayMenu.setEnabled(False)
        self.overlayMenu.setStyleSheet(Layout.QMenu())
        geogMenu = self.overlayMenu.addMenu('&Geography')
        geogMenu.setStyleSheet(Layout.QMenu())
        coast = QAction("Coastlines",self)
        coast.triggered.connect(lambda: self.slbplt[self.currentPlot].plotCoastlines())
        coast.setStatusTip('Overlay Coastlines')
        ctries = QAction("Countries",self)
        ctries.triggered.connect(lambda: self.slbplt[self.currentPlot].plotCountries())
        ctries.setStatusTip('Overlay Country Boundaries')
        states = QAction("States",self)
        states.triggered.connect(lambda: self.slbplt[self.currentPlot].plotStates())
        states.setStatusTip('Overlay State Boundaries')
        county = QAction("Counties",self)
        county.triggered.connect(lambda: self.slbplt[self.currentPlot].plotCounties())
        city = QAction("Cities",self)
        city.triggered.connect(lambda: self.slbplt[self.currentPlot].plotCities())
        city.setStatusTip('Overlay City Boundaries')
        llgrid = QAction("Lat/Lon Grid",self)
        llgrid.triggered.connect(lambda: self.slbplt[self.currentPlot].plotLatLonGrid())
        llgrid.setStatusTip('Overlay Lat/Lon Grid')
        wb = QAction("Wind Barbs",self)
        wb.triggered.connect(lambda: self.slbplt[self.currentPlot].plotWindBarbs())
        wb.setShortcut("Ctrl+B")
        wb.setStatusTip('Overlay Wind Barbs')
        wv = QAction("Wind Vectors",self)
        wv.triggered.connect(lambda: self.slbplt[self.currentPlot].plotWindVectors())
        wv.setShortcut("Ctrl+V")
        wv.setStatusTip('Overlay Wind Vectors')
        c2 = QAction("2nd Contour",self)
        c2.triggered.connect(lambda: self.slbplt[self.currentPlot].plotContour2())
        c2.setShortcut("Ctrl+C")
        c2.setStatusTip('Overlay 2nd Contour')
        geogMenu.addAction(coast)
        geogMenu.addAction(ctries)
        geogMenu.addAction(states)
        geogMenu.addAction(county)
        geogMenu.addAction(city)
        geogMenu.addAction(llgrid)
        self.overlayMenu.addAction(c2)
        self.overlayMenu.addAction(wb)
        self.overlayMenu.addAction(wv)

        #Subsample grid options
        self.subsample = QAction("Subsample Grid", self)
        #Disable until a plot is added
        self.subsample.setEnabled(False)
        self.subsample.triggered.connect(lambda: self.slbplt[self.currentPlot].setSubsampleValue())
        self.subsample.setStatusTip("Subsample the plot array")
        self.plotMenu.addAction(self.subsample)
        
        #### End Menu ####
        self.dataSet = [None]
        self.paths = [None]
        self.slbplt = [None]
        self.path = None 
        self.addplot = None
        self.cbar = None
        self.numDset = 0
        self.on_draw(None)

        #Path is needed to define CRTM coefficient data location
        self.main_path = os.path.abspath(__file__).split("main_GUI.py")[0]
 
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
        self.dname = 'GOES Class'
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

    #Select the files
    def goesROpen(self):
        self.dname = 'GOES R'
        #self.selectdataset.close()
        #self.path = QFileDialog.getExistingDirectory(self, 'Select Directory')
        self.files = QFileDialog.getOpenFileNames(self, 'Select Files')
        #self.prefix = 'OR*'
        #files = sorted(glob.glob(os.path.join(str(self.path),self.prefix)))
        #if self.path:
        if (len(self.files) > 0):        
            self.numDset += 1
            self.dataSet.append(GOESRDataset.GOESRDataset())
            self.paths.append("GOES R")
            QApplication.setOverrideCursor(Qt.WaitCursor)
            #self.dataSet[self.numDset].name(path=str(self.path) , prefix=self.prefix)
            self.dataSet[self.numDset].name(files=self.files)
            QApplication.restoreOverrideCursor()
        else:
            self.errorNoFilesFound()

    #Select the directory
    def goesOpen(self):
        self.dname = 'GOES UAH'
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

    #Open the sounding dataset
    def soundingOpen(self):
        self.dname = 'Sounding'
        #self.selectdataset.close()
        self.numDset += 1
        QApplication.setOverrideCursor(Qt.WaitCursor)
        self.dataSet.append(SoundDataset.SoundingDataset())
        QApplication.restoreOverrideCursor()

    #Select the files
    def netcdfOpen(self):
        self.dname = 'netCDF'
        self.files = QFileDialog.getOpenFileNames(self, 'Select Files')
        if (len(self.files) > 0):
            self.numDset += 1
            self.dataSet.append(NetCDFDataset.NetCDFDataset())
            QApplication.setOverrideCursor(Qt.WaitCursor)
            self.dataSet[self.numDset].name(files=self.files)
            QApplication.restoreOverrideCursor()
        else:
            self.errorNoFilesSelected()


    def EOMget(self):
        if self.dataSet[self.numDset].dsetname == "MERRA":
            self.eom = True
            self.recallProjection = False
            self.on_draw(self.currentPlot)
        else:
            self.EOMerror()
            
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

    #Create error message when no files are selected
    def errorNoFilesSelected(self):
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Information)
        msg.setText('No files were selected')
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
        #self.plotTab.append(QWidget())
        #self.plotControlTabs.addTab(self.plotTab[0],'GridView')
        self.cpLayout.addWidget(self.plotControlTabs)
 
        #qscroll = QScrollArea(self.main_frame)
        #qscroll.setStyleSheet(Layout.QScrollArea())

        #Make Plot Area into a Tab control
        self.plotAreaTab = QTabWidget(self.main_frame)
        #qlayout.addWidget(qscroll)
        qlayout.addWidget(self.plotAreaTab)
        qlayout.addWidget(self.controlPanel)

        #qscrollContents = QWidget()
        #qscrollLayout = QVBoxLayout(qscrollContents)

        #qscroll.setWidget(qscrollContents)
        #qscroll.setWidgetResizable(True)
        
        #qfigWidget = QWidget(qscrollContents)
        
        self.pltList = []
        #self.plotTab.addTab(self.gbox,boxTitleString)
        #self.plotLayout = QVBoxLayout()
        #qfigWidget.setLayout(self.plotLayout)
        #qscrollLayout.addWidget(qfigWidget)
        #qscrollContents.setLayout(qscrollLayout)
        self.setCentralWidget(self.main_frame)

        self.move(self.screenx*.05,self.screeny*.025)
        self.setStyleSheet(Layout.QMain())

        #Connect plot tabs to function
        self.plotControlTabs.currentChanged.connect(self.plotControlSelected)
        self.plotAreaTab.currentChanged.connect(self.plotAreaSelected)

    #Function to change the plot area selected with the
    # plot control change
    def plotControlSelected(self,i):
        self.plotAreaTab.setCurrentIndex(i)
        self.currentPlot = i + 1

    #Function to change the plot control selected with the
    # plot area change
    def plotAreaSelected(self,i):
        self.plotControlTabs.setCurrentIndex(i)
        self.currentPlot = i + 1

    def get_data2(self):
        return np.arange(20).reshape([4, 5]).copy()
    
    def addButtonAction(self):
        if (len(self.dataSet) >= 2):
            #Activate menu items and change layout
            if (self.plotControlTabs.count() == 0):
                self.mapMenu.setEnabled(True)
                self.pMenu.setEnabled(True)
                self.overlayMenu.setEnabled(True)
                self.subsample.setEnabled(True)
                self.getEOM.setEnabled(True)
                self.plotMenu.setStyleSheet(Layout.QMenu())
                self.calcMenu.setStyleSheet(Layout.QMenu())

            self.plotCount = self.plotCount + 1
            self.currentPlot = self.currentPlot + 1
            self.cw = CanvasWidget(self,self.plotCount)
            self.slbplt.append(PlotSlab(dset=self.dataSet,AppWid=self))
            self.slbplt[self.currentPlot].setConnection(self.on_draw,self.plotCount)
            self.slbplt[self.currentPlot].plotCount = self.plotCount
            #self.plotTab.append(QWidget())
            plotTab = QWidget()
            self.tabLayout = QVBoxLayout()
            #self.plotTab[self.plotCount].setLayout(self.tabLayout)
            plotTab.setLayout(self.tabLayout)
            self.slbplt[self.currentPlot].getControlBar()
            self.cw.setPlot(self.slbplt[self.currentPlot])
            #self.plotLayout.addWidget(self.cw)
            self.pltList.append(self.cw)
            tabName = 'Plot #'+ str(self.plotCount)
            self.plotAreaTab.addTab(self.cw,tabName)
            #self.plotControlTabs.addTab(self.plotTab[self.plotCount],tabName)
            self.plotControlTabs.addTab(plotTab,tabName)
            #Tab index starts at 0 so subtract 1 from currentPlot which starts at 1
            self.plotControlTabs.setCurrentIndex(self.currentPlot-1)
            self.plotAreaTab.setCurrentIndex(self.currentPlot-1)
            self.on_draw(self.plotCount)
        else:
            self.errorDset()

    def deleteButtonAction(self):
        if (self.plotControlTabs.count() != 0):
            #Delete the PlotObj from the list
            delete_ind = self.plotControlTabs.currentIndex()
            del self.slbplt[self.currentPlot]
            #Get number of tabs
            numtabs = self.plotAreaTab.count()
            self.plotAreaTab.removeTab(delete_ind)
            self.plotControlTabs.removeTab(delete_ind)
            #After tab removal, set the current index to the last tab if possible
            #If all tabs have been removed, disable some menu features
            if (self.plotControlTabs.count() != 0): 
                self.plotControlTabs.setCurrentIndex(self.plotControlTabs.count()-1)
                self.plotAreaTab.setCurrentIndex(self.plotAreaTab.count()-1)
            else:
                self.mapMenu.setEnabled(False)
                self.pMenu.setEnabled(False)
                self.overlayMenu.setEnabled(False)
                self.subsample.setEnabled(False)
                self.getEOM.setEnabled(False)
                self.plotMenu.setStyleSheet(Layout.DisQMenu())
                self.calcMenu.setStyleSheet(Layout.DisQMenu())

    def on_draw(self,plotnum):
        print("\n On draw: ", plotnum, "\n")
        if len(self.pltList) > 0:
            #if self.dname != 'NEXRAD Radar':
            #    self.pltList[plotnum-1].drawPlot() 
            #else:
            #    self.pltList[plotnum-1].drawPlot()
            self.pltList[plotnum-1].drawPlot()

        #self.main_frame.show()
        
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

class DataSelector():
    def testfxn(self,event):
        print ('leave_figure', event.canvas.figure)
        event.canvas.figure.patch.set_facecolor('blue')
    def __init__(self, plotobj, line):
        self.plotObj = plotobj
        self.points = 0
        self.marker = None
        self.markerEvent = None
        self.startMark = None
        self.visible = False
        self.line = line
        self.userSelects = []
        self.numSelection = 0
        self.mode = 'None'
        self.modeDefs = {'s':'Select', 'S':'Select','d':'Delete','D':'Delete',
                         'v':'View','V': 'View', 'm': 'Marker', 'escape': 'None'}
        self.subMode = 'None'
        self.subModeDefs = {'l':'Line', 'L':'Line','p':'Polygon','P':'Polygon'}
        #self.xs = list(line.get_xdata())
        #self.ys = list(line.get_ydata())
        self.xs = []
        self.ys = []
        #print('init',self.xs)
        #print("data selector")
        line.figure.canvas.setFocusPolicy( Qt.ClickFocus )
        line.figure.canvas.setFocus()
        #self.cid = line.figure.canvas.mpl_connect('button_press_event',self)
        self.cid = line.figure.canvas.mpl_connect('key_press_event', self.keyPress)

    def placeMarker(self,event):
        ix, iy = event.xdata, event.ydata
        print(ix,iy)
        self.marker = self.plotObj.dataSet.map[self.plotObj.currentGrid-1].plot(ix,iy,'bo',markersize=6)      
        self.plotObj.pltFxn(self.plotObj.pNum)
        
    def keyPress(self, event):
        print('in keypress')
        if event.key in self.modeDefs:
            self.mode = self.modeDefs[event.key]
            print('Mode = ',self.mode)
            if (self.mode == 'Marker'):
                self.markerEvent = self.line.figure.canvas.mpl_connect('button_press_event',self.placeMarker)
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
