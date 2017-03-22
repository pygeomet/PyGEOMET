import numpy as np
import pyart.graph
import pyart.io
import tempfile
import boto
import netCDF4
import os
import datetime
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import PyGEOMET.utils.NEXRADsites as NEXRADsites
import matplotlib.pyplot as plt
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import PyGEOMET.utils.LayoutFormat as Layout

class Radardataset:

    # The constructor can be initialized by specifying the #
    # directory and file prefix for WRF output files.  If  #
    # the user does not provide this information when the  #
    # object is initalized, then set these variables to    #
    # None.                                                #

    def __init__(self,path = None, prefix = None):
        self.path = None
        self.grid = None
        self.year = None
        self.month = None
        self.day = None
        self.hour = None
        self.mmss = None
        self.var = None
        self.timeObj = None
        self.swp = None
        self.currentGrid = 1
        self.currentVarIndex = 1
        self.currentYearIndex = 0
        self.currentMonthIndex = -1
        self.currentDayIndex = -1
        self.currentHourIndex = -1
        self.currentMinuteIndex = -1
        self.currentTimeIndex = -1
        self.currentSweep = 0
        self.currentGridIndex = 68 #69 is KHTX
        self.description = None
        self.dsetname = "NEXRAD Radar"
        self.gridList, latitudes, longitudes = NEXRADsites.get_sites()
        self.grid = self.gridList[self.currentGridIndex]
        self.projectionType = "lcc"
        self.resolution = 'i'
        self.cmap = 'pyart_NWSRef'
        self.range = [0,80]
        self.lon0 = [None]
        self.lat0 = [None]
        self.NEXRADfile(update=True)

    def setURL(self,update=None):
        #if no valid path, do nothing:
        if self.grid == None or self.year == None or \
           self.month == None:
            self.path = None
        else:
            self.path = self.year + '/' + self.month + '/' + self.day +\
                        '/' + self.grid + '/' + self.grid

    def NEXRADfile(self, update=None):

        if not hasattr(self,'yearList'):
            self.now = datetime.datetime.now()
            tmp = np.arange(self.now.year,2010-1,-1)
            self.yearList = ["%04d" % x for x in tmp]
            if self.year == None:
                self.year = self.yearList[self.currentYearIndex]

        if not hasattr(self, 'monList'):
            if int(self.year) == int(self.now.year):
                maxmon = self.now.month
            else:
                maxmon = 12
            tmp = np.arange(1,maxmon+1,1)
            self.monList = ["%02d" % x for x in tmp]
            if self.month == None:
                self.month = self.monList[self.currentMonthIndex]

        if not hasattr(self,'dayList'):
            if int(self.month) == int(self.now.month):
                maxday = self.now.day
            elif self.month == '02' and (((int(self.year)-1900) % 4) == 0):
                maxday = 29
            elif self.month == '02':
                maxday = 28
            elif self.month == '04' or self.month == '06' or \
                 self.month == '09' or self.month == 'self.11':
                maxday = 30
            else:
                maxday = 31
            tmp = np.arange(1,maxday+1,1)
            self.dayList = ["%02d" % x for x in tmp]

            if self.day == None:
                self.day = self.dayList[self.currentDayIndex]

        self.setURL()
        s3conn = boto.connect_s3()
        bucket = s3conn.get_bucket('noaa-nexrad-level2')
        flist = bucket.get_all_keys(prefix=self.path)
        self.ntimes = len(flist)
        hharr = []
        mmssarr = []
        self.farr = []
        j = -1
        for i in range(0,len(flist)):
            time = flist[i].name.split("_")[1]
            hh = time[0:2]
            mmss = time[2:6]
            self.farr.append(time[0:6])
            if hh not in hharr:
                mmssarr.append([])
                hharr.append(hh)
                j += 1
            mmssarr[j].append(mmss)

        self.hourList = hharr
        if self.hour == None:
            self.hour = self.hourList[self.currentHourIndex]
        self.mmssList = mmssarr
        if self.mmss == None:
            self.mmss = self.mmssList[self.currentMinuteIndex]

        self.setTimeIndex()

        s3key = bucket.get_key(flist[self.currentTimeIndex])
        print(s3key)
        localfile = tempfile.NamedTemporaryFile()
        s3key.get_contents_to_filename(localfile.name)
        self.radar = pyart.io.read_nexrad_archive(localfile.name)
 
        tmp = list(self.radar.sweep_number['data'])
        self.swpList = ["%02d" % x for x in tmp]

        if self.swp == None:
            self.swp = self.swpList[self.currentSweep]
        self.variableList = list(self.radar.fields.keys())
        if self.var == None:
            self.var = 'reflectivity'
            ind = np.where(np.asarray(self.variableList) == 'reflectivity')[0][0]
            self.currentVarIndex = ind
     
        self.lon0 = [self.radar.longitude['data'][0]]
        self.lat0 = [self.radar.latitude['data'][0]]

    def getTime(self):
        utc = self.hourList[self.currentHourIndex] + ":" +\
              self.mmssList[self.currentHourIndex][self.currentMinuteIndex][0:2] + ":" +\
              self.mmssList[self.currentHourIndex][self.currentMinuteIndex][2:4]
        self.timeString = self.month + '/' + self.day + '/' + self.year +\
                          ' ' + utc + ' UTC'
        return self.timeString

    def getTimeIndex(self):
        ind = np.where(np.array(self.hourList).astype(int) == int(self.farr[self.currentTimeIndex][0:2]))
        self.currentHourIndex = ind[0][0]
        self.hour = self.hourList[self.currentHourIndex]
        ind2 = np.where(np.array(self.mmssList[self.currentHourIndex]).astype(int) == int(self.farr[self.currentTimeIndex][2:6]))
        self.currentMinuteIndex = ind2[0][0]
        self.mmss = self.mmssList[self.currentHourIndex][self.currentMinuteIndex]

    def setTimeIndex(self):
        hr = self.hourList[self.currentHourIndex]
        ms = self.mmssList[self.currentHourIndex][self.currentMinuteIndex]
        self.currentTimeIndex = np.where(np.array(self.farr).astype(int) == np.array(hr+ms).astype(int))[0][0]
        self.getTime()

    def setGrid(self, Indx):
        self.currentGridIndex = Indx
        self.grid = self.gridList[self.currentGridIndex]

    def readNCVariable(self,vname, barbs=None, vectors=None,contour2=None):
        #if vname == 'velocity' and self.currentSweep == 0:
        #    sweep = self.currentSweep+1
        #else:
        #    sweep = self.currentSweep
        variable = np.squeeze(self.radar.get_field(sweep=self.currentSweep,field_name=vname))
        if vname == 'reflectivity':
            self.units = 'dBZ'
            self.description = 'Reflectivity'
            self.cmap = 'pyart_NWSRef'
            self.range = [0,80]
        if vname == 'velocity':
            self.units = 'm s$^-1$'
            self.description = 'Radial Velocity'
            self.cmap = 'pyart_NWSVel'
            self.range = [-30,30]

        if vname == 'differential_phase':
            self.units = 'deg'
            self.description = 'Differential Phase'
            self.cmap = 'pyart_RefDiff'
            self.range = [0,360]

        if vname == 'differential_reflectivity':
            self.units = 'dBZ'
            self.description = 'Differential Reflectivity'
            self.cmap = 'pyart_RefDiff'
            self.range = [-2,8]

        if vname == 'spectrum_width':
            self.units = 'm s$^{-1}$'
            self.description = 'Spectrum Width'
            self.cmap = 'pyart_NWS_SPW'
            self.range = [0,14]

        if vname == 'cross_correlation_ratio':
            self.units = 'ratio'
            self.description = 'Cross-Polar Correlation Coefficient'
            self.cmap = 'pyart_RefDiff'
            self.range = [0.5,1]
        print(vname)
        print(np.nanmax(variable),np.nanmin(variable))
        return variable 

    def setProjection(self,gid=None,axs=None):
        self.map = [None]*1
        self.nx = [None]
        self.ny = [None]
        self.glons = [None]*1
        self.glats = [None]*1
        self.ll_lon = self.lon0[0] - 2.25
        self.ll_lat = self.lat0[0] - 1.75
        self.ur_lon = self.lon0[0] + 2.25
        self.ur_lat = self.lat0[0] + 1.75
        self.map[0] = Basemap(projection=self.projectionType,
                      lat_0 = self.lat0[0],
                      lon_0 = self.lon0[0], llcrnrlon = self.ll_lon,
                      llcrnrlat = self.ll_lat, urcrnrlon = self.ur_lon,
                      urcrnrlat = self.ur_lat, resolution=self.resolution,
                      ax = axs)

        display = pyart.graph.RadarMapDisplay(self.radar)
        if ( self.var == 'velocity' or self.var == 'spectrum_width' ) and self.currentSweep == 0:
            self.currentSweep += 1
        if (self.var == 'differential_reflectivity' or 
            self.var == 'differential_phase' or 
            self.var == 'cross_correlation_ratio') and self.currentSweep == 1:
            self.currentSweep -= 1
        x,y = display._get_x_y(self.currentSweep,True,None)
        x0,y0 = self.map[0](self.lon0[0],self.lat0[0])
        self.glons[0],self.glats[0] = self.map[0]((x0+x*1000.),(y0+y*1000.),inverse=True)

    def getDsetControlBar(self, plotObj):
        self.plothook = plotObj
        self.tabbing = QWidget()
        plotObj.tabbingLayout = QVBoxLayout()
        self.tabbing.setLayout(plotObj.tabbingLayout)
        self.qscroll = QScrollArea()

        qscrollContents = QWidget()
        plotObj.qscrollLayout = QVBoxLayout(qscrollContents)

        self.qscroll.setWidget(qscrollContents)
        self.qscroll.setWidgetResizable(True)

        self.gbox = QGroupBox()
        self.gboxLayout = QVBoxLayout()
        self.gbox.setLayout(self.gboxLayout)

        #Dataset selection
        selectDsetWidget = QWidget()
        selectDsetWidgetLayout = QVBoxLayout()
        selectDsetWidget.setLayout(selectDsetWidgetLayout)
        selectDsetLabel = QLabel()
        selectDsetLabel.setText('Dataset:')

        self.selectDset = QComboBox()
        self.selectDset.setStyleSheet(Layout.QComboBox())
        count = 0
        for i in plotObj.dSet:
            if count == 0:
                count += 1
            else:
                self.selectDset.addItem(self.dsetname)
                count += 1
        self.selectDset.setCurrentIndex(plotObj.currentDset-1)
        self.selectDset.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.selectDset.currentIndexChanged.connect(plotObj.selectionChangeDset)
        selectDsetWidgetLayout.addWidget(selectDsetLabel)
        selectDsetWidgetLayout.addWidget(self.selectDset)
        self.gboxLayout.addWidget(selectDsetWidget)

        #Plot selection
        pltControlLabel = QLabel()
        pltControlLabel.setText('Plot Control:')
        self.gboxLayout.addWidget(pltControlLabel)

        selectPlotWidget = QWidget()
        selectPlotWidgetLayout = QHBoxLayout()
        selectPlotWidget.setLayout(selectPlotWidgetLayout)
        selectPlotLabel = QLabel()
        selectPlotLabel.setText('Plot Type:')
        self.selectPlotType = QComboBox()
        self.selectPlotType.setStyleSheet(Layout.QComboBox())
        ptypes = ['1-panel']#,'2-panel','4-panel']
        self.selectPlotType.addItems(ptypes)
        self.selectPlotType.setSizeAdjustPolicy(QComboBox.AdjustToContents)
#       self.selectPlotType.currentIndexChanged.connect(plotObj.selectionChangePlot)
        selectPlotWidgetLayout.addWidget(selectPlotLabel)
        selectPlotWidgetLayout.addWidget(self.selectPlotType)
        self.gboxLayout.addWidget(selectPlotWidget)

        #Time Control
        timeControlLabel = QLabel()
        timeControlLabel.setText('Time Control:')
        self.gboxLayout.addWidget(timeControlLabel)

        yearbar = QWidget()
        yearbarLayout = QHBoxLayout()
        yearbar.setLayout(yearbarLayout)
        yearWidgetLabel = QLabel()
        yearWidgetLabel.setText('Year:')
        plotObj.selectYear = QComboBox()
        plotObj.selectYear.setStyleSheet(Layout.QComboBox())
        plotObj.selectYear.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectYear.addItems(self.yearList)
        plotObj.selectYear.activated.connect(self.selectionChangeYear)
        yearbarLayout.addWidget(yearWidgetLabel)
        yearbarLayout.addWidget(plotObj.selectYear)
        self.gboxLayout.addWidget(yearbar)

        selectMDbar = QWidget()
        selectMDbarLayout = QHBoxLayout()
        selectMDbar.setLayout(selectMDbarLayout)
        monthWidgetLabel = QLabel()
        monthWidgetLabel.setText('Mon:')
        plotObj.selectMonth = QComboBox()
        plotObj.selectMonth.setStyleSheet(Layout.QComboBox())
        plotObj.selectMonth.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectMonth.addItems(self.monList)
        plotObj.selectMonth.activated.connect(self.selectionChangeMonth)
        selectMDbarLayout.addWidget(monthWidgetLabel)
        selectMDbarLayout.addWidget(plotObj.selectMonth)

        dayWidgetLabel = QLabel()
        dayWidgetLabel.setText('Day:')
        plotObj.selectDay = QComboBox()
        plotObj.selectDay.setStyleSheet(Layout.QComboBox())
        plotObj.selectDay.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectDay.addItems(self.dayList)
        plotObj.selectDay.activated.connect(self.selectionChangeDay)
        selectMDbarLayout.addWidget(dayWidgetLabel)
        selectMDbarLayout.addWidget(plotObj.selectDay)
        self.gboxLayout.addWidget(selectMDbar)

        selectMSbar = QWidget()
        selectMSbarLayout = QHBoxLayout()
        selectMSbar.setLayout(selectMSbarLayout)
        hourWidgetLabel = QLabel()
        hourWidgetLabel.setText('Hour:')
        plotObj.selectHour = QComboBox()
        plotObj.selectHour.setStyleSheet(Layout.QComboBox())
        plotObj.selectHour.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectHour.addItems(self.hourList)
        plotObj.selectHour.activated.connect(self.selectionChangeHour)
        selectMSbarLayout.addWidget(hourWidgetLabel)
        selectMSbarLayout.addWidget(plotObj.selectHour)

        minuteWidgetLabel = QLabel()
        minuteWidgetLabel.setText('MMSS:')
        plotObj.selectMMSS = QComboBox()
        plotObj.selectMMSS.setStyleSheet(Layout.QComboBox())
        plotObj.selectMMSS.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectMMSS.addItems(self.mmssList[self.currentHourIndex])
        plotObj.selectMMSS.activated.connect(self.selectionChangeMMSS)
        selectMSbarLayout.addWidget(minuteWidgetLabel)
        selectMSbarLayout.addWidget(plotObj.selectMMSS)
        self.gboxLayout.addWidget(selectMSbar)


        gridbar = QWidget()
        gridbarLayout = QHBoxLayout()
        gridbar.setLayout(gridbarLayout)
        gridWidgetLabel = QLabel()
        gridWidgetLabel.setText('Site:')
        plotObj.selectGrid = QComboBox()
        plotObj.selectGrid.setStyleSheet(Layout.QComboBox())
        plotObj.selectGrid.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectGrid.addItems(self.gridList)
        plotObj.selectGrid.activated.connect(self.selectionChangeGrid)
        gridbarLayout.addWidget(gridWidgetLabel)
        gridbarLayout.addWidget(plotObj.selectGrid)
        self.gboxLayout.addWidget(gridbar)

        swpbar = QWidget()
        swpbarLayout = QHBoxLayout()
        swpbar.setLayout(swpbarLayout)
        swpWidgetLabel = QLabel()
        swpWidgetLabel.setText('Sweep:')
        plotObj.selectSweep = QComboBox()
        plotObj.selectSweep.setStyleSheet(Layout.QComboBox())
        plotObj.selectSweep.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectSweep.addItems(self.swpList)
        plotObj.selectSweep.activated.connect(self.selectionChangeSweep)
        swpbarLayout.addWidget(swpWidgetLabel)
        swpbarLayout.addWidget(plotObj.selectSweep)
        self.gboxLayout.addWidget(swpbar)

        #Variable Control
        varbar = QWidget()
        varbarLayout = QHBoxLayout()
        varbar.setLayout(varbarLayout)
        varWidgetLabel = QLabel()
        varWidgetLabel.setText('Variable:')
        plotObj.selectVar = QComboBox()
        plotObj.selectVar.setStyleSheet(Layout.QComboBox())
        plotObj.selectVar.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectVar.addItems(self.variableList)
        plotObj.selectVar.activated.connect(self.selectionChangeVar)
        varbarLayout.addWidget(varWidgetLabel)
        varbarLayout.addWidget(plotObj.selectVar)
        self.gboxLayout.addWidget(varbar)

        plotb = QWidget()
        plotbLayout = QHBoxLayout()
        plotb.setLayout(plotbLayout)
        plotObj.plotButton = QPushButton()
        plotObj.plotButton.setStyleSheet(Layout.QPushButton2())
        plotObj.plotButton.setText('Plot')
        plotObj.plotButton.setFixedWidth(200)
        plotObj.plotButton.clicked.connect(self.plotButtonAction)
        plotbLayout.addWidget(plotObj.plotButton)
        self.gboxLayout.addWidget(plotb)

        cpanel = QWidget()
        cpanelLayout = QHBoxLayout()
        cpanel.setLayout(cpanelLayout)
        nxtButton = QPushButton()
        nxtButton.setStyleSheet(Layout.QPushButton2())
        nxtButton.setText('&Next')
        nxtButton.setFixedWidth(75)
        nxtButton.clicked.connect(self.nxtButtonAction)
        prevButton = QPushButton()
        prevButton.setStyleSheet(Layout.QPushButton2())
        prevButton.setText('&Prev')
        prevButton.setFixedWidth(75)
        prevButton.clicked.connect(self.prevButtonAction)
        cpanelLayout.addWidget(prevButton)
        cpanelLayout.addWidget(nxtButton)
        self.gboxLayout.addWidget(cpanel)
        plotObj.tabbingLayout.addWidget(self.gbox)
        plotObj.qscrollLayout.addWidget(self.tabbing,Qt.AlignTop)        

        return self.qscroll

    def selectionChangeYear(self,i):
        self.currentYearIndex = i
        self.year = self.yearList[self.currentYearIndex]
        if self.year == self.yearList[0]:
            self.currentYearIndex = 0
            self.currentMonthIndex = -1
            self.currentDayIndex = -1
            self.month = self.monList[self.currentMonthIndex]
            self.day = self.dayList[self.currentDayIndex]
        self.plothook.selectMonth.clear()
        self.plothook.selectDay.clear()
        self.plothook.selectHour.clear()
        self.plothook.selectMMSS.clear()

        if int(self.year) == int(self.now.year):
            maxmon = self.now.month
        else:
            maxmon = 12
        tmp = np.arange(1,maxmon+1,1)
        self.monList = ["%02d" % x for x in tmp]
        self.month = self.monList[self.currentMonthIndex]
        if int(self.month) == int(self.now.month):
            maxday = self.now.day
        elif self.month == '02' and (((int(self.year)-1900) % 4) == 0):
            maxday = 29
        elif self.month == '02':
            maxday = 28
        elif self.month == '04' or self.month == '06' or \
             self.month == '09' or self.month == 'self.11':
            maxday = 30
        else:
            maxday = 31
        tmp = np.arange(1,maxday+1,1)
        self.dayList = ["%02d" % x for x in tmp]
        self.day = self.dayList[self.currentDayIndex]
        self.plothook.selectMonth.addItems(self.monList)
        self.plothook.selectDay.addItems(self.dayList)

    def selectionChangeMonth(self,i):
        self.month = self.monList[i]
        if int(self.month) == int(self.now.month) and int(self.year) == int(self.now.year):
            maxday = self.now.day
        elif self.month == '02' and (((int(self.year)-1900) % 4) == 0):
            maxday = 29
        elif self.month == '02':
            maxday = 28
        elif self.month == '04' or self.month == '06' or \
             self.month == '09' or self.month == 'self.11':
            maxday = 30
        else:
            maxday = 31
        tmp = np.arange(1,maxday+1,1)
        self.dayList = ["%02d" % x for x in tmp]
        self.day = self.dayList[self.currentDayIndex]
        self.plothook.selectDay.clear()
        self.plothook.selectDay.addItems(self.dayList)

    def selectionChangeDay(self,i):
        self.day = self.dayList[i]
        self.NEXRADfile(update=True)
        self.plothook.selectVar.clear()
        self.plothook.selectHour.clear()
        self.plothook.selectMMSS.clear()
        self.plothook.selectVar.addItems(self.variableList)
        self.plothook.selectHour.addItems(self.hourList)
        self.plothook.selectMMSS.addItems(self.mmssList[self.currentHourIndex])

    def selectionChangeHour(self,i):
        self.currentHourIndex = i
        self.hour = self.hourList[self.currentHourIndex]
        self.plothook.selectMMSS.clear()
        self.plothook.selectSweep.clear()
        self.plothook.selectMMSS.addItems(self.mmssList[self.currentHourIndex])
        self.plothook.selectSweep.addItems(self.swpList)

    def selectionChangeMMSS(self,i):
        self.currentMinuteIndex = i
        self.minute = self.mmssList[self.currentHourIndex][self.currentMinuteIndex]

    def selectionChangeVar(self,i):
        self.currentVarIndex = i
        self.var = self.variableList[self.currentVarIndex]
        self.plothook.currentLevel = 0
        self.plothook.colormax = None
        self.plothook.colormin = None

    def selectionChangeGrid(self,i):
        self.plothook.selectVar.clear()
        self.plothook.selectSweep.clear()
        self.plothook.selectHour.clear()
        self.plothook.selectMMSS.clear()
        self.currentGridIndex = i
        self.grid = self.gridList[i]
        self.currentTimeIndex = 0
        self.currentVarIndex = 0
        self.var = None
        self.plothook.currentLevel = 0
        self.NEXRADfile(update=True)
        self.var = self.variableList[self.currentVarIndex]
        self.plothook.selectMMSS.addItems(self.mmssList[self.currentHourIndex])
        self.plothook.selectHour.addItems(self.hourList)
        self.plothook.selectVar.addItems(self.variableList)
        self.plothook.selectSweep.addItems(self.swpList)

    def selectionChangeSweep(self,i):
        self.plothook.currentSweep=i
        self.plothook.colormax = None
        self.plothook.colormin = None

    def plotButtonAction(self):
        self.NEXRADfile(update=True)
        self.setProjection(axs=self.plothook.appobj.axes1[self.plothook.pNum-1])
        self.plothook.nz = 1
        self.plothook.currentVar = self.currentVarIndex
        self.plothook.currentTime = self.currentTimeIndex
        self.plothook.readField()
        self.plothook.pltFxn(self.plothook.pNum)
#        if len(self.plothook.appobj.axes1) != 0:
#            self.plothook.appobj.axes1.remove(self.plothook.appobj.axes1[self.plothook.pNum-1])
#            self.plothook.figure.clear()
#        self.plotRadarData()#self.plothook.pNum)

    def nxtButtonAction(self):
        self.currentTimeIndex+=1
        errorflag = False
        if self.currentTimeIndex == self.ntimes:
            self.currentTimeIndex = 0
            self.currentYearIndex -= 1
            if self.currentYearIndex == -1:
                errorflag = True
                self.errorInvalidYear()
            else:
                self.year = self.yearList[self.currentYearIndex]
                self.selectionChangeYear(self.currentYearIndex)
                self.setURL(update=True)
                self.currentMonthIndex = 0
                self.month = self.monList[self.currentMonthIndex]
                self.currentDayIndex = 0
                self.day = self.dayList[self.currentDayIndex]
                self.currentHourIndex = 0
                self.hour = self.hourList[self.currentHourIndex]
                self.setTimeIndex()
        if not errorflag:
            self.getTimeIndex()
            self.plothook.currentTime = self.currentTimeIndex
            self.NEXRADfile(update=True)
            self.plothook.readField()
            self.plothook.pltFxn(self.plothook.pNum)
#            self.plotRadarData()
        else:
            pass

    def prevButtonAction(self):
        self.currentTimeIndex-=1
        errorflag = False
        if self.currentTimeIndex == -1:
            self.currentTimeIndex = -1
            self.currentYearIndex += 1
            if self.currentYearIndex == len(self.yearList):
                errorflag = True
                self.errorInvalidYear()
            else:
                self.year = self.yearList[self.currentYearIndex]
                self.selectionChangeYear(self.currentYearIndex)
                self.setURL(update=True)
                self.currentMonthIndex = -1
                self.month = self.monList[self.currentMonthIndex]
                self.currentDayIndex = -1
                self.day = self.dayList[self.currentDayIndex]
                self.currentHourIndex = -1
                self.hour = self.hourList[self.currentHourIndex]
                self.setTimeIndex()
        if not errorflag:
            self.getTimeIndex()
            self.plothook.currentTime = self.currentTimeIndex
            self.NEXRADfile(update=True)
            self.plothook.readField()
            self.plothook.pltFxn(self.plothook.pNum)
#            self.plotRadarData()
        else:
            pass

    def errorInvalidYear(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText('The AWS NEXRAD GUI is only set up from ' +\
                     self.yearList[-1] + '-' + self.yearList[0])
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def plotRadarData(self):
        self.pparams = [
        # variable-name in pyart, display-name that we want, sweep-number of radar (0=lowest ref, 1=lowest velocity)
	['reflectivity', 'Reflectivity (dBZ)', self.currentSweep, 0, 80, 10, 111]
#       ['velocity', 'Radial Velocity (m s-1)',self.currentSweep+1,-20,20,5, 212]
        ]
        self.tilt = "%.1f" % self.radar.fixed_angle['data'][self.currentSweep] + '$\mathsf{^{o}}$'
        for plotno, params in enumerate(self.pparams,start=1):
            figure = self.plothook.figure  
            if len(self.plothook.appobj.axes1) >= self.plothook.pNum:
                 self.plothook.appobj.axes1[self.plothook.pNum-1] = figure.add_subplot(params[6])
                 print("equal")
            else:
                 self.plothook.appobj.axes1.append(figure.add_subplot(params[6]))
                 print("not equal")
#            m = Basemap(projection=self.projectionType,
#                        lat_0 = self.lat0[0],
#                        lon_0 = self.lon0[0], llcrnrlon = self.ll_lon,
#                        llcrnrlat = self.ll_lat, urcrnrlon = self.ur_lon,
#                        urcrnrlat = self.ur_lat, resolution=self.resolution,
#                        ax = self.plothook.appobj.axes1[self.plothook.pNum-1])
#
            self.setProjection(axs=self.plothook.appobj.axes1[self.plothook.pNum-1])
            display = pyart.graph.RadarMapDisplay(self.radar)
            cs = display.plot_ppi_map(params[0],params[2], 
                                 colorbar_label=params[1],
                                 colorbar_flag=True,
                                 basemap = self.map[0], 
                                 vmin = params[3], vmax = params[4],
                                 lon_lines=[0], lat_lines=[0],
                                 ax=self.plothook.appobj.axes1[self.plothook.pNum-1])

            divider = make_axes_locatable(self.plothook.appobj.axes1[self.plothook.pNum-1])
            cax = divider.append_axes("right",size="5%",pad=0.1)
#            figure.colorbar(cs, cax=cax)
            cb = display.cbs[-1]
            cs =cb.mappable
            cb1 = figure.colorbar(cs, cax=cax)
            cb1.set_ticks(np.arange(params[3],params[4]+1,params[5]))

            self.map[0].drawcounties(color='darkslategray', linewidth=1)
            self.map[0].drawcoastlines(linewidth=2)
            self.map[0].drawstates(linewidth=2)
            display.plot_point(self.radar.longitude['data'][0], self.radar.latitude['data'][0],color='k')
        
        #self.plothook.appobj.axes1[self.plothook.pNum-1] = self.plothook.figure.add_subplot(111)
        self.plothook.appobj.axes1[self.plothook.pNum-1].set_title(self.grid + ' ' + self.tilt + ' ' + self.timeString, fontsize=16)
        print('apparently done with plotting')       

 
