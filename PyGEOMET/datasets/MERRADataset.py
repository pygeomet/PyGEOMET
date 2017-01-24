#MERRA dataset

from mpl_toolkits.basemap import Basemap
import numpy as np
import datetime
import netCDF4
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import os
import PyGEOMET.utils.wrf_functions as wrf
import PyGEOMET.utils.LayoutFormat as Layout

class MERRAdataset:

    def __init__(self, path = None, prefix=None):
        self.path = "http://goldsmr5.gesdisc.eosdis.nasa.gov/dods"
        self.grid = None
        self.year = None
        self.month = None
        self.day = None
        self.hour = None
        self.var = None
        self.timeObj = None
        self.currentGrid = 1
        self.ext = None
        self.currentVarIndex = 1
        self.currentYearIndex = 0
        self.currentMonthIndex = 0
        self.currentDayIndex = 0
        self.currentHourIndex = 0
        self.currentTimeIndex = 0
        self.currentGridIndex = 0
        self.currentExtIndex = 0
        self.description = None
        self.dsetname = "MERRA"
        self.gridList = ['M2I3NPASM', 'M2I3NVAER', 'M2I3NVASM', 'M2I3NVCHM', 
                         'M2I3NVGAS', 'M2I6NPANA', 'M2I6NVANA', 'M2IMNPANA',
                         'M2IMNPASM', 'M2IUNPANA', 'M2IUNPASM', 'M2T3NEMST',
                         'M2T3NENAV', 'M2T3NETRB', 'M2T3NPCLD', 'M2T3NPMST',
                         'M2T3NPODT', 'M2T3NPQDT', 'M2T3NPRAD', 'M2T3NPTDT',
                         'M2T3NPTRB', 'M2T3NPUDT', 'M2T3NVASM', 'M2T3NVCLD',
                         'M2T3NVMST', 'M2T3NVRAD', 'M2TMNPCLD', 'M2TMNPMST',
                         'M2TMNPODT', 'M2TMNPQDT', 'M2TMNPRAD', 'M2TMNPTDT',
                         'M2TMNPTRB', 'M2TMNPUDT', 'M2TUNPCLD','M2TUNPMST',
                         'M2TUNPODT', 'M2TUNPQDT', 'M2TUNPRAD', 'M2TUNPTDT',
                         'M2TUNPTRB', 'M2TUNPUDT']
        self.grid = self.gridList[self.currentGridIndex]
        self.glons = [None]*1
        self.glats = [None]*1
        self.projectionType=None
        self.resolution = "l"
        self.setURL(update=True)

    def setURL(self,update=None):
        #if no valid path, do nothing:
        if self.path == None or self.grid == None:
            self.url = None
        else:
            self.url = self.path + '/' + self.grid

        if update == True:
            self.MERRAfile(update=True)

    def MERRAfile(self, update=None):
        if hasattr(self,'ncId'):
            self.ncId.close()
        self.ncId = netCDF4.Dataset(self.url, mode='r')
        self.variables = self.ncId.variables
        self.dimensions = self.ncId.dimensions
        self.attributes = self.ncId.__dict__
        self.variableList = list(sorted(self.variables.keys()))
        if self.var == None:
            self.var = self.variableList[self.currentVarIndex]
        if 'lev' in self.variables[self.var].dimensions:
            tmp = np.squeeze(self.variables['lev'])
            if tmp.size > 1:
                self.nz = [tmp.size]
            else:
                self.nz = [1]
            if tmp.size == 1:
                self.levelList=[str(tmp)]
            else:
                self.levelList = ["%7.2f" % x for x in tmp]
                if self.grid[5] == "V":
                    self.levelList.reverse()
        else:
            self.levelList = ['NA']
        
        if not hasattr(self, 'oday') or not hasattr(self,'ohr'):
            #get time offset from datetime conversion
            reftime = datetime.datetime(1,1,1,0,0,0)
            fmt = "%Hz%d%b%Y"
            self.init = datetime.datetime.strptime(self.variables['time'].minimum, fmt)
            dtkday = (self.init - reftime).days  
            dtkhr = ((self.init-reftime).seconds/3600.)
            time0 = self.variables['time'][0]
            offset = time0 - (dtkday + (self.init - reftime).seconds/86400.)
            oday = int(np.floor(offset))
            ohr = int((offset - np.floor(offset))*24.)
            self.offset = datetime.timedelta(days = oday, hours = ohr)

        if not hasattr(self,'yearList'):
            lasttime = self.variables['time'][-1]
            dt = datetime.timedelta(days = lasttime)
            timeobj = reftime + dt - self.offset
            tmp = np.arange(timeobj.year-1, self.init.year,-1)
            self.yearList = ["%04d" % x for x in tmp] 
            if self.year == None:
                self.year = self.yearList[self.currentYearIndex]
        if not hasattr(self, 'monList'):
            maxmon = 12
            tmp = np.arange(1,maxmon+1,1)
            self.monList = ["%02d" % x for x in tmp]
            if self.month == None:
                self.month = self.monList[self.currentMonthIndex]

        if not hasattr(self,'dayList'):
            dayList = np.arange(1,31+1,1)
            if self.month == '02':
                tmp = dayList[0:28]
            elif self.month == '02' and (((int(self.year)-1900) % 4) == 0):
                tmp = dayList[0:29]
            elif self.month == '04' or self.month == '06' or self.month =='09' or self.month == '11':
                tmp = dayList[0:30]
            else:
                tmp = dayList
                self.dayList = ["%02d" % x for x in tmp]
            if self.day == None:
                self.day = self.dayList[self.currentDayIndex]
        if '3' in self.grid:
            tmp = np.arange(0,21+1,3)
            self.hourList = ["%02d" % x for x in tmp]
        elif '6' in self.grid:
            tmp = np.arange(0,18+1,6)
            self.hourList = ["%02d" % x for x in tmp]
        else:
            self.hourList = ["00"]
            self.hour = "00"
        if self.hour == None:
            self.hour = self.hourList[self.currentHourIndex]

        if update == True:
            self.readNCVariable(self.var)
            self.setProjection()

    def getTime(self):
        times = self.variables['time'][self.currentTimeIndex]
        self.ntimes = self.variables['time'].grads_size
        reftime = datetime.datetime(1,1,1,0,0,0)
        dt = datetime.timedelta(days = times)
        self.timeObj = reftime + dt - self.offset
        if self.timeObj.second == 59:
            addsecond = datetime.timedelta(seconds=1)
            self.timeObj += addsecond
        self.timeString = self.timeObj.strftime("%d %b %Y %H:%M:%S UTC")
        return self.timeString

    def setTimeIndex(self,update=None):
        reftime = datetime.datetime(1,1,1,0,0,0)
        selection = datetime.datetime(int(self.year),
                    int(self.month),int(self.day),
                    int(self.hour))
        dt = selection + self.offset -reftime
        days = dt.days + dt.seconds/86400.
        self.currentTimeIndex = np.where(np.squeeze(self.variables['time']) == days)[0][0]
        self.getTime()

    def setGrid(self, Indx,update=None):
        self.currentGridIndex = Indx
        self.grid = self.gridList[self.currentGridIndex]

    def readNCVariable(self,vname, barbs=None, vectors=None, contour2=None):
        variable = np.squeeze(self.variables[vname][self.currentTimeIndex])
        if (hasattr(self.variables[vname],'units')):
            self.units = self.variables[self.var].units
        else:
            self.units = ''
        if (hasattr(self.variables[vname],'long_name')):
            self.description = self.variables[vname].long_name

        if barbs == True or vectors == True:
            if self.grid[5] == 'P':
                inst = 'P'
            elif self.grid[5] == 'V':
                inst = 'V'
            else:
                print("This one hasn't been accounted for. Fix it")

            if (len(self.hourList) ==  8):
                gd = 'M2I3N' + inst + 'ASM'
            elif (len(self.hourList) == 4):
                gd = 'M2I6N' + inst + 'ANA'
            elif (len(self.hourList) == 1):
                gd = 'M2IMN' + inst + 'ASM'
            else:
                print("This one has not yet been accounted for. Fix it")
            uwindurl = self.path + "/" + gd
            tmp =netCDF4.Dataset(uwindurl)
            self.u10 = np.squeeze(tmp.variables['u'][self.currentTimeIndex])
            vwindurl = self.path + "/" + gd
            tmp =netCDF4.Dataset(vwindurl)
            self.v10 = np.squeeze(tmp.variables['v'][self.currentTimeIndex])

        return variable

    def getVertVars(self):
        if self.grid[5] == 'P':
            inst = 'P'
        elif self.grid[5] == 'V':
            inst = 'V'
        else:
            print("This one hasn't been accounted for. Fix it")
        if (len(self.hourList) ==  8):
            gd = 'M2I3N' + inst + 'ASM'
        elif (len(self.hourList) == 4):
            gd = 'M2I6N' + inst + 'ANA'
        elif (len(self.hourList) == 1):
            gd = 'M2IMN' + inst + 'ASM'
        else:
            print("This one has not yet been accounted for. Fix it")
        uwindurl = self.path + "/" + gd
        tmp =netCDF4.Dataset(uwindurl)
        u = np.squeeze(tmp.variables['u'][self.currentTimeIndex])
        v = np.squeeze(tmp.variables['v'][self.currentTimeIndex])
        tmp1 = np.squeeze(tmp.variables['omega'][self.currentTimeIndex])
        t = np.squeeze(tmp.variables['t'][self.currentTimeIndex])
        height = np.squeeze(tmp.variables['h'][self.currentTimeIndex])/1000.
        if inst == 'P':
            press = np.squeeze(tmp.variables['lev'])
        elif inst == 'V':
            press = np.squeeze(tmp.variables['pl'][self.currentTimeIndex])/100.

        rho = (press*100.)/(287.05*t)
        w = -tmp1/(rho*9.81)
        
        return u, v, w, press, height
 
    def setProjection(self,gid=None,axs=None):
        self.map = [None]*1
        lon = self.variables['lon']
        lat = self.variables['lat']
        self.glons[0], self.glats[0] = np.meshgrid(lon,lat)
        self.nx = [self.glons[0].shape[1]]
        self.ny = [self.glons[0].shape[0]]
        self.dx = [50000.]
        self.dy = [62500.]
        if self.projectionType == None or self.projectionType == 'robin':
            self.projectionType = 'robin'
            self.lon0 = [0.0]
            self.lat0 = [0.0]
            self.ll_lon = [-179.375]
            self.ll_lat = [-90]
            self.ur_lon = [180]
            self.ur_lat = [90]
            self.map[0] = Basemap(ax=axs, projection=self.projectionType, 
                          lat_0 = self.lat0[0],
	                  lon_0 = self.lon0[0], llcrnrlon = self.ll_lon[0],
	    	          llcrnrlat = self.ll_lat[0], urcrnrlon = self.ur_lon[0],
	    	          urcrnrlat = self.ur_lat[0], resolution=self.resolution)
        else:
            self.lon0 = [270.0]
            self.lat0 = [0.0]
            self.map[0] = Basemap(ax=axs, projection=self.projectionType,
                          lon_0 = self.lon0[0], boundinglat = self.lat0[0],
                          resolution=self.resolution)


    def getDsetControlBar(self, plotObj):
        self.plothook = plotObj
        self.tabbing = QWidget()
        self.tabbingLayout = QVBoxLayout()
        self.tabbing.setLayout(self.tabbingLayout)

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
        ptypes = ['Horizontal Slice', 'Vertical Slice', 'SkewT','Vertical Profile',
                  'Time Series', 'Difference Plot', 'Hovmoller Diagram']
        self.selectPlotType.addItems(ptypes)
        self.selectPlotType.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.selectPlotType.currentIndexChanged.connect(plotObj.selectionChangePlot)
        selectPlotWidgetLayout.addWidget(selectPlotLabel)
        selectPlotWidgetLayout.addWidget(self.selectPlotType)
        self.gboxLayout.addWidget(selectPlotWidget)

        #Change projection controls
        projectionControlLabel = QLabel()
        projectionControlLabel.setText("Projection Control:")
        self.gboxLayout.addWidget(projectionControlLabel)

        selectProjectionWidget = QWidget()
        selectProjectionWidgetLayout = QHBoxLayout()
        selectProjectionWidget.setLayout(selectProjectionWidgetLayout)
        selectDefault = QRadioButton("Default")
        selectDefault.setChecked(True)
        selectDefault.clicked.connect(lambda:self.DefaultProj())
        selectNH = QRadioButton("NH")
        selectNH.clicked.connect(lambda:self.NHProj())
        selectSH = QRadioButton("SH")
        selectSH.clicked.connect(lambda:self.SHProj())
        selectProjectionWidgetLayout.addWidget(selectDefault)
        selectProjectionWidgetLayout.addWidget(selectNH)
        selectProjectionWidgetLayout.addWidget(selectSH)
        self.gboxLayout.addWidget(selectProjectionWidget)

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

        selectMDHbar = QWidget()
        selectMDHbarLayout = QHBoxLayout()
        selectMDHbar.setLayout(selectMDHbarLayout)
        monthWidgetLabel = QLabel()
        monthWidgetLabel.setText('Mon:')
        plotObj.selectMonth = QComboBox()
        plotObj.selectMonth.setStyleSheet(Layout.QComboBox())
        plotObj.selectMonth.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectMonth.addItems(self.monList)
        plotObj.selectMonth.activated.connect(self.selectionChangeMonth)
        selectMDHbarLayout.addWidget(monthWidgetLabel)
        selectMDHbarLayout.addWidget(plotObj.selectMonth)

        dayWidgetLabel = QLabel()
        dayWidgetLabel.setText('Day:')
        plotObj.selectDay = QComboBox()
        plotObj.selectDay.setStyleSheet(Layout.QComboBox())
        plotObj.selectDay.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectDay.addItems(self.dayList)
        plotObj.selectDay.activated.connect(self.selectionChangeDay)
        selectMDHbarLayout.addWidget(dayWidgetLabel)
        selectMDHbarLayout.addWidget(plotObj.selectDay)

        hourWidgetLabel = QLabel()
        hourWidgetLabel.setText('Hour:')
        plotObj.selectHour = QComboBox()
        plotObj.selectHour.setStyleSheet(Layout.QComboBox())
        plotObj.selectHour.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectHour.addItems(self.hourList)
        plotObj.selectHour.activated.connect(self.selectionChangeHour)
        selectMDHbarLayout.addWidget(hourWidgetLabel)
        selectMDHbarLayout.addWidget(plotObj.selectHour)
        self.gboxLayout.addWidget(selectMDHbar)

        gridbar = QWidget()
        gridbarLayout = QHBoxLayout()
        gridbar.setLayout(gridbarLayout)
        gridWidgetLabel = QLabel()
        gridWidgetLabel.setText('Grid:')
        plotObj.selectGrid = QComboBox()
        plotObj.selectGrid.setStyleSheet(Layout.QComboBox())
        plotObj.selectGrid.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectGrid.addItems(self.gridList)
        plotObj.selectGrid.activated.connect(self.selectionChangeGrid)
        gridbarLayout.addWidget(gridWidgetLabel)
        gridbarLayout.addWidget(plotObj.selectGrid)
        self.gboxLayout.addWidget(gridbar)       

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

        levWidgetLabel = QLabel()
        levWidgetLabel.setText('Level:')
        plotObj.selectLev = QComboBox()
        plotObj.selectLev.setStyleSheet(Layout.QComboBox())
        plotObj.selectLev.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectLev.addItems(self.levelList)
        plotObj.selectLev.activated.connect(self.selectionChangeLevel)
        varbarLayout.addWidget(levWidgetLabel)
        varbarLayout.addWidget(plotObj.selectLev)
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
        self.tabbingLayout.addWidget(self.gbox)
        self.optionTabs = QTabWidget()
        self.optionTabs.setMaximumHeight(plotObj.appobj.screeny*.8*.4)
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

    def selectionChangeMonth(self,i):
        self.month = self.monList[i]
        dayList = np.arange(1,31+1,1)
        if self.month == '02' and (((int(self.year)-1900) % 4) == 0):
            tmp = dayList[0:29]
        elif self.month == '02':
            tmp = dayList[0:28]
        elif self.month == '04' or self.month == '06' or\
             self.month =='09' or self.month == '11':
            tmp = dayList[0:30]
        else:
            tmp = dayList
        self.dayList = ["%02d" % x for x in tmp]
        self.day = self.dayList[self.currentDayIndex]
        self.plothook.selectDay.clear()
        self.plothook.selectDay.addItems(self.dayList)

    def selectionChangeDay(self,i):
        self.day = self.dayList[i]

    def selectionChangeHour(self,i):
        self.hour = self.hourList[i]

    def selectionChangeVar(self,i):
        self.currentVarIndex = i
        self.var = self.variableList[self.currentVarIndex]
        self.plothook.currentLevel = 0
        self.plothook.colormax = None
        self.plothook.colormin = None

    def selectionChangeGrid(self,i):
        self.plothook.selectVar.clear()
        self.plothook.selectLev.clear()
        self.plothook.selectHour.clear()
        self.currentGridIndex = i
        self.grid = self.gridList[i]
        self.currentTimeIndex = 0
        self.currentVarIndex = 0
        self.var = None
        self.plothook.currentLevel = 0
        self.setURL(update=False)
        self.MERRAfile(update=False)
        self.var = self.variableList[self.currentVarIndex]
        self.plothook.selectHour.addItems(self.hourList)
        self.plothook.selectVar.addItems(self.variableList)
        self.plothook.selectLev.addItems(self.levelList)

    def selectionChangeLevel(self,i):
        self.plothook.currentLevel=i
        self.plothook.colormax = None
        self.plothook.colormin = None

    def plotButtonAction(self):
        self.setTimeIndex()
        self.setURL(update=True)
        self.plothook.nz = len(self.levelList)
        self.plothook.currentVar = self.currentVarIndex
        self.plothook.currentTime = self.currentTimeIndex
        self.plothook.readField()
        self.plothook.pltFxn(self.plothook.pNum)

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
            self.getTime()
            self.plothook.currentTime = self.currentTimeIndex
            self.plothook.readField()
            self.plothook.pltFxn(self.plothook.pNum)
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
            self.getTime()
            self.plothook.currentTime = self.currentTimeIndex
            self.plothook.readField()
            self.plothook.pltFxn(self.plothook.pNum)
        else:
            pass

    def DefaultProj(self):
        self.projectionType = None
        self.plothook.ColorBar = None
        self.plothook.appobj.cs = None
        self.plothook.appobj.cs2 = None
        self.plothook.appobj.barbs = None
        self.plothook.appobj.vectors = None
        self.plothook.appobj.vectorkey = None
        self.plothook.appobj.cs2label = None
        self.plothook.coasts = None
        self.plothook.countries = None
        self.plothook.states = None
        self.plothook.appobj.recallProjection = True
        self.plothook.appobj.axes1.remove(self.plothook.appobj.axes1[self.plothook.pNum-1])
        self.plothook.figure.clear()
        self.plothook.pltFxn(self.plothook.pNum)

    def NHProj(self):
        self.projectionType = 'npstere'
        self.plothook.ColorBar = None
        self.plothook.appobj.cs = None
        self.plothook.appobj.cs2 = None
        self.plothook.appobj.barbs = None
        self.plothook.appobj.vectors = None
        self.plothook.appobj.vectorkey = None
        self.plothook.appobj.cs2label = None
        self.plothook.coasts = None
        self.plothook.countries = None
        self.plothook.states = None
        self.plothook.appobj.recallProjection = True
        self.plothook.appobj.axes1.remove(self.plothook.appobj.axes1[self.plothook.pNum-1])
        self.plothook.figure.clear()
        self.plothook.pltFxn(self.plothook.pNum)
    
    def SHProj(self):
        self.projectionType = 'spstere'
        self.plothook.ColorBar = None
        self.plothook.appobj.cs = None
        self.plothook.appobj.cs2 = None
        self.plothook.appobj.barbs = None
        self.plothook.appobj.vectors = None
        self.plothook.appobj.vectorkey = None
        self.plothook.appobj.cs2label = None
        self.plothook.coasts = None
        self.plothook.countries = None
        self.plothook.states = None
        self.plothook.appobj.recallProjection = True
        self.plothook.appobj.axes1.remove(self.plothook.appobj.axes1[self.plothook.pNum-1])
        self.plothook.figure.clear()
        self.plothook.pltFxn(self.plothook.pNum)

    def errorInvalidYear(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText('The MERRA Reanalysis Dataset is only valid from ' +\
                     self.yearList[-1] + '-' + self.yearList[0])
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def calcPGF(self,col=None, row=None,level=None):
        if self.grid[5] == 'P':
            inst = 'P'
        elif self.grid[5] == 'V':
            inst = 'V'
        else:
            print("This one hasn't been accounted for. Fix it")
        if (len(self.hourList) ==  8):
            gd = 'M2I3N' + inst + 'ASM'
        elif (len(self.hourList) == 4):
            gd = 'M2I6N' + inst + 'ANA'
        elif (len(self.hourList) == 1):
            gd = 'M2IMN' + inst + 'ASM'
        else:
            print("This one has not yet been accounted for. Fix it")
        phi = np.squeeze(self.variables['h'][self.currentTimeIndex])*9.81
        if inst == 'P':
            press = np.squeeze(self.variables['lev'])*100.
        elif inst == 'V':
            press = np.squeeze(self.variables['pl'][self.currentTimeIndex])
            press = press[:,row,col]
        tmp = np.gradient(phi)
        if press[level] == press.min():
            dphidp = tmp[0][level][row][col]/(press[level-1] - press[level])
        else:
            dphidp = tmp[0][level][row][col]/((press[level] - press[level+1]))
        print(self.glons[0].shape, self.glats[0].shape)
        dy = wrf.get_distance(self.glats[0][row][col], self.glons[0][row][col],
                              self.glats[0][row+1][col], self.glons[0][row+1][col])
        dx = wrf.get_distance(self.glats[0][row][col], self.glons[0][row][col],
                              self.glats[0][row][col+1], self.glons[0][row][col+1])
        print("dx",dx)
        print("phi-x",tmp[2][level][row][col])
        print("dy",dy)
        print("phi-y",tmp[1][level][row][col])

        dphidy = tmp[1][level][row][col]/(dy*1000.)
        dphidx = tmp[2][level][row][col]/(dx*1000.)
        print("dphidy",dphidy)
        print("dphidx",dphidx)        


        return dphidx, dphidy, dphidp

    def calcCOR(self,col=None, row=None,level=None):
        if self.grid[5] == 'P':
            inst = 'P'
        elif self.grid[5] == 'V':
            inst = 'V'
        else:
            print("This one hasn't been accounted for. Fix it")
        if (len(self.hourList) ==  8):
            gd = 'M2I3N' + inst + 'ASM'
        elif (len(self.hourList) == 4):
            gd = 'M2I6N' + inst + 'ANA'
        elif (len(self.hourList) == 1):
            gd = 'M2IMN' + inst + 'ASM'
        else:
            print("This one has not yet been accounted for. Fix it")
        u = np.squeeze(self.variables['u'][self.currentTimeIndex])
        v = np.squeeze(self.variables['v'][self.currentTimeIndex])

        om = 360*np.pi/180./86164
        corx = 2*om*v[level][row][col]*np.sin(self.glats[0][row][col]*np.pi/180.)
        cory = -2*om*u[level][row][col]*np.sin(self.glats[0][row][col]*np.pi/180.)
        corz = 2*om*u[level][row][col]*np.cos(self.glats[0][row][col]*np.pi/180.)
        print("corx",corx)
        print("cory",cory)
        print("corz", corz)

        return corx, cory, corz

