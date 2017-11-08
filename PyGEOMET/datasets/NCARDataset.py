from mpl_toolkits.basemap import Basemap
import numpy as np
import datetime
import netCDF4
import requests
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import os
import PyGEOMET.utils.LayoutFormat as Layout

class NCARdataset:
    def __init__(self, path = None, prefix=None):
#http://www.esrl.noaa.gov/psd/thredds/fileServer/Datasets/ncep.reanalysis2/pressure/hgt.2002.nc
        self.path = "http://www.esrl.noaa.gov/psd/thredds/catalog/"
        self.path = "https://www.esrl.noaa.gov/psd/thredds/dodsC/"
        self.path += "Datasets/ncep.reanalysis2"
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
        self.currentMonthIndex = -1
        self.currentDayIndex = -1
        self.currentHourIndex = -1
        self.currentTimeIndex = -1
        self.currentGridIndex = 1
        self.currentExtIndex = 0
        self.description = None
        self.dsetname = "NCEP"
        self.gridList = ['gaussian_grid', 'pressure', 'surface']
        self.glons = [None]*1
        self.glats = [None]*1
        self.projectionType = None
        self.resolution = "l"
        self.NCARgetDDlists()

        #Define plot type available for the dataset within the GUI
        self.ptypes = ['Horizontal Slice']

    def setURL(self,update=None):
        #if no valid path, do nothing:
        if self.path == None or self.grid == None or self.var == None or\
        self.ext == None or self.year == None:
            self.url = None
        else:
            self.url = self.path + '/' + self.grid + '/' + self.var + "." +\
            self.ext + self.year + ".nc"

        #Check if the file exists
        fexist = requests.head(self.url+'.html')
        if (fexist.status_code >= 400):
            self.yearList = self.yearList[1:]
            self.year = self.yearList[self.currentYearIndex+1]
            self.url = self.path + '/' + self.grid + '/' + self.var + "." +\
            self.ext + self.year + ".nc"
        if update == True:
            self.NCARfile()

    def NCARgetDDlists(self):
        now = datetime.datetime.now()
        tmp = np.arange(now.year,1979-1,-1)
        self.yearList = ["%04d" % x for x in tmp]
        if self.year == None:
            self.year = self.yearList[self.currentYearIndex]

        if self.grid == None and self.var == None and self.timeObj == None:
            self.currentGridIndex = 1
            self.currentVarIndex = 1
            self.currentTimeIndex = -1
            self.grid = self.gridList[self.currentGridIndex]
        if self.grid == "gaussian_grid":
            self.variableList = ["air","cprat", "dlwrf", "dswrf", "ulwrf","uswrf",
                            "gflux", "icec", "lhtfl", "shtfl", "pevpr", "prate",
                            "pres", "runof","shum","skt","soilw", "tcdc", "tmax",
                            "tmin", "tmp", "uflx", "vflx", "ugwd", "vgwd",
                            "uwnd", "vwnd", "weasd"]
            self.var = self.variableList[self.currentVarIndex]

            if self.var == "air" or self.var == "shum" or\
               self.var == "tmax" or self.var == "tmin":
                self.extList = ["2m.gauss."]
            if self.var == "uwnd" or self.var == "vwnd":
                self.extList = ["10m.gauss."]
            if self.var == "cprat" or self.var == "dlwrf" or\
               self.var == "gflux" or self.var == "icec" or\
               self.var == "lhtfl" or self.var == "pevpr" or\
               self.var == "prate" or self.var == "runof" or\
               self.var == "shtfl" or self.var == "skt" or\
               self.var == "uflx" or self.var == "ugwd" or\
               self.var == "vflx" or self.var == "vgwd" or self.var == "weasd":

                self.extList = ["sfc.gauss."]
            if self.var == "tcdc":
                self.extList = ["eatm.gauss."]
            if self.var == "dswrf" or self.var == "ulwrf" or self.var == "uswrf":
                self.extList = ["sfc.gauss.","ntat.gauss."]
            if self.var == "tmp" or self.var == "soilw":
                self.extList = ["0-10cm.gauss.","10-200cm.gauss."]
            if self.var == "pres":
                self.extList = ["hcb.gauss.","hct.gauss.","lcb.gauss.","lct.gauss.",\
                                "mcb.gauss.", "mct.gauss.","sfc.gauss."]

        if self.grid == "pressure":
            self.variableList = ["air", "hgt", "omega", "shum", "tke", "uwnd", "vwnd"]
            self.var = self.variableList[self.currentVarIndex]
            self.extList = [""]

        if self.grid == "surface":
            self.variableList = ["mslp", "pr_wtr", "pres"]
            self.var = self.variableList[self.currentVarIndex]
            if self.var == "mslp":
                self.extList = [""]
            elif self.var == "pr_wtr":
                self.extList = ["eatm."]
            else:
                self.extList = ["sfc."]

        if self.ext == None:
            self.currentExtIndex = 0
        self.ext = self.extList[self.currentExtIndex]
        self.setURL(update=True)

    def NCARfile(self):
        if hasattr(self,'ncId'):
            self.ncId.close()
        self.ncId = netCDF4.Dataset(self.url)
        self.variables = self.ncId.variables
        self.attributes = self.ncId.__dict__
        if 'level' in self.variables.keys():
            tmp = np.squeeze(self.variables['level'])
            if tmp.size > 1:
                self.nz = [tmp.size]
            else:
                self.nz = [1]
            if tmp.size == 1:
                self.levelList=[str(tmp)]
            else:
                self.levelList = ["%4d" % x for x in tmp]
        else:
            self.levelList = ['NA']
        if not hasattr(self,'monList'):
            maxmon = datetime.datetime.strptime(self.getTime(),"%d %b %Y %H:%M:%S UTC").month
            tmp = np.arange(1,maxmon+1,1)
            self.monList = ["%02d" % x for x in tmp]
            if self.month == None:
                self.month = self.monList[self.currentMonthIndex]

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

            tmp = np.arange(0,18+1,6)
            self.hourList = ["%02d" % x for x in tmp]
            if self.hour == None:
                self.hour = self.hourList[self.currentHourIndex]

        self.readNCVariable(self.var)
        self.setProjection()

    def getTime(self):
        times = self.variables['time']
        self.ntimes = len(times)
        reftime = datetime.datetime(1800,1,1,0,0,0)
        hrs = times[self.currentTimeIndex]
        dt = datetime.timedelta(hours = hrs)
        self.timeObj = reftime + dt
        self.timeString = self.timeObj.strftime("%d %b %Y %H:%M:%S UTC")
        return self.timeString

    def setTimeIndex(self,update=None):
        reftime = datetime.datetime(1800,1,1,0,0,0)
        selection = datetime.datetime(int(self.year),
                    int(self.month),int(self.day),int(self.hour))
        dt = selection-reftime
        hrs = dt.days*24. + dt.seconds/3600.
        self.currentTimeIndex = np.where(np.asarray(self.variables['time']) == hrs)[0][0]
        self.getTime()

    def setGrid(self, Indx, update=None):
        self.currentGridIndex = Indx
        self.grid = self.gridList[self.currentGridIndex]
        self.NCARgetDDlists()

    def readNCVariable(self,vname,barbs=None, vectors=None, contour2=None):
        variable = np.squeeze(self.variables[vname][self.currentTimeIndex])
        if (hasattr(self.variables[vname],'var_desc')):
            self.description = self.variables[self.var].var_desc
        if (hasattr(self.variables[vname],'units')):
            self.units = self.variables[self.var].units
        if (hasattr(self.variables[vname],'long_name')):
            self.description = self.variables[vname].long_name

        if barbs == True or vectors == True:
            if (self.variables[vname].shape[1] != 1 and\
            len(self.variables[vname].shape) == 4):
                uwindurl = self.path + "/pressure/"
                uwindurl += "uwnd." + self.year + ".nc"
                tmp =netCDF4.Dataset(uwindurl)
                self.u10 = np.squeeze(tmp.variables['uwnd'][self.currentTimeIndex])
                vwindurl = self.path + "/pressure/"
                vwindurl += "vwnd." + self.year + ".nc"
                tmp =netCDF4.Dataset(vwindurl)
                self.v10 = np.squeeze(tmp.variables['vwnd'][self.currentTimeIndex])

            else:
                uwindurl = self.path + "/gaussian_grid/"
                uwindurl += "uwnd.10m.gauss." + self.year + ".nc"
                tmp =netCDF4.Dataset(uwindurl)
                self.u10 = np.squeeze(tmp.variables['uwnd'][self.currentTimeIndex])
                vwindurl = self.path + "/gaussian_grid/"
                vwindurl += "vwnd.10m.gauss." + self.year + ".nc"
                tmp =netCDF4.Dataset(vwindurl)
                self.v10 = np.squeeze(tmp.variables['vwnd'][self.currentTimeIndex])
        return variable

    def setProjection(self,gid=None,axs=None,west=None,north=None,east=None,south=None):
        self.map = [None]*1
        lon = self.variables['lon']
        lat = self.variables['lat']
        self.glons[0], self.glats[0] = np.meshgrid(lon,lat)
        self.nx = [self.glons[0].shape[1]]
        self.ny = [self.glons[0].shape[0]]
        self.dx = [100000.]
        self.dy = [100000.]
        if self.projectionType == None or self.projectionType == 'robin':
            self.projectionType = 'robin'
            self.lon0 = [0.0]
            self.lat0 = [0.0]
            self.ll_lon = [-180]
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
        self.selectPlotType.addItems(self.ptypes)
        self.selectPlotType.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.selectPlotType.currentIndexChanged.connect(plotObj.selectionChangePlot)
        selectPlotWidgetLayout.addWidget(selectPlotLabel)
        selectPlotWidgetLayout.addWidget(self.selectPlotType)
        self.gboxLayout.addWidget(selectPlotWidget)

        #Projection Control
        projectionControlLabel = QLabel()
        projectionControlLabel.setText("Projection Control:")
        self.gboxLayout.addWidget(projectionControlLabel)

        selectProjectionWidget = QWidget()
        selectProjectionWidgetLayout = QHBoxLayout()
        selectProjectionWidget.setLayout(selectProjectionWidgetLayout)
        selectDefault = QRadioButton("Default")
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
        hourWidgetLabel.setText('Hr:')
        plotObj.selectHour = QComboBox()
        plotObj.selectHour.setStyleSheet(Layout.QComboBox())
        plotObj.selectHour.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectHour.addItems(self.hourList)
        plotObj.selectHour.activated.connect(self.selectionChangeHour)
        selectMDHbarLayout.addWidget(hourWidgetLabel)
        selectMDHbarLayout.addWidget(plotObj.selectHour)
        self.gboxLayout.addWidget(selectMDHbar)

        #Grid Control
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

        extbar = QWidget()
        extbarLayout = QHBoxLayout()
        extbar.setLayout(extbarLayout)
        extWidgetLabel = QLabel()
        extWidgetLabel.setText('Surface:')
        plotObj.selectExt = QComboBox()
        plotObj.selectExt.setStyleSheet(Layout.QComboBox())
        plotObj.selectExt.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectExt.addItems(self.extList)
        plotObj.selectExt.activated.connect(self.selectionChangeExt)
        extbarLayout.addWidget(extWidgetLabel)
        extbarLayout.addWidget(plotObj.selectExt)
        self.gboxLayout.addWidget(extbar)

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

        self.plothook = plotObj

        return self.qscroll

    def selectionChangeYear(self,i):
        self.currentYearIndex = i
        self.year = self.yearList[self.currentYearIndex]
        self.setURL(update=False)
        if self.year == self.yearList[0]:
            self.currentYearIndex = 0
            self.currentMonthIndex = -1
            self.currentDayIndex = -1
            self.month = self.monList[self.currentMonthIndex]
            self.day = self.dayList[self.currentDayIndex]
            maxmon = datetime.datetime.strptime(self.getTime(),"%d %b %Y %H:%M:%S UTC").month
        else:
            maxmon = 12
        tmp = np.arange(1,maxmon+1,1)
        self.monList = ["%02d" % x for x in tmp]
        self.month = self.monList[self.currentMonthIndex]
        self.plothook.selectMonth.clear()
        self.plothook.selectMonth.addItems(self.monList)

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
        self.plothook.selectExt.clear()
        self.plothook.selectLev.clear()
        self.currentVarIndex = i
        self.var = self.variableList[self.currentVarIndex]
        if self.grid == "gaussian_grid":
            if self.var == "air" or self.var == "shum" or\
               self.var == "tmax" or self.var == "tmin":
                self.extList = ["2m.gauss."]
            if self.var == "uwnd" or self.var == "vwnd":
                self.extList = ["10m.gauss."]
            if self.var == "cprat" or self.var == "dlwrf" or\
               self.var == "gflux" or self.var == "icec" or\
               self.var == "lhtfl" or self.var == "pevpr" or\
               self.var == "prate" or self.var == "runof" or\
               self.var == "shtfl" or self.var == "skt" or\
               self.var == "uflx" or self.var == "ugwd" or\
               self.var == "vflx" or self.var == "vgwd" or self.var == "weasd":
                self.extList = ["sfc.gauss."]
            if self.var == "tcdc":
                self.extList = ["eatm.gauss."]
            if self.var == "dswrf" or self.var == "ulwrf" or self.var == "uswrf":
                self.extList = ["sfc.gauss.","ntat.gauss."]
            if self.var == "tmp" or self.var == "soilw":
                self.extList = ["0-10cm.gauss.","10-200cm.gauss."]
            if self.var == "pres":
                self.extList = ["hcb.gauss.","hct.gauss.","lcb.gauss.",\
                                "lct.gauss.","mcb.gauss.", "mct.gauss.",\
                                "sfc.gauss."]

        if self.grid == "pressure":
            self.extList = [""]

        if self.grid == "surface":
            if self.var == "mslp":
                self.extList = [""]
            elif self.var == "pr_wtr":
                self.extList = ["eatm."]
            else:
                self.extList = ["sfc."]
        self.currentExtIndex = 0
        self.ext = self.extList[self.currentExtIndex]
        self.plothook.colormax = None
        self.plothook.colormin = None
        self.plothook.selectExt.addItems(self.extList)
        self.plothook.selectLev.addItems(self.levelList)

    def selectionChangeExt(self,i):
        self.plothook.selectLev.clear()
        self.currentExtIndex = i
        self.ext = self.extList[self.currentExtIndex]
        self.setURL(update=False)
        self.plothook.selectLev.addItems(self.levelList)

    def selectionChangeGrid(self,i):
        self.plothook.selectVar.clear()
        self.plothook.selectExt.clear()
        self.plothook.selectLev.clear()
        self.currentGridIndex = i
        self.grid = self.gridList[i]
        self.currentVarIndex = 0
        self.currentextIndex = 0
        self.plothook.currentLevel=0
        if self.grid == "gaussian_grid":
            self.variableList = ["air","cprat","dlwrf", "dswrf", 
                                 "ulwrf","uswrf","gflux", "icec", 
                                 "lhtfl", "shtfl", "pevpr", "prate",
                                 "pres", "runof","shum","skt","soilw", 
                                 "tcdc", "tmax","tmin", "tmp", "uflx", 
                                 "vflx", "ugwd", "vgwd","uwnd", "vwnd", "weasd"]
            self.var = self.variableList[self.currentVarIndex]
            if self.var == "air" or self.var == "shum" or\
               self.var == "tmax" or self.var == "tmin":
                self.extList = ["2m.gauss."]
            if self.var == "uwnd" or self.var == "vwnd":
                self.extList = ["10m.gauss."]
            if self.var == "cprat" or self.var == "dlwrf" or\
               self.var == "gflux" or self.var == "icec" or\
               self.var == "lhtfl" or self.var == "pevpr" or\
               self.var == "prate" or self.var == "runof" or\
               self.var == "shtfl" or self.var == "skt" or\
               self.var == "uflx" or self.var == "ugwd" or\
               self.var == "vflx" or self.var == "vgwd" or self.var == "weasd":
                self.extList = ["sfc.gauss."]
            if self.var == "tcdc":
                self.extList = ["eatm.gauss."]
            if self.var == "dswrf" or self.var == "ulwrf" or self.var == "uswrf":
                self.extList = ["sfc.gauss.","ntat.gauss."]
            if self.var == "tmp" or self.var == "soilw":
                self.extList = ["0-10cm.gauss.","10-200cm.gauss."]
            if self.var == "pres":
                self.extList = ["hcb.gauss.","hct.gauss.","lcb.gauss.",\
                                "lct.gauss.","mcb.gauss.", "mct.gauss.",\
                                "sfc.gauss."]
            self.ext = self.extList[self.currentExtIndex]
        if self.grid == "pressure":
            self.variableList = ["air", "hgt", "omega", "shum", "tke", "uwnd", "vwnd"]
            self.var = self.variableList[self.currentVarIndex]
            self.extList = [""]
            self.ext = self.extList[self.currentExtIndex]

        if self.grid == "surface":
            self.levelList = ["sfc"]
            self.variableList = ["mslp", "pr_wtr", "pres"]
            self.var = self.variableList[self.currentVarIndex]
            if self.var == "mslp":
                self.extList = [""]
            elif self.var == "pr_wtr":
                self.extList = ["eatm."]
            else:
                self.extList = ["sfc."]
            self.ext = self.extList[self.currentExtIndex]

        self.var = self.variableList[self.currentVarIndex]
        self.plothook.selectVar.addItems(self.variableList)
        self.plothook.selectExt.addItems(self.extList)
        self.plothook.selectLev.addItems(self.levelList)

    def selectionChangeLevel(self,i):
        self.plothook.currentLevel=i
        self.plothook.colormax = None
        self.plothook.colormin = None

    def plotButtonAction(self):
        self.setURL(update=True)
        self.setTimeIndex()
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
        self.projectionType = 'robin'
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
        self.plothook.appobj.domain_average = None
        self.plothook.appobj.recallProjection = True
        self.plothook.appobj.axes1.remove(self.plothook.appobj.axes1[self.plothook.pNum-1])
        self.plothook.figure.clear()
        self.setProjection()
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
        self.plothook.appobj.domain_average = None
        self.plothook.appobj.recallProjection = True
        self.plothook.appobj.axes1.remove(self.plothook.appobj.axes1[self.plothook.pNum-1])
        self.plothook.figure.clear()
        self.setProjection()
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
        self.plothook.appobj.domain_average = None
        self.plothook.appobj.recallProjection = True
        self.plothook.appobj.axes1.remove(self.plothook.appobj.axes1[self.plothook.pNum-1])
        self.plothook.figure.clear()
        self.setProjection()
        self.plothook.pltFxn(self.plothook.pNum)

    def errorInvalidYear(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText('The NCAR Reanalysis Dataset is only valid from ' +\
                     self.yearList[-1] + '-' + self.yearList[0])
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

