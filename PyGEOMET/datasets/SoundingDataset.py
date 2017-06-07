import numpy as np
import netCDF4
import os
import datetime
import wget
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import PyGEOMET.utils
import PyGEOMET.utils.LayoutFormat as Layout

class SoundingDataset:

    # The constructor can be initialized by specifying the #
    # directory and file prefix for WRF output files.  If  #
    # the user does not provide this information when the  #
    # object is initalized, then set these variables to    #
    # None.                                                #

    def __init__(self):
        self.path = None
        self.grid = None
        #self.year = 2011
        #self.month = 04
        #self.day = 27
        #self.hour = 12
        self.var = None
        self.timeObj = None
        self.currentGrid = 1
        self.currentVarIndex = 1
        self.currentYearIndex = 0
        self.currentMonthIndex = -1
        self.currentDayIndex = -1
        self.currentHourIndex = -1
        self.currentMinuteIndex = -1
        self.currentTimeIndex = -1
        self.currentSweep = 0
        self.currentGridIndex = 67
        self.description = None
        self.dsetname = "Sounding"
        self.projectionType = "lcc"
        self.resolution = 'l'
        #Define plot type available for the dataset within the GUI
        self.ptypes = ['SkewT/Hodograph']
        self.stat_id = []
        self.stat_num = []
        self.obs_lon = []
        self.obs_lat = []
        self.getSoundLoc()
        self.setProjection()

    #Determine the sounding locations
    def getSoundLoc(self):
        path = os.path.dirname(PyGEOMET.utils.__file__)
        f = open(os.path.join(path,'sounding_locations.txt'))
        lin_num = 0
        for line in f:
            #Skip header
            if lin_num > 0:
                self.stat_id.append(line[0:4])
                self.obs_lon.append(float(line[12:20]))
                self.obs_lat.append(float(line[24:30]))
                self.stat_num.append(line[36:42])
            lin_num += 1
        f.close()
        #info = np.genfromtxt(f,skip_header=1,dtype=None) 
        #for i in range(len(info)):    
        #    self.stat_id.append(info[i][0])
        #    self.obs_lon.append(info[i][1])
        #    self.obs_lat.append(info[i][2])  
        #    self.stat_num.append(info[i][3])
 
        self.glons = [self.obs_lon]
        self.glats = [self.obs_lat]

    #Pull the selected sounding
    def getObsFile(self,ind,year=None,month=None,hour=None,station=None):
        if year == None or month == None or day == None:
            print("Using default time of: year=2011, month=04, day=28")
            year = '2011'
            month = '04'
            day = '28'
        if hour == None:
            hour = '00'

        #Define url
        url = "http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR="+year+"&MONTH="+month+"&FROM="+day+hour+"&TO="+day+hour+"&STNM="+self.stat_num[ind]

        #Get the file
        filename = wget.download(url)

        #Read in file
        p_o = []
        h_o = []
        t_o = []
        dew_o = []
        rh_o = []
        q_o = []
        direct_o = []
        wind_o = []
        theta_o = []

        f = open(filename,'r')
        count = 0
        count1 = -1
        count2 = -1
        switch = 0
        for line in f:
            #Check if there is any data
            if count == 1:
                if line[0:9] == "Can't get":
                    p_o.append(float('NaN'))
                    h_o.append(float('NaN'))
                    t_o.append(float('NaN'))
                    dew_o.append(float('NaN'))
                    q_o.append(float('NaN')) 
                    rh_o.append(float('NaN'))
                    direct_o.append(float('NaN'))
                    wind_o.append(float('NaN'))
                    theta_o.append(float('NaN'))
                    print("Missing data...")
                    break
            #Check if there is any data
            if count == 4:
                if line[3:10] == "Descrip":
                    p_o.append(float('NaN'))
                    h_o.append(float('NaN'))
                    t_o.append(float('NaN'))
                    dew_o.append(float('NaN'))
                    q_o.append(float('NaN'))
                    rh_o.append(float('NaN'))
                    direct_o.append(float('NaN'))
                    wind_o.append(float('NaN'))
                    theta_o.append(float('NaN'))
                    print("Missing data...")
                    break
        
            if count > 10:
                if line[0:6] == '</PRE>':
                    break
                p_o.append(float(line[0:7] if line[0:7] != '       '
                                else 'NaN'))
                h_o.append(float(line[9:14] if line[9:14] != '     '
                                else 'NaN'))
                t_o.append(float(line[16:21] if line[16:21] != '     '
                                else 'NaN'))
                dew_o.append(float(line[23:28] if line[23:28] != '     '
                                else 'NaN'))
                q_o.append(float(line[37:42] if line[37:42] != '     '
                                else 'NaN'))
                rh_o.append(float(line[32:35] if line[32:35] != '   '
                                else 'NaN'))
                direct_o.append(float(line[46:49] if line[46:49] != '   ' 
                                else 'NaN'))
                wind_o.append(float(line[53:56] if line[53:56] != '   ' 
                                else 'NaN'))
                theta_o.append(float(line[58:63] if line[58:63] != '     '
                                else 'NaN'))
            count += 1
        f.close()    

        os.remove(filename)

        self.p = np.array(p_o,dtype=np.float32)
        self.h = np.array(h_o,dtype=np.float32) 
        self.t = np.array(t_o,dtype=np.float32)
        self.dew = np.array(dew_o,dtype=np.float32)
        self.q = np.array(q_o,dtype=np.float32)
        self.rh = np.array(rh_o,dtype=np.float32)
        #Covert wind direction to "math" angle
        # from meteorological angle
        direct = 270. - np.array(direct_o,dtype=np.float32)
        #Convert wind to m/s from knots
        wind = np.array(wind_o,dtype=np.float32)*0.514444
        #Get the components of the wind
        self.u = wind * np.cos(direct*(np.pi/180.))
        self.v = wind * np.sin(direct*(np.pi/180.))
        self.theta = np.array(theta_o,dtype=np.float32)    

    def getTime(self):
        #utc = self.hourList[self.currentHourIndex] + ":" +\
        #      self.mmssList[self.currentHourIndex][self.currentMinuteIndex][0:2] + ":" +\
        #      self.mmssList[self.currentHourIndex][self.currentMinuteIndex][2:4]
        utc = '00:00'
        #self.timeString = self.month + '/' + self.day + '/' + self.year +\
        #                  ' ' + utc + ' UTC'
        self.timeString = '04/28/2011 00:00 UTC'
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

    def readNCVariable(self,vname, barbs=None, contour2=None):
        variable = np.squeeze(self.variables[vname][self.currentTimeIndex])
        if (hasattr(self.variables[vname],'units')):
            self.units = self.variables[self.var].units
        else:
            self.units = ''
        if (hasattr(self.variables[vname],'long_name')):
            self.description = self.variables[vname].long_name

    def setProjection(self,gid=None,axs=None):
        self.map = [None]*1
        self.nx = [None]
        self.ny = [None]
        #Hardcode for U.S.
        self.lat0 = [38.0]
        self.lon0 = [-95.0] 
        self.ll_lon = -122.0
        self.ll_lat = 20.0
        self.ur_lon = -63.0
        self.ur_lat = 50.0
        self.map[0] = Basemap(projection=self.projectionType,
                      lat_0 = self.lat0[0],
                      lon_0 = self.lon0[0], llcrnrlon = self.ll_lon,
                      llcrnrlat = self.ll_lat, urcrnrlon = self.ur_lon,
                      urcrnrlat = self.ur_lat, resolution=self.resolution,
                      ax = axs)

    def getDsetControlBar(self, plotObj):
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
        self.ptypes = ['SkewT']
        self.selectPlotType.addItems(self.ptypes)
        self.selectPlotType.setSizeAdjustPolicy(QComboBox.AdjustToContents)
#       self.selectPlotType.currentIndexChanged.connect(plotObj.selectionChangePlot)
        selectPlotWidgetLayout.addWidget(selectPlotLabel)
        selectPlotWidgetLayout.addWidget(self.selectPlotType)
        self.gboxLayout.addWidget(selectPlotWidget)

        #Time Control
        timeControlLabel = QLabel()
        timeControlLabel.setText('Time Control:')
        self.gboxLayout.addWidget(timeControlLabel)

        cpanel = QWidget()
        cpanelLayout = QHBoxLayout()
        cpanel.setLayout(cpanelLayout)
        nxtButton = QPushButton()
        nxtButton.setStyleSheet(Layout.QPushButton2())
        nxtButton.setText('&Next')
        nxtButton.setFixedWidth(75)
        nxtButton.clicked.connect(plotObj.nxtButtonAction)
        prevButton = QPushButton()
        prevButton.setStyleSheet(Layout.QPushButton2())
        prevButton.setText('&Prev')
        prevButton.setFixedWidth(75)
        prevButton.clicked.connect(plotObj.prevButtonAction)
        cpanelLayout.addWidget(prevButton)
        cpanelLayout.addWidget(nxtButton)
        self.gboxLayout.addWidget(cpanel)
        self.tabbingLayout.addWidget(self.gbox)
        self.optionTabs = QTabWidget()
        plotObj.qscrollLayout.addWidget(self.tabbing,Qt.AlignTop)        

        return self.qscroll

