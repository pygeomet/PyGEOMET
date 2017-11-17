import numpy as np
import netCDF4
import os
import datetime
from calendar import monthrange
import wget
import sys
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import PyGEOMET.utils
import PyGEOMET.utils.LayoutFormat as Layout

class SoundingDataset:

    def __init__(self):
        self.path = None
        self.grid = None
        self.var = None
        self.timeObj = None
        self.currentGrid = 1
        self.currentVarIndex = 1
        self.currentGridIndex = 67
        self.dsetname = "Sounding"
        self.projectionType = "lcc"
        self.resolution = 'l'
        #Define plot type available for the dataset within the GUI
        self.ptypes = ['SkewT/Hodograph']
        self.monthList = ['01','02','03','04','05','06','07',
                          '08','09','10','11','12']
        self.monNames = ['January', 'February', 'March', 'April', 'May',
                         'June', 'July', 'August', 'September', 'October',
                         'November', 'December']
        self.hourList = ['00','12']
        self.hourNames = ['00Z', '12Z']
        self.stat_id = []
        self.stat_num = []
        self.obs_lon = []
        self.obs_lat = []
        self.getSoundLoc()
        self.setProjection()
        #Flag to indicate station has been selected
        self.availStation = False

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
 
        self.glons = [self.obs_lon]
        self.glats = [self.obs_lat]

    #Pull the selected sounding
    def getObsFile(self,ind=None,year=None,month=None,day=None,hour=None,station=None):
        #Set flag to indicate station has been selected
        self.availStation = True
        #Either ind or station have to be set. 
        #The main_GUI will set ind from clicked selection
        #If used outside of the GUI then station is most likely provided
        if (ind == None and station == None):
            print("Error: pass a station index or ID")
            sys.exit()
        elif (station != None):
            #Determine index in file
            station_ind = np.where(np.array(self.stat_id) == station)[0]
            if (len(station_ind) == 0):
                print("Input station not found: "+station)
                sys.exit()
            stat = self.stat_num[station_ind]
        else:
            stat = self.stat_num[ind]
        if year == None or month == None or day == None:
            print("Using default time of: year=2011, month=04, day=27")
            print("Pass year, month and day to switch from default time")
            year = '2011'
            month = '04'
            day = '27'
        if hour == None:
            print("Using default hour of: 12 Z")
            hour = '12'

        #Define url
        url = "http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR="+year+"&MONTH="+month+"&FROM="+day+hour+"&TO="+day+hour+"&STNM="+stat

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
     
        #In the case of elevated stations check for first temperature
        # value
        ind_real = np.min(np.where(np.isnan(t_o) == False)[0])
        
        self.p = np.array(p_o[ind_real:],dtype=np.float32)
        self.h = np.array(h_o[ind_real:],dtype=np.float32) 
        self.t = np.array(t_o[ind_real:],dtype=np.float32)
        self.dew = np.array(dew_o[ind_real:],dtype=np.float32)
        self.q = np.array(q_o[ind_real:],dtype=np.float32)
        self.rh = np.array(rh_o[ind_real:],dtype=np.float32)
        #Covert wind direction to "math" angle
        # from meteorological angle
        direct = 270. - np.array(direct_o[ind_real:],dtype=np.float32)
        #Convert wind to m/s from knots
        wind = np.array(wind_o[ind_real:],dtype=np.float32)*0.514444
        #Get the components of the wind
        self.u = wind * np.cos(direct*(np.pi/180.))
        self.v = wind * np.sin(direct*(np.pi/180.))
        self.theta = np.array(theta_o[ind_real:],dtype=np.float32)    

    def getTime(self):
        self.timeString = self.month + '/' + self.day + '/' + self.year +\
                          ' ' + self.hour + ':00 UTC'
        return self.timeString

    #Note Update is a dummy variable to fit within GUI framework
    def setGrid(self, Indx, update = None):
        self.currentGridIndex = Indx
        #self.grid = self.gridList[self.currentGridIndex]

    def setProjection(self,gid=None,axs=None,west=None,north=None,east=None,south=None):

        #If there is input on grid extent change the corners
        if (west != None and north != None and east != None and south != None):
            llcrnrlon = west
            llcrnrlat = south
            urcrnrlon = east
            urcrnrlat = north
        else:
            #Hardcode for U.S. for now
            self.ll_lon = -122.0
            self.ll_lat = 20.0
            self.ur_lon = -63.0
            self.ur_lat = 50.0
            llcrnrlon = self.ll_lon
            llcrnrlat = self.ll_lat
            urcrnrlon = self.ur_lon
            urcrnrlat = self.ur_lat

        self.map = [None]*1
        self.nx = [None]
        self.ny = [None]
        #Hardcode for U.S.
        self.lat0 = [38.0]
        self.lon0 = [-95.0] 
        self.map[0] = Basemap(projection=self.projectionType,
                      lat_0 = self.lat0[0],
                      lon_0 = self.lon0[0], llcrnrlon = llcrnrlon,
                      llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon,
                      urcrnrlat = urcrnrlat, resolution=self.resolution,
                      ax = axs)

############# GUI Related Functions ##################

    #Change the year

    #Sets the year 
    def setYear(self,i):
        
        self.currentYearIndex = i
        #From the index get the year
        self.year = self.yearList[self.currentYearIndex]
        #If the year is the current year change the month
        # list to match the current available months
        old_month = self.month
        if (self.year == str(self.cyear)):
            old_month = self.month
            self.selectMonth.clear()
            for i in range(self.cmonth):
                self.selectMonth.addItem(self.monNames[i])
            #Try to set the month back to the same month before
            # switching the year. If it is in the future, start
            # at Jan.
            if (self.currentMonthIndex >= self.selectMonth.count()):     
                self.currentMonthIndex = 0
        else:
            self.selectMonth.clear()
            self.selectMonth.addItems(self.monNames) 

        #From the index get the month
        self.month = self.monthList[self.currentMonthIndex]
        self.selectMonth.setCurrentIndex(self.currentMonthIndex)

        #Check if the day list is still valid
        old_day = self.day
        if (self.month != old_month or self.month == str(self.cmonth).zfill(2)):
            self.selectDay.clear()
            #Reset the day list
            self.dayList = []
            #If the month is the current month - get current number
            # of days
            if (self.month == str(self.cmonth).zfill(2) and 
                self.year == str(self.cyear)):
                for i in range(self.cday):
                    self.dayList.append(str(i+1).zfill(2))
                self.selectDay.addItems(self.dayList)
            else: #Different month but not current month
                days = monthrange(int(self.year), int(self.month))[1]
                for i in range(days):
                    self.dayList.append(str(i+1).zfill(2))
                self.selectDay.addItems(self.dayList)
 
            #Check if the previous day selection is still available
            if (self.currentDayIndex >= self.selectDay.count()):
                self.currentDayIndex = 0
            #From the index get the day
            self.day = self.dayList[self.currentDayIndex]
            self.selectDay.setCurrentIndex(self.currentDayIndex)            

        #Check if the hour list is still valid
        if (self.year == str(self.cyear) and
            self.month == str(self.cmonth).zfill(2) and
            self.day == str(self.cday).zfill(2)):
            self.selectHour.clear()
            if (self.chour >= 12):
                self.selectHour.addItems(self.hourNames)
            else:
                self.selectHour.addItem(self.hourNames[0])            
            #Check if the previous hour selection is still available
            if (self.currentHourIndex >= self.selectHour.count()):
                self.currentHourIndex = 0
            #From the index get the hour
            self.hour = self.hourList[self.currentHourIndex]              
            self.selectHour.setCurrentIndex(self.currentHourIndex)

    #Sets the Month
    def setMonth(self,i):

        self.currentMonthIndex = i
        #From the index get the month
        self.month = self.monthList[self.currentMonthIndex]
        #If the month and year are the current month/year 
        # change the day list to match the current available days
        if (self.year == str(self.cyear) and self.month == str(self.cmonth).zfill(2)):
            self.selectDay.clear()
            #Create day list
            self.dayList = []
            #Initial list is based on current days in the current month
            for i in range(self.cday):
                self.dayList.append(str(i+1).zfill(2))
            self.selectDay.addItems(self.dayList)
        else: #Get the number of days in the month
            days = monthrange(int(self.year), int(self.month))[1]
            if (days != len(self.dayList)):
                self.selectDay.clear()
                #Create day list - need it to compare lengths
                self.dayList = []
                for i in range(days):
                    self.dayList.append(str(i+1).zfill(2))
                self.selectDay.addItems(self.dayList)
 
        #Check if the previous day selection is still available
        if (self.currentDayIndex >= self.selectDay.count()):
            self.currentDayIndex = 0
        #From the index get the day
        self.day = self.dayList[self.currentDayIndex]
        self.selectDay.setCurrentIndex(self.currentDayIndex)

        #Check if the hour list is still valid
        if (self.year == str(self.cyear) and
            self.month == str(self.cmonth).zfill(2) and
            self.day == str(self.cday).zfill(2)):
            self.selectHour.clear()
            if (self.chour >= 12):
                self.selectHour.addItems(self.hourNames)
            else:
                self.selectHour.addItem(self.hourNames[0])
            #Check if the previous hour selection is still available
            if (self.currentHourIndex >= self.selectHour.count()):
                self.currentHourIndex = 0
            #From the index get the hour
            self.hour = self.hourList[self.currentHourIndex]
            self.selectHour.setCurrentIndex(self.currentHourIndex)

    #Sets the Day
    def setDay(self,i):

        self.currentDayIndex = i
        #From the index get the day
        self.day = self.dayList[self.currentDayIndex]
        #If the month, year and day are the current month/day/year
        # change the hour list to match the current available hours
        if (self.year == str(self.cyear) and 
            self.month == str(self.cmonth).zfill(2) and 
            self.day == str(self.cday).zfill(2)):
            self.selectHour.clear()
            if (self.chour >= 12):
                self.selectHour.addItems(self.hourNames)
            else:
                self.selectHour.addItem(self.hourNames[0])
            #Check if the previous hour selection is still available
            if (self.currentHourIndex >= self.selectHour.count()):
                self.currentHourIndex = 0
            #From the index get the hour
            self.hour = self.hourList[self.currentHourIndex]
            self.selectHour.setCurrentIndex(self.currentHourIndex)

    #Sets the Hour
    def setHour(self,i):

        self.currentHourIndex = i
        #From the index get the hour
        self.hour = self.hourList[self.currentHourIndex]

    #Handles clicking of next button in the GUI
    def nxtButtonAction(self):

        #Check if you are already at the last available time
        if (self.year == str(self.cyear) and 
            self.month == str(self.cmonth).zfill(2) and
            self.day == str(self.cday).zfill(2) and
            self.hour == str(self.chour).zfill(2)):
            #At the end of the time period
            #Pop-up an error message notify the user
            self.errorTime('nxtButton')
        else: 
            #First change the hour index
            self.currentHourIndex += 1
            #Check if we need to switch to the next day
            if (self.currentHourIndex >= self.selectHour.count()):
                self.currentHourIndex = 0
                #Move the day forward
                self.currentDayIndex += 1
                #Check if the next day is available
                if (self.currentDayIndex >= self.selectDay.count()):
                    self.currentDayIndex = 0
                    #Move the month forward 
                    self.currentMonthIndex += 1
                    #Check if the next month is available
                    if (self.currentMonthIndex >= self.selectMonth.count()):
                        self.currentMonthIndex = 0
                        #Also set the year forward - dont need to check becuase
                        # it is captured in the initial if - reverse order
                        self.currentYearIndex -= 1
                        #Set the year
                        self.year = self.yearList[self.currentYearIndex]
                        self.selectYear.setCurrentIndex(self.currentYearIndex)
                        self.setYear(self.currentYearIndex)
                        #From the index get the hour
                        self.hour = self.hourList[self.currentHourIndex]
                        self.selectHour.setCurrentIndex(self.currentHourIndex)
                    else:                        
                        #Set the Month index
                        self.month = self.monthList[self.currentMonthIndex]
                        self.selectMonth.setCurrentIndex(self.currentMonthIndex)
                        #Reset the lists
                        self.setMonth(self.currentMonthIndex)
                        #From the index get the hour
                        self.hour = self.hourList[self.currentHourIndex]
                        self.selectHour.setCurrentIndex(self.currentHourIndex) 
                else:
                    #Set the day index
                    self.day = self.dayList[self.currentDayIndex]
                    self.selectDay.setCurrentIndex(self.currentDayIndex)
                    self.setDay(self.currentDayIndex)
                    #Set the index
                    self.hour = self.hourList[self.currentHourIndex]
                    self.selectHour.setCurrentIndex(self.currentHourIndex)
            else: 
                #Set the hour index
                self.hour = self.hourList[self.currentHourIndex]
                self.selectHour.setCurrentIndex(self.currentHourIndex) 

        #Create a new plot if station has been selected
        if (self.availStation):
            self.plotObj.appobj.cw.drawPlot()
 
    #Handles clicking of previous button in the GUI
    def prevButtonAction(self):

        #Check if you are already at the first available time
        if (self.year == '1985' and self.month == '01' and
            self.day == '01' and self.hour == '00'):
            #At the beginning of the time period
            #Pop-up an error message notify the user
            self.errorTime('prevButton')
        else:
            #First change the hour index
            self.currentHourIndex -= 1
            #Check if we need to switch to the next day
            if (self.currentHourIndex < 0):
                self.currentHourIndex = len(self.hourList)-1
                #Move the day back
                self.currentDayIndex -= 1
                #Check if the previous day is available
                if (self.currentDayIndex < 0):
                    #Move the month back
                    self.currentMonthIndex -= 1
                    #Check if the previous month is available
                    if (self.currentMonthIndex < 0):
                        self.currentMonthIndex = len(self.monthList)-1
                        #Can set the day index to 30 - always 31 days in Dec.
                        self.currentDayIndex = 30
                        #Also set the year back - dont need to check becuase
                        # it is captured in the initial if - reverse order
                        self.currentYearIndex += 1
                        #Set the year
                        self.year = self.yearList[self.currentYearIndex]
                        self.selectYear.setCurrentIndex(self.currentYearIndex)
                        self.setYear(self.currentYearIndex)
                        #From the index get the hour
                        self.hour = self.hourList[self.currentHourIndex]
                        self.selectHour.setCurrentIndex(self.currentHourIndex)
                    else:
                        #Set the Month index
                        self.month = self.monthList[self.currentMonthIndex]
                        self.selectMonth.setCurrentIndex(self.currentMonthIndex)
                        #From the month and year determine the number of days
                        # in the month
                        days = monthrange(int(self.year), int(self.month))[1]
                        #Set the index to the last day of the month
                        self.currentDayIndex = days - 1
                        #Reset the lists
                        self.setMonth(self.currentMonthIndex)
                        #From the index get the hour
                        self.hour = self.hourList[self.currentHourIndex]
                        self.selectHour.setCurrentIndex(self.currentHourIndex)
                else:
                    #Set the day index
                    self.day = self.dayList[self.currentDayIndex]
                    self.selectDay.setCurrentIndex(self.currentDayIndex)
                    self.setDay(self.currentDayIndex)
                    #Set the index
                    self.hour = self.hourList[self.currentHourIndex]
                    self.selectHour.setCurrentIndex(self.currentHourIndex)
            else:
                #Set the hour index
                self.hour = self.hourList[self.currentHourIndex]
                self.selectHour.setCurrentIndex(self.currentHourIndex)

        #Create a new plot if station has been selected
        if (self.availStation):
            self.plotObj.appobj.cw.drawPlot()

    #Error message
    def errorTime(self,location):
        msg = QMessageBox(self.plotObj.appobj)
        msg.setIcon(QMessageBox.Information)
        if (location == 'nxtButton'):
            msg.setText('At the end of time : \n' +\
                     self.month + '/' + self.day + '/' +\
                     self.year + '   ' + self.hour + ':00Z')
        elif (location == 'prevButton'):
            msg.setText('At the first available of time : \n' +\
                     self.month + '/' + self.day + '/' +\
                     self.year + '   ' + self.hour + ':00Z')            
        msg.setWindowTitle("Warning")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    #Control bar for GUI
    def getDsetControlBar(self, plotObj):
        self.plotObj = plotObj
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
                dsetname = os.path.basename(str(plotObj.dSet[count].path))
                self.selectDset.addItem(dsetname+' [Dataset '+str(count)+']')
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
        plotObj.selectPlotType = QComboBox()
        plotObj.selectPlotType.setStyleSheet(Layout.QComboBox())
        plotObj.selectPlotType.addItems(self.ptypes)
        plotObj.selectPlotType.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        selectPlotWidgetLayout.addWidget(selectPlotLabel)
        selectPlotWidgetLayout.addWidget(plotObj.selectPlotType)
        self.gboxLayout.addWidget(selectPlotWidget)

        #Time Control
        #Determine the current UTC time
        # Subtract 3 hours to account for delay in data from sounding
        #  launch time
        ctime = datetime.datetime.utcnow() - datetime.timedelta(hours=3)
        self.cyear = ctime.year
        self.cmonth = ctime.month
        self.cday = ctime.day
        self.chour = ctime.hour
        #Set the initial default time
        self.year = str(self.cyear).zfill(4)
        self.month = str(self.cmonth).zfill(2)
        self.day = str(self.cday).zfill(2)
        if (self.chour >= 12):
            #Reset chour - easier to compare to later
            self.chour = 12
            self.hour = '12'
        else:
            self.chour = 0
            self.hour = '00'

        timeControlLabel = QLabel()
        timeControlLabel.setText('Time Control [UTC]:')
        self.gboxLayout.addWidget(timeControlLabel)

        yearbar = QWidget()
        yearbarLayout = QHBoxLayout()
        yearbar.setLayout(yearbarLayout)
        yearWidgetLabel = QLabel()
        yearWidgetLabel.setText('Year:')
        self.selectYear = QComboBox()
        self.selectYear.setStyleSheet(Layout.QComboBox())
        self.selectYear.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        #Create year list
        self.yearList = []
        #Start year of 1985
        for i in range(self.cyear-1985+1):
            self.yearList.append(str(self.cyear-i))
        self.selectYear.addItems(self.yearList)
        self.selectYear.activated.connect(self.setYear)
        yearbarLayout.addWidget(yearWidgetLabel)
        yearbarLayout.addWidget(self.selectYear)
        self.gboxLayout.addWidget(yearbar)

        monthbar = QWidget()
        monthbarLayout = QHBoxLayout()
        monthbar.setLayout(monthbarLayout)
        monthWidgetLabel = QLabel()
        monthWidgetLabel.setText('Mon:')
        self.selectMonth = QComboBox()
        self.selectMonth.setStyleSheet(Layout.QComboBox())
        self.selectMonth.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        #Initial list is based on current days in the current month
        for i in range(self.cmonth):
            self.selectMonth.addItem(self.monNames[i])
        self.selectMonth.activated.connect(self.setMonth)
        monthbarLayout.addWidget(monthWidgetLabel)
        monthbarLayout.addWidget(self.selectMonth)
        self.gboxLayout.addWidget(monthbar)

        #Day
        daybar = QWidget()
        daybarLayout = QHBoxLayout()
        daybar.setLayout(daybarLayout)
        dayWidgetLabel = QLabel()
        dayWidgetLabel.setText('Day:')
        self.selectDay = QComboBox()
        self.selectDay.setStyleSheet(Layout.QComboBox())
        self.selectDay.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        #Create day list
        self.dayList = []
        #Initial list is based on current days in the current month
        for i in range(self.cday):
            self.dayList.append(str(i+1).zfill(2))        
        self.selectDay.addItems(self.dayList)
        self.selectDay.activated.connect(self.setDay)
        daybarLayout.addWidget(dayWidgetLabel)
        daybarLayout.addWidget(self.selectDay)
        self.gboxLayout.addWidget(daybar)

        #Hour
        hourbar = QWidget()
        hourbarLayout = QHBoxLayout()
        hourbar.setLayout(hourbarLayout)
        hourWidgetLabel = QLabel()
        hourWidgetLabel.setText('Hour:')
        self.selectHour = QComboBox()
        self.selectHour.setStyleSheet(Layout.QComboBox())
        self.selectHour.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        #Initial list is based on current hour in the current day
        if (self.chour >= 12):        
            self.selectHour.addItems(self.hourNames)
        else:
            self.selectHour.addItem(self.hourNames[0])
        self.selectHour.activated.connect(self.setHour)

        hourbarLayout.addWidget(hourWidgetLabel)
        hourbarLayout.addWidget(self.selectHour)
        self.gboxLayout.addWidget(hourbar)

        #Find the index to set the combobox displayed values
        self.currentYearIndex = np.where(np.array(self.yearList) == self.year)[0][0]
        self.currentMonthIndex = np.where(np.array(self.monthList) == self.month)[0][0]
        self.currentDayIndex = np.where(np.array(self.dayList) == self.day)[0][0]
        self.currentHourIndex = np.where(np.array(self.hourList) == self.hour)[0][0]
        #Set the index for each combobox
        self.selectYear.setCurrentIndex(self.currentYearIndex)
        self.selectMonth.setCurrentIndex(self.currentMonthIndex)
        self.selectDay.setCurrentIndex(self.currentDayIndex)
        self.selectHour.setCurrentIndex(self.currentHourIndex)

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

