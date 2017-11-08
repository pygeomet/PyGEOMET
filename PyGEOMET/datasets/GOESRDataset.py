import numpy as np
import os
import glob
import netCDF4
import datetime as dt
from mpl_toolkits.basemap import Basemap
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import PyGEOMET.utils.LayoutFormat as Layout

class GOESRDataset:

    # The constructor can be initialized by specifying the #
    # directory and file prefix for WRF output files.  If  #
    # the user does not provide this information when the  #
    # object is initalized, then set these variables to    #
    # None.                                                #

    def __init__(self,files = None):

        # Intialize a list to store list containing filenames #
        # for each domain.                                    #

        self.fileList = []

        self.timeList = [None]
        # Specify the WRF dataset by calling the function name()#

        self.name(files)

        # Current grid to which all operations apply !

        self.currentGrid = 1

        self.ncId = None

        self.variables = None

        self.attributes = None

        self.varName    = None

        self.ntimes = None

        self.currentTimeIndex = 0

        self.currentFileIndex = 1

        self.units='None'
	
        self.projectionType = 'cyl'
	
        self.pObj = None
          
        self.raw = None

        self.derived = None

        self.dvarlist = [None]
 
        self.dsetname = "GOES R"

        self.resolution = "l"

        #Define plot type available for the dataset within the GUI
        self.ptypes = ['Horizontal Slice']

    # name() is the function for specifying the WRF dataset. #
    # Based on the directory name and prefix for WRF output  #
    # files, this function will compile a list of file names #
    # for each WRF domain. 

    def name(self, files) :

        #Create dummy path - dataset name in GUI
        self.path = "/home/GOES_R"
        if files == None :

            self.numGrids = 0

        else :

            totalFiles = len(files)
            i = 1

            while(True):

                numFiles = len(files)
                times = []
                append_time = times.append 

                if numFiles != 0:
                    self.fileList = self.fileList + files[0]
                    self.setGrid(self.numGrids+1)
                    if i == 1:
                        dvar = ['Temperature']#, 'Effective Albedo']
                    else:
                        dvar = ['Temperature']#, 'Temperature']
                    self.dvarlist.append(dvar)
                    for ii in range(len(files[0])): 
                        self.setTimeIndex(ii)
                        append_time(self.getTime())
                    self.timeList.append(times)
                    self.numGrids += 1
                totalFiles = totalFiles-numFiles
                i = i+1
                if totalFiles == 0:
                    break;
                self.setTimeIndex(0,update=False)
            self.setTimeIndex(0,update=False)
            self.setGrid(1)
            self.variableReadFxn = { }
            for varname in self.variableList :
                self.variableReadFxn[varname] = self.readNCVariable
            self.setGridDefinition()
            self.rawvarlist = ['Radiance']

    def readNCVariable(self, vname) :
        variable = self.dsets[vname]
        if( hasattr(self.dsets[vname],'description')):
            self.description = self.dsets[vname].description
        if( hasattr(self.dsets[vname],'units')):
            self.units = self.dsets[vname].units
        return variable

    def getVariable(self,vname):
        return self.variableReadFxn[vname](vname)

    def addReadFxn(self,vname,fxn):
        self.variableReadFxn[vname] = fxn
        self.variableList.append(vname)

    def getAttribute(attname):
        return self.ncId.__dict__[attname]

    def setGridDefinition(self):

        self.nx  = []
        self.ny  = []
        self.band_num  = [None]
        self.band_wave = [None]
	#Needed for the GUI
        self.lat0  = [None]*self.numGrids
        self.lon0  = [None]*self.numGrids
        self.maplon = [None]*self.numGrids
        self.maplat = [None]*self.numGrids	
        self.map   = [None]*self.numGrids
        self.ur_lon = [0.0]*self.numGrids
        self.ur_lat = [0.0]*self.numGrids
        self.ll_lon = [0.0]*self.numGrids
        self.ll_lat = [0.0]*self.numGrids
        self.glats  = [None]*self.numGrids
        self.glons  = [None]*self.numGrids

        cg = self.currentGrid
        indx = self.currentFileIndex
        for i in range(0,self.numGrids):
            if i == 0:
                self.setTimeIndex(0,update=False)
                self.setGrid(i+1)
            else:
                self.setGrid(i+1,update=False)
                self.setTimeIndex(indx)

            self.setGrid(i+1)
            self.nx.append(self.ncId.dimensions['x'].size)
            self.ny.append(self.ncId.dimensions['y'].size)
            self.band_num.append(self.ncId.variables['band_id'][0])
            self.band_wave.append(self.ncId.variables['band_wavelength'][0])

        if len(self.nx) == self.numGrids:
            self.currentGrid = cg
            self.currentFileIndex = indx
	    	    
        self.setGridCorners()

    # getNumGrids() is a function that returns the maximum     #
    # number of grids detected by the object for the specified #
    # directory and file prefix.                               #


    def setGrid(self, GridNo, update = None):
        self.currentGrid = GridNo
        if update is None:
            self.updateData()

    def getNumFiles(self):
        return len(self.fileList)

    def getFileList(self):
        return(self.fileList[self.currentGrid])

    def updateData(self):
        print("Update data", self.fileList[self.currentFileIndex])
        self.ncId = netCDF4.Dataset(self.fileList[self.currentFileIndex],'r')
        self.variables = self.ncId.variables
        self.attributes = self.ncId.__dict__
        self.variableList = list(sorted(self.variables.keys()))
        self.dsets = self.variables
        self.dsets.keys()

    def setTimeIndex(self, Indx, update = None):
        self.currentTimeIndex = 0
        self.currentFileIndex = Indx
        print("File Index", self.currentFileIndex)
        if update is None:
            self.updateData()

    def getTime(self) :
    
        #Variable time is seconds past 1/1/2000 @ 12:00:00
        #Select the first index in the variable - scan start time
        tm = self.readNCVariable('time_bounds')
        self.timeObj  = dt.datetime(2000,1,1,12) + dt.timedelta(seconds=int(tm[0]))
        self.timeString = self.timeObj.strftime("%d %b %Y, %H:%M:%S UTC")
	
        return self.timeString

    def variableExist(self,varname):
        if varname in self.variableList :
            return True
        else:
            return False

    def doesExist(self, GridNo):
        if(self.fileList[GridNo]):
            return 1
        else:
            return None

    def setProjection(self,gid,axs=None,west=None,north=None,east=None,south=None):
        i = gid-1

        #If there is input on grid extent change the corners
        if (west != None and north != None and east != None and south != None):
            llcrnrlon = west
            llcrnrlat = south
            urcrnrlon = east
            urcrnrlat = north
        else:
            llcrnrlon = self.ll_lon[i]
            llcrnrlat = self.ll_lat[i]
            urcrnrlon = self.ur_lon[i]
            urcrnrlat = self.ur_lat[i]
	    
	#GOES is not in a specific projection so using cyl    
#        self.map[i] = Basemap(ax=axs,projection=self.projectionType,
#                                   llcrnrlon=self.ll_lon[i],llcrnrlat=self.ll_lat[i],
#                                   urcrnrlat = self.ur_lat[i],urcrnrlon = self.ur_lon[i],
#                                   resolution=self.resolution)
        self.map[i] = Basemap(ax=axs, projection=self.projectionType,
                      lat_ts = self.lat0[i], llcrnrlon = llcrnrlon,
                      llcrnrlat = llcrnrlat, urcrnrlon = urcrnrlon,
                      urcrnrlat = urcrnrlat, resolution=self.resolution)

        self.maplon[i],self.maplat[i] = self.map[i].makegrid(self.nx[i],self.ny[i])

    #This function navigate the input GOES file
    # i.e. determines the lat/long for each pixel
    # Based on the GOES-R Product Definition and User's Guide (PUG)
    def createLatLong(self):

       #Conversion
       rad2deg = 180./np.pi
        
       #Get input x and y fixed grid projection coordinates
       x1 = self.ncId.variables['x'][:] 
       y1 = self.ncId.variables['y'][:]
       #Get projection attributes
       req = self.ncId.variables['goes_imager_projection'].semi_major_axis
       rpol = self.ncId.variables['goes_imager_projection'].semi_minor_axis
       H = self.ncId.variables['goes_imager_projection'].perspective_point_height + req
       lon0 = self.ncId.variables['goes_imager_projection'].longitude_of_projection_origin
       
       #Create the projection coordinates meshgrid
       x, y = np.meshgrid(x1,y1)
      
       #Calculate variables
       a = (np.sin(x))**2 + (np.cos(x)**2)*((np.cos(y)**2) + (req**2/rpol**2)*(np.sin(y)**2))
       b = -2.*H*np.cos(x)*np.cos(y)
       c = H**2 - req**2
       rs = (-1.*b - np.sqrt(b**2 - 4*a*c))/(2*a)
       sx = rs*np.cos(x)*np.cos(y)
       sy = -1.*rs*np.sin(x)
       sz = rs*np.cos(x)*np.sin(y) 
     
       #Calculate lat/lon
       lat = np.arctan((req**2/rpol**2)*(sz/np.sqrt((H-sx)**2 + sy**2))) * rad2deg
       lon = lon0 - np.arctan(sy/(H-sx)) * rad2deg
     
       return lat, lon


    def setGridCorners(self):

        # Save current state.

        cg = self.currentGrid
        indx = self.currentTimeIndex

        for i in range(0,self.numGrids):

            if i == 0:
                self.setTimeIndex(0,update=False)
                self.setGrid(i+1)
            else:
                self.setGrid(i+1,update=False)
                self.setTimeIndex(0)

            #Get lat/long
            self.glats[i], self.glons[i] = self.createLatLong() 
 
            #Get geospatial information from the input file 
            self.lat0[i] = self.ncId.variables['geospatial_lat_lon_extent'].geospatial_lat_center
            self.lon0[i] = self.ncId.variables['geospatial_lat_lon_extent'].geospatial_lon_center
            
            self.ll_lon[i] = self.ncId.variables['geospatial_lat_lon_extent'].geospatial_westbound_longitude
            self.ur_lon[i] = self.ncId.variables['geospatial_lat_lon_extent'].geospatial_eastbound_longitude
            self.ll_lat[i] = self.ncId.variables['geospatial_lat_lon_extent'].geospatial_southbound_latitude
            self.ur_lat[i] = self.ncId.variables['geospatial_lat_lon_extent'].geospatial_northbound_latitude


            self.setProjection(i+1)

        self.setGrid(cg,update=False)
        self.setTimeIndex(indx)
    
    #Read in the radiance data and plot if required
    def getRadiance(self, plot = None):
        
        self.data = self.readNCVariable('Rad')
        #Check if plot is not None then check for connection to GUI for plotting
        if plot is not None:
            if self.pObj is not None:
                #Pass to plot object
                self.pObj.var = self.data
                varTitle = "Radiance (mW/(m$^2$ sr cm$^{-1}$))"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + "\nBand=" + str(self.band_num[self.currentGrid])
                varTitle = varTitle +', Wavelength=' + str(self.band_wave[self.currentGrid])
                self.pObj.varTitle = varTitle
                #Need to call plotting function to plot when variables are selected in the GUI
                # If statement is needed to prevent double plotting when time is changed
                if (self.varChange):
                    #Reset variable
                    self.varChange = False
                    self.pObj.pltFxn(self.pObj.pNum)

    #Calculate the brightness temperature
    def calculateTemperature(self):
        
        #Read in the constants from the netCDF file
        fk1 = self.ncId.variables['planck_fk1'][0]
        fk2 = self.ncId.variables['planck_fk2'][0]
        bc1 = self.ncId.variables['planck_bc1'][0]
        bc2 = self.ncId.variables['planck_bc2'][0]
 
        #Calculate brightness temperature [K]
        bt = (fk2 / (np.log((fk1/self.data)+1)) - bc1) / bc2

        #Pass to plot object
        if self.pObj is not None:
            self.pObj.var = bt
            varTitle = "Brightness Temperature (K)"
            varTitle = varTitle + "\n"+self.getTime()
            varTitle = varTitle + "\nBand=" + str(self.band_num[self.currentGrid])
            varTitle = varTitle +', Wavelength=' + str(self.band_wave[self.currentGrid])
            self.pObj.varTitle = varTitle
            print("GOES:", self.pObj.varTitle)
            #Need to call plotting function to plot when variables are selected in the GUI
            # If statement is needed to prevent double plotting when time is changed
            if (self.varChange):
                #Reset variable
                self.varChange = False
                self.pObj.pltFxn(self.pObj.pNum)
       
#####################  End of function setGridCorners() #######################
#####################  Start connection to GUI #######################

    def selectionChangeVar(self,i):
        self.varChange = True
        #Clear the color bar
        self.pObj.colormax = None
        self.pObj.colormin = None
        self.pObj.ncontours = None
        if self.pObj.colorbox is not None:
            self.pObj.colorbox.setParent(None)
            self.pObj.colorbox = None
        #Set dummy values for GUI plot object
        self.pObj.currentVar = i
        self.pObj.nz = 1
        #Get radiance
        self.getRadiance(plot = True)
        #Set variable for when time is changed
        self.raw = 1  
        self.derved = None
         
    def selectionChangeDVar(self,i):
        self.varChange = True
        #Clear the color bar
        self.pObj.colormax = None
        self.pObj.colormin = None
        self.pObj.ncontours = None
        if self.pObj.colorbox is not None:
            self.pObj.colorbox.setParent(None)
            self.pObj.colorbox = None
        self.pObj.currentVar = i
        self.pObj.nz = 1
        #Calculate raw counts
        self.getRadiance()
        self.calculateTemperature()
        #Set variable for when time is changed
        self.raw = None
        self.derived = 1

    def advanceTime(self,pobj):
        print("Advance Time")
        self.pObj = pobj
        if self.raw == 1:
            #Calculate raw counts
            self.getRadiance(plot = True)
        elif self.derived == 1:
            print("derived")
            #Calculate raw counts
            self.getRadiance()
            #Calculate derived variable
            self.calculateTemperature()

    def getDsetControlBar(self, plotObj):
        #Define GUI needed variables
        self.varChange = False
        self.ntimes = 1

        self.pObj = plotObj
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
                self.selectDset.addItem(dsetname)
                count += 1
        self.selectDset.setCurrentIndex(plotObj.currentDset-1)
        self.selectDset.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.selectDset.currentIndexChanged.connect(plotObj.selectionChangeDset)
        selectDsetWidgetLayout.addWidget(selectDsetLabel)
        selectDsetWidgetLayout.addWidget(self.selectDset)
        self.gboxLayout.addWidget(selectDsetWidget)

        # Create combo box for selecting variables
        varControlLabel = QLabel()
        varControlLabel.setText('Variable Control:')
        self.gboxLayout.addWidget(varControlLabel)
        plotObj.selectVar = QComboBox()
        plotObj.selectVar.setStyleSheet(Layout.QComboBox())
        plotObj.selectVar.addItems(self.rawvarlist)
        selectVarWidget = QWidget()
        selectVarLabel = QLabel()
        selectVarLabel.setText('Variable:')
        selectVarLayout = QHBoxLayout()
        selectVarLayout.addWidget(selectVarLabel)
        selectVarLayout.addWidget(plotObj.selectVar)
        selectVarWidget.setLayout(selectVarLayout)
        plotObj.selectVar.activated.connect(self.selectionChangeVar)

        #Derived Variables
        self.selectdVar = QComboBox()
        self.selectdVar.setStyleSheet(Layout.QComboBox())
        self.selectdVar.addItems(self.dvarlist[self.currentGrid])
        selectdVarWidget = QWidget()
        selectdVarLabel = QLabel()
        selectdVarLabel.setText('Derived Var:')
        selectdVarLayout = QHBoxLayout()
        selectdVarLayout.addWidget(selectdVarLabel)
        selectdVarLayout.addWidget(self.selectdVar)
        selectdVarWidget.setLayout(selectdVarLayout)
        self.selectdVar.activated.connect(self.selectionChangeDVar)


        self.gboxLayout.addWidget(selectVarWidget)
        self.gboxLayout.addWidget(selectdVarWidget)


        timeControlLabel = QLabel()
        timeControlLabel.setText('Time Control:')
        self.gboxLayout.addWidget(timeControlLabel)

        timebar = QWidget()
        timebarLayout = QHBoxLayout()
        timebar.setLayout(timebarLayout)
        timeWidgetLabel = QLabel()
        timeWidgetLabel.setText('Time:')
        plotObj.selectTime = QComboBox()
        plotObj.selectTime.setStyleSheet(Layout.QComboBox())
        plotObj.selectTime.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectTime.addItems(plotObj.dataSet.timeList[plotObj.dataSet.currentGrid])
        plotObj.selectTime.activated.connect(plotObj.selectionChangeTime)
        timebarLayout.addWidget(timeWidgetLabel)
        timebarLayout.addWidget(plotObj.selectTime)
        self.gboxLayout.addWidget(timebar)

        cpanel = QWidget()
        cpanelLayout = QHBoxLayout()
        cpanel.setLayout(cpanelLayout)
        nxtButton = QPushButton()
        nxtButton.setStyleSheet(Layout.QPushButton2())
        nxtButton.setText('&Next')
        nxtButton.setFixedWidth(75)
        nxtButtonLayout = QHBoxLayout()
        nxtButton.clicked.connect(plotObj.nxtButtonAction)
        prevButton = QPushButton()
        prevButton.setStyleSheet(Layout.QPushButton2())
        prevButton.setText('&Prev')
        prevButton.setFixedWidth(75)
        prevButtonLayout = QHBoxLayout()
        prevButton.clicked.connect(plotObj.prevButtonAction)
        cpanelLayout.addWidget(prevButton)
        cpanelLayout.addWidget(nxtButton)
        self.gboxLayout.addWidget(cpanel)
        plotObj.tabbingLayout.addWidget(self.gbox)
        plotObj.qscrollLayout.addWidget(self.tabbing,Qt.AlignTop)

        return self.qscroll

###############################################################################
####                          End WRFDataset() Object                      ####
###############################################################################

