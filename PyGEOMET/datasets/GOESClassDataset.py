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

class GOESClassDataset:

    # The constructor can be initialized by specifying the #
    # directory and file prefix for WRF output files.  If  #
    # the user does not provide this information when the  #
    # object is initalized, then set these variables to    #
    # None.                                                #

    def __init__(self,path = None, prefix = None):

        # Intialize a list to store list containing filenames #
        # for each domain.                                    #

        self.fileList = [None]

        self.timeList = [None]
        # Specify the WRF dataset by calling the function name()#

        self.name(path, prefix)

        # Current grid to which all operations apply !

        self.currentGrid = 1

        self.ncId = None

        self.variables = None

        self.attributes = None

        self.varName    = None

        self.ntimes = None

        self.currentTimeIndex = 0

        self.currentFileIndex = 0

        self.units='None'
	
        self.projectionType = 'cyl'
	
        self.pObj = None
          
        self.raw = None

        self.derived = None

        self.dvarlist = [None]
 
        self.dsetname = "GOES Class"

        self.resolution = "l"

        #Define plot type available for the dataset within the GUI
        self.ptypes = ['Horizontal Slice']

    # name() is the function for specifying the WRF dataset. #
    # Based on the directory name and prefix for WRF output  #
    # files, this function will compile a list of file names #
    # for each WRF domain. 

    def name(self, path, prefix) :

        # If no valid path and prefix are specified, then do #
        # nothing and set number of grids = 0                #

        self.path = path

        if path == None :

            self.numGrids = 0

        else :

            # If a valid path and prefix has been specified  #
            # then figure out the total number of files with #
            # the same prefix.                               #

            self.dsPrefix = os.path.join(path,prefix)
            files = sorted(glob.glob(self.dsPrefix+'*'))

            totalFiles = len(files)
            i = 1

            while(True):
                files = sorted(glob.glob(self.dsPrefix+'*BAND_0'+str(i)+'*'))

                numFiles = len(files)
                times = []
                append_time = times.append
                

                if numFiles != 0:
                    self.fileList.append(files)
                    self.setGrid(self.numGrids+1)
                    if i == 1:
                        dvar = ['Radiance', 'Effective Albedo']
                    else:
                        dvar = ['Radiance', 'Temperature']
                    self.dvarlist.append(dvar)
                    for ii in range(len(files)*self.ntimes):
                        self.setTimeIndex(ii)
                        append_time(self.getTime())
                    self.timeList.append(times)
                    self.numGrids += 1
                    self.setTimeIndex(0,update=False)
                totalFiles = totalFiles-numFiles
                i = i+1
                if totalFiles == 0:
                    break;
                #self.setTimeIndex(0,update=False)
            self.setTimeIndex(0,update=False)
            self.setGrid(1)
            self.variableReadFxn = { }
            for varname in self.variableList :
                self.variableReadFxn[varname] = self.readNCVariable
            self.setGridDefinition()
            self.rawvarlist = ['Raw Counts']

    def readNCVariable(self, vname) :
        variable = self.dsets[vname][self.currentTimeIndex]
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
        self.i_bad  = [0.0]*self.numGrids
        self.j_bad  = [0.0]*self.numGrids

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
            self.nx.append(self.ncId.dimensions['xc'].size)
            self.ny.append(self.ncId.dimensions['yc'].size)
            self.band_num.append(self.ncId.variables['bands'][:])
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
        return len(self.fileList[self.currentGrid])

    def getFileList(self):
        return(self.fileList[self.currentGrid])

    def updateData(self):	
        self.ncId = netCDF4.Dataset(self.fileList[self.currentGrid][self.currentFileIndex],'r')
        self.variables = self.ncId.variables
        self.attributes = self.ncId.__dict__
        self.ntimes = self.ncId.dimensions['time'].size
        self.variableList = list(sorted(self.variables.keys()))
        self.dsets = self.variables
        self.dsets.keys()

    def setTimeIndex(self, Indx, update = None):
   
        if self.ntimes == 1:
            self.currentTimeIndex = 0
            self.currentFileIndex = Indx
        else:
            self.currentFileIndex = int(Indx/self.ntimes)
            self.currentTimeIndex = Indx % self.ntimes
        if update is None:
            self.updateData()

    def getTime(self) :
    
        #Variable time is seconds past 1/1/1970
        tm = self.readNCVariable('time')
        self.timeObj  = dt.datetime(1970,1,1) + dt.timedelta(seconds=int(tm))
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
                      urcrnrlat = llcrnrlon, resolution=self.resolution)

        self.maplon[i],self.maplat[i] = self.map[i].makegrid(self.nx[i],self.ny[i])

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
                self.setTimeIndex(indx)

            #Read in lat/lon
            self.glons[i] = self.ncId.variables['lon'][:]
            self.glats[i] = self.ncId.variables['lat'][:]
	    
            #Determine the locations of bad lat/lon data
            self.j_bad[i], self.i_bad[i] = np.where((self.glons[i] > 180.) | 
                                                    (self.glons[i] < -180) |
                                                    (self.glats[i] > 90) | 
                                                    (self.glats[i] < -90))	    

            #Set bad values to NaN
            lon_new = self.glons[i]
            lat_new = self.glats[i]
            lon_new[self.j_bad[i],self.i_bad[i]] = np.nan
            lat_new[self.j_bad[i],self.i_bad[i]] = np.nan	    
             
            self.lat0[i] = np.nanmean(lat_new)
            self.lon0[i] = np.nanmean(lon_new) 
            
            self.ll_lon[i] = np.floor(np.nanmin(lon_new))
            self.ur_lon[i] = np.ceil(np.nanmax(lon_new))

            self.ll_lat[i] = np.floor(np.nanmin(lat_new))
            self.ur_lat[i] = np.ceil(np.nanmax(lat_new))

            self.setProjection(i+1)

        self.setGrid(cg,update=False)
        self.setTimeIndex(indx)

    def convertToCounts(self, plot = None, directPlot = None):
        
        self.data = self.readNCVariable('data')/32.
        self.data[self.j_bad[self.currentGrid-1],self.i_bad[self.currentGrid-1]] = np.nan
     
        #Check if plot is not None then check for connection to GUI for plotting
        if plot is not None:
            if self.pObj is not None:
                #Pass to plot object
                self.pObj.var = self.data
                varTitle = "Raw Counts (0-1023)"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])
                self.pObj.varTitle = varTitle
                #Needed because advance time plots in Plot Object but
                # plotting when changing variables is handled in this dataset
                if (directPlot is not None):
                    self.pObj.pltFxn(self.pObj.pNum)                


    def calculateVariables(self,vaname,directPlot = None):
        print(self.band_num[self.currentGrid],vaname)
        if self.band_num[self.currentGrid] == 1:
            #Radiance is needed for all calculations so calculate it first            
            m1 = 0.6120196  #[W/(m2 sr um count)]
            b1 = -17.749    #[W/(m2 sr um)]
            radiance = abs(m1*self.data + b1)  #[W/(m2 sr um)]
            if vaname == "Radiance":
                pltvar = radiance            
                varTitle = "Radiance (W/(m$^2$ sr $\mu$m))"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])
            elif vaname == "Effective Albedo":
                k = 1.86544e-3   #[(m2 sr um)/W]
                a_eff = k*radiance*100
                pltvar = a_eff
                varTitle = "Effective Albedo (%)"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])

        elif self.band_num[self.currentGrid] == 2:
            #Radiance is needed for all calculations so calculate it first
            m2 = 227.3889
            b2 = 68.2167
            radiance = (self.data - b2)/m2 * 10.  #[W/(m2 sr um)]
            if vaname == "Radiance":
                pltvar = radiance
                varTitle = "Radiance (W/(m$^2$ sr $\mu$m))"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])
            elif vaname == "Temperature":
                c1 = 1.19100e-5      #[mW/m2 sr cm-4]
                c2 = 1.438833        #[K/cm-1]
                v2a = 2591.74        #cm-1
                alpha2a = -1.437204
                beta2a = 1.002562
                #Effective temperature in (K)
                t_eff = (c2 * v2a) / np.log(1+ (c1*v2a**3)/(radiance/10.))

                #Convert to actual temperature (K)
                t_act = alpha2a + beta2a*t_eff

                pltvar = t_act
                varTitle = "Temperature (K)"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])

        elif self.band_num[self.currentGrid] == 3:
            #Radiance is needed for all calculations so calculate it first
            m3 = 38.8383
            b3 = 29.1287
            radiance = (self.data - b3)/m3 * 10.  #[W/(m2 sr um)]
            if vaname == "Radiance":
                pltvar = radiance
                varTitle = "Radiance (W/(m$^2$ sr $\mu$m))"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])
            elif vaname == "Temperature":
                c1 = 1.19100e-5      #[mW/m2 sr cm-4]
                c2 = 1.438833        #[K/cm-1]
                v3a = 1522.52        #cm-1
                alpha3a = -3.607841
                beta3a = 1.0010018
                #Effective temperature in (K)
                t_eff = (c2 * v3a) / np.log(1+ (c1*v3a**3)/(radiance/10.))

                #Convert to actual temperature (K)
                t_act = alpha3a + beta3a*t_eff

                pltvar = t_act
                varTitle = "Temperature (K)"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])

        elif self.band_num[self.currentGrid] == 4:
            #Radiance is needed for all calculations so calculate it first
            m4 = 5.2285     
            b4 = 15.6854    
            radiance = (self.data - b4)/m4 * 10.  #[W/(m2 sr um)]
            if vaname == "Radiance":
                pltvar = radiance
                varTitle = "Radiance (W/(m$^2$ sr $\mu$m))"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])
            elif vaname == "Temperature":
                c1 = 1.19100e-5      #[mW/m2 sr cm-4]
                c2 = 1.438833        #[K/cm-1]
                v4a = 937.23         #cm-1
                alpha4a = -0.386043
                beta4a = 1.001298
                #Effective temperature in (K)
                t_eff = (c2 * v4a) / np.log(1+ (c1*v4a**3)/(radiance/10.))

                #Convert to actual temperature (K)
                t_act = alpha4a + beta4a*t_eff

                pltvar = t_act 
                varTitle = "Temperature (K)"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])

        elif self.band_num[self.currentGrid] == 6:
            #Radiance is needed for all calculations so calculate it first
            m6 = 5.5297
            b6 = 16.5892
            radiance = (self.data - b6)/m6 * 10.  #[W/(m2 sr um)]
            if vaname == "Radiance":
                pltvar = radiance
                varTitle = "Radiance (W/(m$^2$ sr $\mu$m))"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])
            elif vaname == "Temperature":
                c1 = 1.19100e-5      #[mW/m2 sr cm-4]
                c2 = 1.438833        #[K/cm-1]
                v6a = 751.93         #cm-1
                alpha6a = -0.134688
                beta6a = 1.000481
                #Effective temperature in (K)
                t_eff = (c2 * v6a) / np.log(1+ (c1*v6a**3)/(radiance/10.))

                #Convert to actual temperature (K)
                t_act = alpha6a + beta6a*t_eff

                pltvar = t_act
                varTitle = "Temperature (K)"
                varTitle = varTitle + "\n"+self.getTime()
                varTitle = varTitle + ', Band=' + str(self.band_num[self.currentGrid])

        #Pass to plot object
        if self.pObj is not None:
            self.pObj.var = pltvar
            self.pObj.varTitle = varTitle
            #Needed because advance time plots in Plot Object but
            # plotting when changing variables is handled in this dataset
            if directPlot is not None:
                self.pObj.pltFxn(self.pObj.pNum)

        

#####################  End of function setGridCorners() #######################
#####################  Start connection to GUI #######################

    def selectionChangeVar(self,i):
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
        #Calculate raw counts
        self.convertToCounts(plot = True, directPlot = True)
        #Set variable for when time is changed
        self.raw = 1  
        self.derved = None
         
    def selectionChangeDVar(self,i):
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
        self.convertToCounts()
        self.calculateVariables(self.dvarlist[self.currentGrid][i], directPlot = True)
        #Set variable for when time is changed
        self.raw = None
        self.derived = 1

    def selectionChangeTime(self,i):
        self.pObj.currentTime = i
        self.setTimeIndex(self.pObj.currentTime)
        if self.raw == 1:
            #Calculate raw counts
            self.convertToCounts(plot = True)
        elif self.derived == 1:
             #Calculate raw counts
             self.convertToCounts()
             #Calculate derived variable
             self.calculateVariables(self.dvarlist[self.currentGrid][self.pObj.currentVar])

    def selectionChangeBand(self,i):
        #Clear the color bar
        self.pObj.colormax = None
        self.pObj.colormin = None
        self.pObj.ncontours = None
        if self.pObj.colorbox is not None:
            self.pObj.colorbox.setParent(None)
            self.pObj.colorbox = None
        #Clear the derived variable list
        self.selectdVar.clear()
        self.pObj.currentGrid = i+1
        self.pObj.currentTime = 0
        self.setTimeIndex(self.pObj.currentTime)
        self.setGrid(self.pObj.currentGrid)
        #Clear the time list and the add current grid list
        self.pObj.selectTime.clear()
        self.pObj.selectTime.addItems(self.timeList[self.currentGrid])
        #Add the derived variables based on the selected band
        if self.band_num[self.currentGrid] == 1:
            self.selectdVar.addItems(self.dvarlist[self.currentGrid])
        else:
            self.selectdVar.addItems(self.dvarlist[self.currentGrid])
        if self.raw == 1:
            self.selectionChangeVar(1)
        elif self.derived == 1:
            self.selectionChangeDVar(self.pObj.currentVar)
        else:
            self.pObj.pltFxn(self.pObj.pNum)        

    def advanceTime(self,pobj):
         self.pObj = pobj
         if self.raw == 1:
             #Calculate raw counts
             self.convertToCounts(plot = True)
         elif self.derived == 1:
             #Calculate raw counts
             self.convertToCounts()
             #Calculate derived variable
             self.calculateVariables(self.dvarlist[self.currentGrid][self.pObj.currentVar])
#         else:
#             self.pObj.pltFxn(self.pObj.pNum)

    def nxtButtonAction(self):
         self.pObj.currentTime+=1
         if self.pObj.currentTime == self.ntimes*self.getNumFiles():
             self.pObj.currentTime = 0
         self.pObj.selectTime.setCurrentIndex(self.pObj.currentTime)
         self.setTimeIndex(self.pObj.currentTime)
         if self.raw == 1:
             #Calculate raw counts
             self.convertToCounts(plot = True)
         elif self.derived == 1:
             #Calculate raw counts
             self.convertToCounts()
             #Calculate derived variable
             self.calculateVariables(self.dvarlist[self.currentGrid][self.pObj.currentVar])
         else:
             self.pObj.pltFxn(self.pObj.pNum)

    def prevButtonAction(self):
         self.pObj.currentTime-=1
         if self.pObj.currentTime == -1:
             self.pObj.currentTime = self.getNumFiles()*self.ntimes-1
         self.pObj.selectTime.setCurrentIndex(self.pObj.currentTime)
         self.setTimeIndex(self.pObj.currentTime)
         if self.raw == 1:
             #Calculate raw counts
             self.convertToCounts(plot = True)
         elif self.derived == 1:
             #Calculate raw counts
             self.convertToCounts()
             #Calculate derived variable
             self.calculateVariables(self.dvarlist[self.currentGrid][self.pObj.currentVar])
         else:
             self.pObj.pltFxn(self.pObj.pNum)

    def getDsetControlBar(self, plotObj):
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

        #Band selection
        self.selectGrid = QComboBox()
        self.selectGrid.setStyleSheet(Layout.QComboBox())
        self.selectGrid.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        for i in range(0,plotObj.dataSet.numGrids):
            self.selectGrid.addItem(str(self.band_num[i+1]))
        self.selectGridWidget = QWidget()
        selectGridLabel = QLabel()
        selectGridLabel.setText('Select Band:')
        selectGridLayout = QVBoxLayout()
        selectGridLayout.addWidget(selectGridLabel)
        selectGridLayout.addWidget(self.selectGrid)
        self.selectGridWidget.setLayout(selectGridLayout)
        self.selectGrid.currentIndexChanged.connect(self.selectionChangeBand)
        self.gboxLayout.addWidget(self.selectGridWidget)
        
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
        plotObj.selectTime.activated.connect(self.pObj.selectionChangeTime)
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
        nxtButton.clicked.connect(self.pObj.nxtButtonAction)
        prevButton = QPushButton()
        prevButton.setStyleSheet(Layout.QPushButton2())
        prevButton.setText('&Prev')
        prevButton.setFixedWidth(75)
        prevButtonLayout = QHBoxLayout()
        prevButton.clicked.connect(self.pObj.prevButtonAction)
        cpanelLayout.addWidget(prevButton)
        cpanelLayout.addWidget(nxtButton)
        self.gboxLayout.addWidget(cpanel)
        plotObj.tabbingLayout.addWidget(self.gbox)
        plotObj.qscrollLayout.addWidget(self.tabbing,Qt.AlignTop)

        return self.qscroll

###############################################################################
####                          End WRFDataset() Object                      ####
###############################################################################

