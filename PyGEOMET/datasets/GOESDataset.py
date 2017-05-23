import numpy as np
import os
import glob
import netCDF4
from datetime import datetime
from mpl_toolkits.basemap import Basemap
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import PyGEOMET.utils.LayoutFormat as Layout

class GOESDataset:

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

        self.dsetname = "GOES UAH"
 
        self.dvarlist = []

        self.resolution = "l"

        #Define plot type available for the dataset within the GUI
        self.ptypes = ['Horizontal Slice', 'Time Series']

    # name() is the function for specifying the GOES dataset. #
    # Based on the directory name and prefix for WRF output  #
    # files, this function will compile a list of file names #
    # for each WRF domain.                                   #


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
            files = sorted(glob.glob(self.dsPrefix+'*'),reverse=True)
            for ind in files:
                self.fileList.append(ind)
            totalFiles = len(self.fileList)
            times = []            
            append_time = times.append
            
            for i in range(totalFiles-1):

                self.setGrid(i+1)
                for ii in range(self.ntimes):
                    self.setTimeIndex(ii)
                    append_time(self.getTime())
                self.timeList.append(times)
                self.numGrids += 1

            self.setTimeIndex(0,update=False)
            self.setGrid(1)
            self.variableReadFxn = { }
            for varname in self.variableList :
                self.variableReadFxn[varname] = self.readNCVariable
            self.setGridDefinition()

    def readNCVariable(self, vname,barbs=None,vectors=None,contour2=None,varonly=False) :
        variable = self.dsets[vname][:]
        dims = variable.shape
        if dims[0] == self.ntimes:
            variable = variable[self.currentTimeIndex]
        if( hasattr(self.dsets[vname],'description')):
            self.description = self.dsets[vname].description
        if( hasattr(self.dsets[vname],'units')):
            self.units = self.dsets[vname].units
        if( hasattr(self.dsets[vname],'missing_value')):
            missing = self.dsets[vname].missing_value
            variable[variable == missing] = np.nan
        return variable

    def getAttribute(attname):
        return self.ncId.__dict__[attname]

    def setGridDefinition(self):

        self.nx  = []
        self.ny  = []
        self.dx  = []
        self.dy  = []
        self.lat1 = []
        self.lat2 = []
        self.lon0 = []
        self.lat0 = []
        self.map   = [None]*self.numGrids
        self.wd     = [0.0]*self.numGrids
        self.ht     = [0.0]*self.numGrids
        self.padset    = [-1]*self.numGrids
        self.pad    = [(0.0,0.0,0.0,0.0)]*self.numGrids
        self.ur_lon = [0.0]*self.numGrids
        self.ur_lat = [0.0]*self.numGrids
        self.ll_lon = [0.0]*self.numGrids
        self.ll_lat = [0.0]*self.numGrids
        self.xs = [None]*self.numGrids
        self.ys = [None]*self.numGrids
        self.glats  = [None]*self.numGrids
        self.glons  = [None]*self.numGrids

        cg = self.currentGrid
        indx = self.currentFileIndex
        for i in range(0,self.numGrids):
            self.setTimeIndex(0,update=False)
            self.setGrid(i+1)

            self.nx.append(self.ncId.dimensions['west_east'].size)
            self.ny.append(self.ncId.dimensions['south_north'].size)
            self.dx.append(self.ncId.__dict__['DX'])
            self.dy.append(self.ncId.__dict__['DY'])
            self.lat1.append(self.ncId.__dict__['TRUELAT1'])
            self.lat2.append(self.ncId.__dict__['TRUELAT2'])
            self.lon0.append(self.ncId.__dict__['CEN_LON'])
            self.lat0 .append(self.ncId.__dict__['CEN_LAT'])

        if len(self.nx) == self.numGrids:
            self.currentGrid = cg
            self.currentFileIndex = indx

        proj_opt = ['lcc','npstere','merc']
        self.projectionType = proj_opt[self.ncId.__dict__['MAP_PROJ']-1]
        self.setGridCorners()

    # getNumGrids() is a function that returns the maximum     #
    # number of grids detected by the object for the specified #
    # directory and file prefix. 

    def setGrid(self, GridNo, update = None):
        self.currentGrid = GridNo
        if update is None:
            self.updateData()

    def getNumFiles(self):
        return len(self.fileList[self.currentGrid])

    def getFileList(self):
        return(self.fileList[self.currentGrid])

    def updateData(self):
        self.ncId = netCDF4.Dataset(self.fileList[self.currentGrid],'r')
        self.variables = self.ncId.variables
        self.attributes = self.ncId.__dict__
        self.ntimes = self.ncId.dimensions['Time'].size
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
        tm = self.readNCVariable('Times')
        tstr   = tm[0]
        year   = tm[0]+tm[1]+tm[2]+tm[3]
        month  = tm[5]+tm[6]
        day    = tm[8]+tm[9]
        hour   = tm[11]+tm[12]
        minute = tm[14]+tm[15]
        second = tm[17]+tm[18]

        self.timeObj  = datetime(
                                    int(year),
                                    int(month),
                                    int(day),
                                    int(hour),
                                    int(minute),
                                    int(second)
                                    )
        self.timeString = self.timeObj.strftime("%d %b %Y, %H:%M:%S UTC")
        return self.timeString

    def setMapPad(self,pad, grid=0, ang = True) :
        if(int(grid)) :
            if(ang):
                self.padset[grid-1] = 1
            else:
                self.padset[grid-1] = 2
            self.pad[grid-1] = pad
        else:
            if(ang):
                self.padset[grid-1] = [1]*self.numGrids
            else:
                self.padset[grid-1] = [2]*self.numGrids
            self.pad[:] = [pad]*self.numGrids
        self.setGridCorners()


    def setProjection(self,gid,axs=None,i0=None,i1=None,j0=None, j1=None):
        i = gid-1
        if(i0 != None and j0 != 0):
            self.ll_lon[i] = self.glons[i][i0,j0]
            self.ll_lat[i] = self.glats[i][i0,j0]

        if(i1 != None and j1 != 0):
            self.ur_lon[i] = self.glons[i][i1,j1]
            self.ur_lat[i] = self.glats[i][i1,j1]

        if self.projectionType == 'lcc':
            self.map[i] = Basemap(ax=axs,projection=self.projectionType,lon_0=self.lon0[i],
                                 lat_0 =self.lat0[i], lat_1=self.lat1[i], lat_2=self.lat2[i],
                                 width=self.wd[i],height=self.ht[i],resolution=self.resolution)
        if self.projectionType == 'npstere':
            self.map[i] = Basemap(ax=axs,projection=self.projectionType,
                                       lat_0=self.lat0[i], lon_0=self.lon0[i],
                                       llcrnrlon=self.ll_lon[i],llcrnrlat=self.ll_lat[i],
                                       urcrnrlat = self.ur_lat[i],urcrnrlon = self.ur_lon[i],
                                       resolution=self.resolution)
        if self.projectionType == 'merc':
            self.map[i] = Basemap(ax=axs,projection=self.projectionType,
                                       lat_0=self.lat0[i], lon_0=self.lon0[i],
                                       llcrnrlon=self.ll_lon[i],llcrnrlat=self.ll_lat[i],
                                       urcrnrlat = self.ur_lat[i],urcrnrlon = self.ur_lon[i],
                                       resolution=self.resolution)

        x_ll,y_ll = self.map[0](self.ll_lon[i],self.ll_lat[i])

        # Calculate the lower left coordinates of the outer most grid
        x_ur = x_ll + self.wd[i]
        y_ur = y_ll + self.ht[i]

        if(self.padset[i] > 0) :
            p1,p2,p3,p4 = self.pad[i]

            if(self.padset[i] == 2) :
                ll_lon1,ll_lat1 = self.map[0](x_ll-p1,y_ll-p2,inverse=True)
                ur_lon1,ur_lat1 = self.map[0](x_ur+p3,y_ur+p4,inverse=True)
            else:
                ll_lon1 = self.ll_lon[i]-p1
                ll_lat1 = self.ll_lat[i]-p2
                ur_lon1 = self.ur_lon[i]+p3
                ur_lat1 = self.ur_lat[i]+p4

            self.map[i] = Basemap(ax=axs,projection=self.projectionType,lat_ts=self.lat1,
                                    lat_0=self.lat0, lon_0=self.lon0,
                                    llcrnrlon=ll_lon1,llcrnrlat=ll_lat1,
                                    urcrnrlat = ur_lat1,urcrnrlon = ur_lon1
                                  )
            x_ll,y_ll = self.map[0](self.ll_lon[i],self.ll_lat[i])
            x_ur = x_ll + self.wd[i]
            y_ur = y_ll + self.ht[i]

        self.xs[i] = (x_ll, x_ur, x_ur, x_ll, x_ll)
        self.ys[i] = (y_ll, y_ll, y_ur, y_ur, y_ll)

    def setGridCorners(self):

        # Save current state.

        cg = self.currentGrid
        indx = self.currentTimeIndex

        for i in range(0,self.numGrids):

            self.setTimeIndex(0,update=False)
            self.setGrid(i+1)

            self.glons[i] = self.readNCVariable('XLONG_M')
            dims = self.glons[i].shape
            self.ll_lon[i] = self.glons[i][0,0]
            self.ur_lon[i] = self.glons[i][dims[0]-1,dims[1]-1]

            self.glats[i] = self.readNCVariable('XLAT_M')
            dims = self.glats[i].shape
            self.ll_lat[i] = self.glats[i][0,0]
            self.ur_lat[i] = self.glats[i][dims[0]-1,dims[1]-1]

            self.wd[i] = (self.nx[i]-2)*self.dx[i]
            self.ht[i] = (self.ny[i]-2)*self.dy[i]
            self.setProjection(i+1)

        self.setGrid(cg,update=False)
        self.setTimeIndex(indx)


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
                dsetname = os.path.basename(str(plotObj.dSet[count].path))
                self.selectDset.addItem(dsetname)
                count += 1
        self.selectDset.setCurrentIndex(plotObj.currentDset-1)
        self.selectDset.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.selectDset.currentIndexChanged.connect(plotObj.selectionChangeDset)
        selectDsetWidgetLayout.addWidget(selectDsetLabel)
        selectDsetWidgetLayout.addWidget(self.selectDset)
        self.gboxLayout.addWidget(selectDsetWidget)

        # Create combo box for selecting model variables
        varControlLabel = QLabel()
        varControlLabel.setText('Variable Control:')
        self.gboxLayout.addWidget(varControlLabel)

        plotObj.selectVar = QComboBox()
        plotObj.selectVar.setStyleSheet(Layout.QComboBox())
        plotObj.selectVar.addItems(plotObj.dataSet.variableList)
        selectVarWidget = QWidget()
        selectVarLabel = QLabel()
        selectVarLabel.setText('Variable:')
        selectVarLayout = QHBoxLayout()
        selectVarLayout.addWidget(selectVarLabel)
        selectVarLayout.addWidget(plotObj.selectVar)
        selectVarWidget.setLayout(selectVarLayout)
        plotObj.selectVar.activated.connect(plotObj.selectionChangeVar)

        #Grid selection
        self.selectGrid = QComboBox()
        self.selectGrid.setStyleSheet(Layout.QComboBox())
        self.selectGrid.setSizeAdjustPolicy(QComboBox.AdjustToContents) 
        for i in range(0,plotObj.dataSet.numGrids):
            self.selectGrid.addItem(str(i+1))
        self.selectGridWidget = QWidget()
        selectGridLabel = QLabel()
        selectGridLabel.setText('Grid:')
        selectGridLayout = QHBoxLayout()
        selectGridLayout.addWidget(selectGridLabel)
        selectGridLayout.addWidget(self.selectGrid)

        #Level selection
        selectLevelLabel = QLabel()
        selectLevelLabel.setText('Level:')

        plotObj.selectLevel = QComboBox()
        plotObj.selectLevel.setStyleSheet(Layout.QComboBox())
        plotObj.selectLevel.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        plotObj.selectLevel.addItem('NA')
        selectGridLayout.addWidget(selectLevelLabel)
        selectGridLayout.addWidget(plotObj.selectLevel)
        self.selectGridWidget.setLayout(selectGridLayout)

        self.gboxLayout.addWidget(selectVarWidget)
        self.gboxLayout.addWidget(self.selectGridWidget)
        self.selectGrid.currentIndexChanged.connect(plotObj.selectionChangeGrid)
        plotObj.selectLevel.currentIndexChanged.connect(plotObj.selectionChangeLevel)

        #Time Control
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
        plotObj.selectTime.currentIndexChanged.connect(plotObj.selectionChangeTime)
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

