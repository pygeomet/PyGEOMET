import numpy as np
import os
import glob
import netCDF4
from datetime import datetime,timedelta
from mpl_toolkits.basemap import Basemap
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import PyGEOMET.utils.LayoutFormat as Layout

class CmaqDataset:

    def __init__(self,path = None, prefix = None):

        # Intialize a list to store list containing filenames #
        # for each domain.                                    #

        self.fileList = []

        self.timeList = [None]
        # Specify the CMAQ dataset by calling the function name()#

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

        self.dsetname = "CMAQ"

        self.resolution = "l"

        self.pftypes = sorted(['ACONC','CONC','B3GTS_S','AERODIAM','AEROVIS',
                               'CGRID','DRYDEP','BBB','PHOTDIAG1','PHOTDIAG2',
                               'SSEMIS','WETDEP1'])
      
        self.ftypes = []
 
        self.currentFType = 0
 
        self.varList = []

        #Define plot type available for the dataset within the GUI
        #Should make this file type dependent for CMAQ
        self.ptypes = ['Horizontal Slice', 'Vertical Slice', 'Vertical Profile', 
                       'Time Series', 'Difference Plot']
                       #, 'Hovmoller Diagram']

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

            #Not accounting for multiple CMAQ grids
            self.numGrids = 1
            self.dsPrefix = os.path.join(path,prefix)
            totalFiles = 0
            for ftype in self.pftypes:
                files = sorted(glob.glob(self.dsPrefix+'*.'+ftype+'*')) 
                totalFiles = totalFiles + len(files)
                if len(files) != 0:
                    self.ftypes.append(ftype)
                    self.fileList.append([None])
            i = 1
            while(True):
                files = sorted(glob.glob(self.dsPrefix+'*.'+self.ftypes[i-1]+'*'))
                print(files)
                numFiles = len(files)
                times = []
                append_time = times.append
                if numFiles != 0:
                    self.fileList[i-1].append(files)
                    self.setGrid(self.numGrids)
                    for ii in range(len(files)*self.ntimes):
                        self.setTimeIndex(ii)
                        append_time(self.getTime())
                    self.timeList.append(times)
                    self.setTimeIndex(0)
                    #self.setGrid(1)
                    self.threeDVars = []
                    self.twoDVars = []                    
                    self.variableReadFxn = { }
                    for varname in self.variableList:
                        self.variableReadFxn[varname] = self.readNCVariable
                        var = self.readNCVariable(varname,varonly=True)
                        if len(var.shape) > 2:
                            self.threeDVars.append(varname)
                        else:
                            self.twoDVars.append(varname)
                    i += 1
                    self.currentFType = i-1
                #self.numGrids += 1
                totalFiles = totalFiles-numFiles
                if totalFiles == 0:
                    break;
                self.setTimeIndex(0,update=False)
            
            self.currentFType = 0
            self.setGridDefinition()
                    
            #self.setTimeIndex(0,update=False)
            #self.setGrid(1)
            #self.threeDVars = []
            #self.twoDVars = []
            #self.variableReadFxn = { }
            #for varname in self.variableList:
            #    self.variableReadFxn[varname] = self.readNCVariable
            #    var = self.readNCVariable(varname,varonly=True)
            #    if len(var.shape) > 2:
            #        self.threeDVars.append(varname)
            #    else:
            #        self.twoDVars.append(varname)
            #self.setGridDefinition()
            #self.dvarlist = ['300mb_winds', '500mb_hgt', '500mb_temp']

    def readNCVariable(self,vname,barbs=None, vectors=None, contour2=None,varonly=False):
        variable = self.dsets[vname][self.currentTimeIndex]
        if (variable.shape[0] == 1):
            variable = variable[0]
        if varonly == False:
            if( hasattr(self.dsets[vname],'var_desc')):
                self.description = self.dsets[vname].var_desc
            if( hasattr(self.dsets[vname],'units')):
                self.units = self.dsets[vname].units
            if vname != 'Times':
                if len(self.dsets[vname].dimensions) == 4:
                    if len(self.dsets[vname].shape) >= 3:
                        tmp = np.arange(1,self.nz[self.currentGrid-1]+1,1)
                        self.levelList = ["%2d" % x for x in tmp]
                    else:
                        tmp = [1]
                        self.levelList = ["%2d" % x for x in tmp]
                else:
                    self.levelList = ['NA']

        return variable

    def setGridDefinition(self):

        #self.gridRatio = []
        #self.pID = []
        self.nx  = []
        self.ny  = []
        self.nz  = []
        self.dx  = []
        self.dy  = []
        #self.istart = []
        #self.jstart = []
        self.lat1 = []
        self.lat2 = []
        self.lon0 = []
        self.lat0 = []
        self.x_offset = []
        self.y_offset = []
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
        print(self.numGrids)
        for i in range(0,self.numGrids):
            if i == 0:
                self.setTimeIndex(0,update=False)
                self.setGrid(i+1)
            else:
                self.setGrid(i+1,update=False)
                self.setTimeIndex(indx)

            self.setGrid(i+1)
            self.setTimeIndex(0)
            #self.gridRatio.append(self.ncId.__dict__['PARENT_GRID_RATIO'])
            #self.pID.append(self.ncId.__dict__['PARENT_ID'])
            self.nx.append(self.ncId.__dict__['NCOLS'])
            self.ny.append(self.ncId.__dict__['NROWS'])
            self.nz.append(self.ncId.__dict__['NLAYS'])
            self.dx.append(self.ncId.__dict__['XCELL'])
            self.dy.append(self.ncId.__dict__['YCELL'])
            #self.istart.append(self.ncId.__dict__['I_PARENT_START'])
            #self.jstart.append(self.ncId.__dict__['J_PARENT_START'])
            self.lat1.append(self.ncId.__dict__['P_ALP'])
            self.lat2.append(self.ncId.__dict__['P_BET'])
            self.x_offset.append(self.ncId.XORIG)
            self.y_offset.append(self.ncId.YORIG)            

            #Calculate the CMAQ grid center
            temp_lon0 = self.ncId.XCENT
            temp_lat0 = self.ncId.YCENT
            self.wd[i] = self.ncId.XCELL*(self.ncId.NCOLS-1)
            self.ht[i] = self.ncId.YCELL*(self.ncId.NROWS-1)
            proj_opt = ['npstere','lcc','merc']
            pType = proj_opt[self.ncId.GDTYP-1]
            m_temp = Basemap(width=self.wd[i], height=self.ht[i], resolution = 'c',
                             projection=pType,
                             lat_0=temp_lat0, lon_0=temp_lon0,
                             lat_1=self.lat1[i], lat_2=self.lat2[i])
            
            xx = self.x_offset[i]+self.wd[i]
            yy = self.y_offset[i]+self.ht[i]
            clon,clat = m_temp(xx,yy,inverse=True)
            self.lon0.append(clon)
            print(clat)
            self.lat0.append(clat)

        if len(self.nx) == self.numGrids:
            self.currentGrid = cg
            self.currentFileIndex = indx

        proj_opt = ['npstere','lcc','merc']
        self.projectionType = proj_opt[self.ncId.__dict__['GDTYP']-1]
        self.setGridCorners()

    def setGrid(self, GridNo, update = None):
        self.currentGrid = GridNo
        if update is None:
            self.updateData()

    def getNumFiles(self):
        return len(self.fileList[self.currentFType][self.currentGrid])

    def getFileList(self):
        return(self.fileList[self.currentFType][self.currentGrid])

    def updateData(self):
        #print("FileList:", self.fileList)
        #print("Current FType:", self.currentFType)
        #print("Current Grid:", self.currentGrid)
        #print("Current File Index:", self.currentFileIndex)
        self.ncId = netCDF4.Dataset(self.fileList[self.currentFType][self.currentGrid][self.currentFileIndex],'r')
        self.variables = self.ncId.variables
        self.attributes = self.ncId.__dict__
        self.ntimes = self.ncId.dimensions['TSTEP'].size
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
        tm = self.readNCVariable('TFLAG')[0,:]      
        #Fill with zeros to make sure input format is yyyyddd and hhmmss
        #CMAQ will leave off 0s because the time is stored as integers instead of strings
        yyyyddd = str(tm[0]).zfill(7)
        hhmmss = str(tm[1]).zfill(6) 
        
        year   = yyyyddd[0:4]
        #Days is in Julian day
        days   = yyyyddd[4:7]
        hour   = hhmmss[0:2] 
        minute = hhmmss[2:4]
        second = hhmmss[4:6]
       
        #Set initial day to the 1st of the month then add Julian days
        self.timeObj  = (datetime(int(year),1,1,int(hour),int(minute),int(second))+
                         timedelta(int(days)-1))
        self.timeString = self.timeObj.strftime("%d %b %Y, %H:%M:%S UTC")
        return self.timeString

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
        self.test = axs

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

            #Create lat and lon arrays
            x = np.linspace(0,self.dx[i]*self.nx[i],self.nx[i])
            y = np.linspace(0,self.dy[i]*self.ny[i],self.ny[i])
            cx, cy = np.meshgrid(x,y)
            
            print( self.lon0[i], self.lat0[i])
            m_temp = Basemap(projection='lcc',lon_0=self.lon0[i],
                             lat_0 =self.lat0[i], lat_1=self.lat1[i], lat_2=self.lat2[i],
                             width=self.wd[i],height=self.ht[i],resolution='c')

            self.glons[i],self.glats[i] = m_temp(cx,cy,inverse=True)

            dims = self.glons[i].shape
            self.ll_lon[i] = self.glons[i][0,0]
            self.ur_lon[i] = self.glons[i][dims[0]-1,dims[1]-1]

            dims = self.glats[i].shape
            self.ll_lat[i] = self.glats[i][0,0]
            self.ur_lat[i] = self.glats[i][dims[0]-1,dims[1]-1]

            self.setProjection(i+1)

        self.setGrid(cg,update=False)
        self.setTimeIndex(indx)

#####################  End of function setGridCorners() #######################
#####################  Start connection to GUI #######################

    def selectionChangeFType(self,i):     
        self.currentFType = i
        if self.pObj.colorbox is not None:
            self.pObj.colorbox.setParent(None)
            self.pObj.colorbox = None
        if self.pObj.vslicebox is not None:
            self.pObj.vslicebox.setParent(None)
            self.pObj.vslicebox = None
        self.pObj.appobj.cbar.deleteLater()
        self.pObj.appobj.cbar = None
        self.pObj.appobj.cs = None
        self.pObj.appobj.cs2 = None
        self.pObj.appobj.barbs = None
        self.pObj.appobj.vectors = None
        self.pObj.appobj.vectorkey = None
        self.pObj.appobj.cs2label = None
        self.pObj.ColorBar = None
        self.pObj.appobj.domain_average = None
        #self.pObj.derivedVar = False
        self.pObj.appobj.recallProjection = True
        #Plus 1 is needed because None is technically the first index
        #self.dataSet = self.dSet[self.currentDset]
        self.pObj.getControlBar()
        self.selectFType.setCurrentIndex(self.currentFType)
        #self.currentGrid = 1
        self.currentTime = 0
        self.setTimeIndex(self.currentTime)
        #self.dataSet.setGrid(self.currentGrid, update=None)
        #Make sure all variables are none type
        self.pObj.currentVar = None
        self.pObj.currentdVar = None
        #self.currentdVar = None
        self.pObj.selectVar.clear()
        self.pObj.selectVar.addItems(self.pObj.dataSet.variableList)
        #self.pObj.appobj.axes1.remove(self.pObj.appobj.axes1[self.pObj.pNum-1])
        self.pObj.appobj.axes1[self.pObj.pNum-1] = None
        self.pObj.figure.clear()
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
                self.selectDset.addItem(dsetname+' [Dataset '+str(count)+']')
                count += 1
        self.selectDset.setCurrentIndex(plotObj.currentDset-1)
        self.selectDset.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        #self.selectDset.setMaximumWidth(plotObj.appobj.screenx*.15*.8)
        self.selectDset.currentIndexChanged.connect(plotObj.selectionChangeDset)
        selectDsetWidgetLayout.addWidget(selectDsetLabel)
        selectDsetWidgetLayout.addWidget(self.selectDset)
        self.gboxLayout.addWidget(selectDsetWidget)

        #File type selection
        selectFTypeWidget = QWidget()
        selectFTypeWidgetLayout = QVBoxLayout()
        selectFTypeWidget.setLayout(selectFTypeWidgetLayout)
        selectFTypeLabel = QLabel('File Type:')
        self.selectFType = QComboBox()
        self.selectFType.setStyleSheet(Layout.QComboBox())
        #self.ftypes = ['ACONC','CONC','B3GTS_S']
        self.selectFType.addItems(self.ftypes)
        self.selectFType.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.selectFType.activated.connect(self.selectionChangeFType)
        self.selectFType.setMaximumWidth(plotObj.appobj.screenx*.15*.8)
        selectFTypeWidgetLayout.addWidget(selectFTypeLabel)
        selectFTypeWidgetLayout.addWidget(self.selectFType)
        self.gboxLayout.addWidget(selectFTypeWidget)

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
        self.selectPlotType.setMinimumContentsLength(2)
        selectPlotWidgetLayout.addWidget(selectPlotLabel)
        selectPlotWidgetLayout.addWidget(self.selectPlotType)
        self.gboxLayout.addWidget(selectPlotWidget)

        # Create combo box for selecting model variables
        varControlLabel = QLabel()
        varControlLabel.setText('Variable Control:')
        self.gboxLayout.addWidget(varControlLabel)

        plotObj.selectVar = QComboBox()
        plotObj.selectVar.setStyleSheet(Layout.QComboBox())
        plotObj.selectVar.addItems(self.variableList)
        selectVarWidget = QWidget()
        selectVarLabel = QLabel()
        selectVarLabel.setText('Variable:')
        selectVarLayout = QHBoxLayout()
        selectVarLayout.addWidget(selectVarLabel)
        selectVarLayout.addWidget(plotObj.selectVar)
        selectVarWidget.setLayout(selectVarLayout)
        plotObj.selectVar.activated.connect(plotObj.selectionChangeVar)

        #Derived Variables
        #self.selectdVar = QComboBox()
        #self.selectdVar.addItems(self.dvarlist)
        #selectdVarWidget = QWidget()
        #selectdVarLabel = QLabel()
        #selectdVarLabel.setText('Derived Var:')
        #selectdVarLayout = QHBoxLayout()
        #selectdVarLayout.addWidget(selectdVarLabel)
        #selectdVarLayout.addWidget(self.selectdVar)
        #selectdVarWidget.setLayout(selectdVarLayout)
        ##selectdVarWidget.setFixedWidth(200)
        #self.selectdVar.activated.connect(plotObj.selectionChangedVar)

        #Grid selection -- add later
        #self.selectGrid = QComboBox()
        #self.selectGrid.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        #for i in range(0,plotObj.dataSet.numGrids):
        #    self.selectGrid.addItem(str(i+1))
        self.selectGridWidget = QWidget()
        #selectGridLabel = QLabel()
        #selectGridLabel.setText('Grid:')
        selectGridLayout = QHBoxLayout()
        #selectGridLayout.addWidget(selectGridLabel)
        #selectGridLayout.addWidget(self.selectGrid)

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
        #self.gboxLayout.addWidget(selectdVarWidget)
        self.gboxLayout.addWidget(self.selectGridWidget)
        #self.selectGrid.activated.connect(plotObj.selectionChangeGrid)
        plotObj.selectLevel.currentIndexChanged.connect(plotObj.selectionChangeLevel)

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
     
        #plotObj.qscrollLayout.addWidget(self.gbox,Qt.AlignTop)

        return self.qscroll
###############################################################################
####                          End CMAQDataset() Object                      ####
###############################################################################
