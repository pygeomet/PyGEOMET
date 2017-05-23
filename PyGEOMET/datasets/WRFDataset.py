import numpy as np
import os
import glob
import netCDF4
import multiprocessing
import time
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import PyGEOMET.utils.wrf_functions as wrf
import PyGEOMET.utils.LayoutFormat as Layout

#This function is used to create the time list using multiprocessing
def read_nctime(infile):
    ncId = netCDF4.Dataset(infile,'r')
    variables = ncId.variables
    ntimes = ncId.dimensions['Time'].size
    #timeString = []
    #Read the variable
    if ntimes > 1:
        timeString = []
        for ctime in range(ntimes):
            tm = variables['Times'][ctime]
            #Create time string
            year   = tm[0]+tm[1]+tm[2]+tm[3]
            month  = tm[5]+tm[6]
            day    = tm[8]+tm[9]
            hour   = tm[11]+tm[12]
            minute = tm[14]+tm[15]
            second = tm[17]+tm[18]
            #if ideal, python2 does not accept year=1, so we force it
            if int(year) < 1900:
                year = '1901'
            timeObj  = datetime(
                                int(year),
                                int(month),
                                int(day),
                                int(hour),
                                int(minute),
                                int(second)
                                )
            timeString.append(timeObj.strftime("%d %b %Y, %H:%M:%S UTC"))
            #timeString = timeObj.strftime("%d %b %Y, %H:%M:%S UTC")
    else:
        tm = variables['Times'][0]
        #Create time string
        year   = tm[0]+tm[1]+tm[2]+tm[3]
        month  = tm[5]+tm[6]
        day    = tm[8]+tm[9]
        hour   = tm[11]+tm[12]
        minute = tm[14]+tm[15]
        second = tm[17]+tm[18]
        #if ideal, python2 does not accept year=1, so we force it
        if int(year) < 1900:
            year = '1901'
        timeObj  = datetime(
                            int(year),
                            int(month),
                            int(day),
                            int(hour),
                            int(minute),
                            int(second)
                            )
        #timeString.append(timeObj.strftime("%d %b %Y, %H:%M:%S UTC"))
        timeString = timeObj.strftime("%d %b %Y, %H:%M:%S UTC")

    return timeString

class WrfDataset:

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

        self.runType = None

        self.currentTimeIndex = 0

        self.currentFileIndex = 0

        self.units='None'

        self.dsetname = "WRF"

        self.resolution = "l"
 
        #Define plot type available for the dataset within the GUI
        self.ptypes = ['Horizontal Slice', 'Vertical Slice', 'SkewT/Hodograph',
                       'Vertical Profile', 'Time Series', 'Difference Plot']
                       #, 'Hovmoller Diagram']


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

            #Determine 1/2 the number of processors
            procs = int(multiprocessing.cpu_count()/2)
            #Error check if only 1 processor
            if procs == 0:
                procs = 1
            print("num of procs", procs)

            while(True):
                if(i <= 9):
                    files = sorted(glob.glob(self.dsPrefix+'_d0'+str(i)+'*'))
                else:
                    files = sorted(glob.glob(self.dsPrefix+'_d'+str(i)+'*'))

                numFiles = len(files)
                times = []
                append_time = times.extend

                if numFiles != 0 :
                    self.fileList.append(files)
                    self.setGrid(self.numGrids+1)
                    #Multi-proc read
                    pool = multiprocessing.Pool(processes=procs)
                    r = pool.map(read_nctime,files)
                    pool.close()
                    #Flatten list if necessary
                    if any(isinstance(el, list) for el in r):
                        r = [ent for sublist in r for ent in sublist]
                    self.timeList.append(r)
                    self.numGrids += 1
                    self.setTimeIndex(0,update=False)
                totalFiles = totalFiles-numFiles
                i = i+1
                if totalFiles == 0:
                    break;
                #self.setTimeIndex(0,update=False)
            self.setTimeIndex(0,update=False)
            self.setGrid(1)
            self.threeDVars = []
            self.twoDVars = []
            self.variableReadFxn = { }         
            for varname in self.variableList :
                self.variableReadFxn[varname] = self.readNCVariable
                var = self.readNCVariable(varname,varonly=True)
                if len(var.shape) > 2:
                    self.threeDVars.append(varname)
                else:
                    self.twoDVars.append(varname)
            self.setGridDefinition()

            self.dvarlist = wrf.getDvarList(self.variableList,self.runType)

    def readNCVariable(self,vname,barbs=None, vectors=None, contour2=None,varonly=False):
        variable = self.dsets[vname][self.currentTimeIndex]
        if varonly == False:
            if( hasattr(self.dsets[vname],'description')):
                self.description = self.dsets[vname].description
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

        if barbs == True or vectors == True:
            if len(self.dsets[vname].shape) == 4:
                tmp = self.dsets['U'][self.currentTimeIndex]
                self.u10 = wrf.unstaggerX(tmp)
                tmp = self.dsets['V'][self.currentTimeIndex]
                self.v10 = wrf.unstaggerY(tmp)
            else:
                self.u10 = self.dsets['U10'][self.currentTimeIndex]
                self.v10 = self.dsets['V10'][self.currentTimeIndex]        

        return variable

    def getVertVars(self):
        u = wrf.unstaggerX(np.squeeze(self.readNCVariable('U')))
        v = wrf.unstaggerY(np.squeeze(self.readNCVariable('V')))
        w = wrf.unstaggerZ(np.squeeze(self.readNCVariable('W')))
        ph = np.squeeze(self.readNCVariable('PH'))
        phb = np.squeeze(self.readNCVariable('PHB'))
        p = np.squeeze(self.readNCVariable('P'))
        pb = np.squeeze(self.readNCVariable('PB'))
        press = (p + pb)/100.
        height = wrf.unstaggerZ((ph + phb)/9.81)/1000.
        
        return u, v, w, press, height

    def getVariable(self,vname):
        return self.variableReadFxn[vname](vname)

    def addReadFxn(self,vname,fxn):
        self.variableReadFxn[vname] = fxn
        self.variableList.append(vname)

    def getAttribute(attname):
        return self.ncId.__dict__[attname]

    def setGridDefinition(self):
        self.gridRatio = []
        self.pID = []
        self.nx  = []
        self.ny  = []
        self.nz  = []
        self.dx  = []
        self.dy  = []
        self.istart = []
        self.jstart = []
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
            if i == 0:
                self.setTimeIndex(0,update=False)
                self.setGrid(i+1)
            else:
                self.setGrid(i+1,update=False)
                self.setTimeIndex(indx)

            self.setGrid(i+1)
            self.setTimeIndex(0)
            self.gridRatio.append(self.ncId.__dict__['PARENT_GRID_RATIO'])
            self.pID.append(self.ncId.__dict__['PARENT_ID'])
            self.nx.append(self.ncId.__dict__['WEST-EAST_GRID_DIMENSION'])
            self.ny.append(self.ncId.__dict__['SOUTH-NORTH_GRID_DIMENSION'])
            self.nz.append(self.ncId.__dict__['BOTTOM-TOP_GRID_DIMENSION'])
            self.dx.append(self.ncId.__dict__['DX'])
            self.dy.append(self.ncId.__dict__['DY'])
            self.istart.append(self.ncId.__dict__['I_PARENT_START'])
            self.jstart.append(self.ncId.__dict__['J_PARENT_START'])
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
        #if ideal, python2 does not accept year=1, so we force it
        if int(year) < 1900:
            year = '1901'
        self.timeObj  = datetime(
                                     int(year),
                                     int(month),
                                     int(day),
                                     int(hour),
                                     int(minute),
                                     int(second)
                                     )
        if self.runType == 'REAL':
            self.timeString = self.timeObj.strftime("%d %b %Y, %H:%M:%S UTC")
        else:
            tmp=timedelta(days=int(day),hours=int(hour),
                                   minutes=int(minute),seconds=int(second))
            #stupid datetime timedelta stuff
            tmp1 = tmp.days*24-24
            tmp2 = (tmp.seconds - np.mod(tmp.seconds,3600))/3600.
            tmp3 = (tmp.seconds-(tmp2*3600)-np.mod(tmp.seconds-(tmp2*3600),60))/60.
            tmp4 = tmp.seconds-(tmp2*3600)-(tmp3*60)
            self.timeString = ('Forecast Hour ' +
                              "{0:0>3}".format(int(tmp1+tmp2)) + 
                              ':' + "{0:0>2}".format(int(tmp3)) +
                              ':' + "{0:0>2}".format(int(tmp4)))
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
        self.test = axs
        
    def setGridCorners(self):

        # Save current state.

        cg = self.currentGrid
        indx = self.currentTimeIndex

        if (self.ncId.variables['XLONG'][:].min() ==
            self.ncId.variables['XLONG'][:].max() and
            self.ncId.variables['XLAT'][:].min() ==
            self.ncId.variables['XLAT'][:].max()):
            self.runType = 'IDEAL'
        else:
            self.runType = 'REAL'


        for i in range(0,self.numGrids):

            if i == 0:
                self.setTimeIndex(0,update=False)
                self.setGrid(i+1)
            else:
                self.setGrid(i+1,update=False)
                self.setTimeIndex(indx)
            
            if self.runType == 'REAL':
                self.glons[i] = self.readNCVariable('XLONG')
                dims = self.glons[i].shape
                self.ll_lon[i] = self.glons[i][0,0]
                self.ur_lon[i] = self.glons[i][dims[0]-1,dims[1]-1]

                self.glats[i] = self.readNCVariable('XLAT')
                dims = self.glats[i].shape
                self.ll_lat[i] = self.glats[i][0,0]
                self.ur_lat[i] = self.glats[i][dims[0]-1,dims[1]-1]

                self.wd[i] = (self.nx[i]-2)*self.dx[i]
                self.ht[i] = (self.ny[i]-2)*self.dy[i]
                self.setProjection(i+1)

            elif self.runType == 'IDEAL':
                xx = np.arange(-1*((self.nx[i]-1)/2)*self.dx[i],((self.nx[i]-1)/2)*self.dx[i],self.dx[i])
                yy = np.arange(-1*((self.ny[i]-1)/2)*self.dy[i],((self.ny[i]-1)/2)*self.dy[i],self.dy[i])
                self.glons[i],self.glats[i] = np.meshgrid(xx/1000.,yy/1000.)

            else:
                print('This does seem to be a REAL or IDEAL WRF Simulation.' +
                      ' PyGEOMET has not been configured to work with this file' +
                       '. Exiting Now')
                exit()

        self.setGrid(cg,update=False)
        self.setTimeIndex(indx)

#####################  End of function setGridCorners() #######################

    def getDsetControlBar(self, plotObj):
        #self.DCBar = QWidget()
        #plotObj.MBar = QGroupBox(self.DCBar)
        self.tabbing = QWidget()
        plotObj.tabbingLayout = QVBoxLayout() 
        self.tabbing.setLayout(plotObj.tabbingLayout)
        #self.qscroll = QScrollArea(plotObj.MBar)
        #self.qscroll = QScrollArea(self.DCBar)
        self.qscroll = QScrollArea()
        #self.qscroll.setMinimumSize(310,620)
        #self.qscroll.setMaximumSize(310,620)        
        #print("here", self.qscroll.geometry())


        qscrollContents = QWidget()
        plotObj.qscrollLayout = QVBoxLayout(qscrollContents)

        self.qscroll.setWidget(qscrollContents)
        self.qscroll.setWidgetResizable(True)

        #qfigWidget = QWidget(qscrollContents)
        #plotObj.qscrollLayout.addWidget(qfigWidget)
        #qscrollContents.setLayout(plotObj.qscrollLayout)
        
        #boxTitleString = 'Plot # ' +str(plotObj.plotCount)
        #self.gbox = QGroupBox(boxTitleString,self.DCBar)
        self.gbox = QGroupBox()

        #font = QFont();
        #font.setBold(True);
        #self.gbox.setFont(font);
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
        self.selectPlotType.activated.connect(plotObj.selectionChangePlot)
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
        self.selectdVar = QComboBox()
        self.selectdVar.setStyleSheet(Layout.QComboBox())
        self.selectdVar.addItems(self.dvarlist)
        selectdVarWidget = QWidget()
        selectdVarLabel = QLabel()
        selectdVarLabel.setText('Derived Var:')
        selectdVarLayout = QHBoxLayout()
        selectdVarLayout.addWidget(selectdVarLabel)
        selectdVarLayout.addWidget(self.selectdVar)
        selectdVarWidget.setLayout(selectdVarLayout)
        #selectdVarWidget.setFixedWidth(200)
        self.selectdVar.activated.connect(plotObj.selectionChangedVar)

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
        self.gboxLayout.addWidget(selectdVarWidget)
        self.gboxLayout.addWidget(self.selectGridWidget)
        self.selectGrid.currentIndexChanged.connect(plotObj.selectionChangeGrid)
        plotObj.selectLevel.activated.connect(plotObj.selectionChangeLevel)

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
        #self.gbox.autosize(True)       
     
        #print("QScroll",self.qscroll.geometry())
        #print("GBox",self.gbox.geometry()) 
        #self.gbox.setGeometry(0,0,271,240)

        plotObj.tabbingLayout.addWidget(self.gbox)
        #self.optionTabs = QTabWidget()
        #self.optionTabs.setMaximumHeight(plotObj.appobj.screeny*.8*.4)
        plotObj.qscrollLayout.addWidget(self.tabbing,Qt.AlignTop)
        #plotObj.qscrollLayout.addWidget(self.tabbing,Qt.AlignCenter)
        #tab1 = QWidget()
        #tab2 = QWidget()
        #self.test2.addTab(tab1,"testing1")
        #self.test2.addTab(tab2,"testing2")
        #self.testLayout.addWidget(self.test2)        
        
        return self.qscroll

###############################################################################
####                          End WRFDataset() Object                      ####
###############################################################################

