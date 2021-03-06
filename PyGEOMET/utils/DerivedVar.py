import numpy as np
import PyGEOMET.utils.wrf_functions as wrf
import PyGEOMET.utils.wrf_cython as wrf_cython
try:
    import PyGEOMET.utils.crtm_python as CRTM
except ImportError:
    pass

class WRFDerivedVar:

    def __init__(self, dset = None, var = None, ptype = None, 
                       sensor = None, channel = None, 
                       path = None, req_var = None):
        
        #If input variable is not defined return
        if var == None:
            print("Input varible undefined")
            self.var = None
            self.var2 = None
            self.varTitle = None
        else:
            self.dataSet = dset
            self.ptype = ptype            

            #Set sensor and channel for CRTM
            self.sensor = sensor
            self.channel = channel
            self.main_path = path

            #Set common constants
            self.g = 9.81
            self.rd = 287.04
            self.cp = 1004.
            self.Ep = 18.016/29.87
 
            #Calculate common variables
            ph = self.dataSet.readNCVariable('PH')
            phb = self.dataSet.readNCVariable('PHB')
            p = self.dataSet.readNCVariable('P')
            pb = self.dataSet.readNCVariable('PB')
            t = self.dataSet.readNCVariable('T')
            self.u10 = self.dataSet.readNCVariable('U10')
            self.v10 = self.dataSet.readNCVariable('V10')
            self.qvapor = self.dataSet.readNCVariable('QVAPOR')

            #Make sure all of qvapor is positive
            self.qvapor[np.where(self.qvapor < 0.)] = 0.0

            #Calculate full fields
            self.press = p + pb
            self.height = (ph + phb)/9.81
            self.theta = t + 300.
            self.temp = self.theta * (self.press/100000.)**(self.rd/self.cp)
            self.rho = self.press / (self.rd * (self.temp * (1. + 0.61*self.qvapor)))

            #Determine which derived variable to calculate based on input
            if var == "300mb_winds":
                self.winds_300mb()
            elif var == "500mb_hgt":
                self.hgt_500mb() 
            elif var == "500mb_temp":
                self.temp_500mb()
            elif var == "500mb_vort":
                self.vort_500mb()
            elif var == "700mb_rh":
                self.rh_700mb()
            elif var == "850mb_temp":
                self.temp_850mb()
            elif var == "850mb_vort":
                self.vort_850mb()
            elif var == "acc_pcp":
                self.acc_pcp()
            elif var == "clwp":
                self.clwp()
            elif var == "iwp":
                self.iwp()
            elif var == "twp":
                self.twp()
            elif var == "divergence":
                self.divergence()
            elif var == "froude no":
                self.froude()
            elif var == "lcl_hgt":
                self.lcl_hgt()
            elif var == "lcl_hgt2":
                self.lcl_hgt2()
            elif var == "mslp":
                self.mslp()
            elif var == "rainrate":
                self.precip_rate() 
            elif var == "pwat":
                self.pwat()
            elif var == "pot_vort":
                self.PotentialVorticity()
            elif var == "refl":
                self.refl()
            elif var == "richardson":
                self.richardson()
            elif var == "bulk_rich":
                self.bulk_rich()
            elif var == "shear0_1km":
                self.shear0_1km()
            elif var == "shear0_3km":
                self.shear0_3km()
            elif var == "shear0_6km":
                self.shear0_6km()
            elif var == "t2m":
                self.t2m()
            elif var == "td2m":
                self.td2m()
            elif var == "temp-3d":
                self.temp_3d()
            elif var == "wind10m":
                self.wind10m()
            elif var == "wind-3d":
                self.wind_3d()
            elif var == "CAPE_SB":
                self.cape_SB()
            elif var == "CAPE_ML":
                self.cape_ML()
            elif var == "CAPE_MU":
                self.cape_MU()
            elif var == "CIN_SB":
                self.cin_SB()
            elif var == "CIN_ML":
                self.cin_ML()
            elif var == "CIN_MU":
                self.cin_MU()
            elif var == "theta":
                self.ptheta()
            elif var == "theta-e-3d":
                self.theta_e()
            elif var == "cloud water":
                self.cloud_water()
            elif var == "ice water":
                self.ice_water()
            elif var == "total water":
                self.total_water()
            elif var == "cloud albedo":
                self.cloud_albedo()
            elif var == "temp_ML":
                self.mean_layer_temp()
            elif var == "BrightTemp/Radiance":
                self.crtm_wrapper(req_var)
            elif var == 'RH':
                self.relative_humidity()
            elif var == 'Total Pressure':
                self.total_pressure()
            elif var == 'Sat_QVapor':
                self.sat_qvapor()
            else:
                print("Cannot find input variable")
                self.var = None
                self.var2 = None
                self.varTitle = None

    def total_pressure(self):
        self.var = self.press/100.
        self.var2 = self.var
        self.varTitle = "Total Pressure [hPa]\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Total Pressure [hPa]"

    def winds_300mb(self):
        
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        height = wrf.unstaggerZ(self.height)
        ref_val = 30000.
        #Switched to Cython
        #self.u10 = wrf.loglinear_interpolate(u_corr, self.press, ref_val)
        #self.v10 = wrf.loglinear_interpolate(v_corr, self.press, ref_val)
        self.u10 = np.array(wrf_cython.loglinear_interpolate(u_corr, self.press, ref_val))
        self.v10 = np.array(wrf_cython.loglinear_interpolate(v_corr, self.press, ref_val))
        var1 = wrf.get_bulk_wind(self.u10,self.v10)
        self.var = wrf.convertWind_MStoKT(var1)
        #self.var2 = wrf.hypsometric(height, self.press, ref_val, self.temp)
        self.var2 = np.array(wrf_cython.hypsometric(height, self.press, ref_val, self.temp))
        self.varTitle = "300-mb Wind\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "300-mb Wind"

    def hgt_500mb(self):

        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        height = wrf.unstaggerZ(self.height)
        ref_val = 50000.
        #Switched to Cython
        #self.u10 = wrf.loglinear_interpolate(u_corr, self.press, ref_val)
        #self.v10 = wrf.loglinear_interpolate(v_corr, self.press, ref_val)
        self.u10 = np.array(wrf_cython.loglinear_interpolate(u_corr, self.press, ref_val))
        self.v10 = np.array(wrf_cython.loglinear_interpolate(v_corr, self.press, ref_val))
        #self.var = wrf.hypsometric(height,self.press, ref_val, self.temp)
        self.var = np.array(wrf_cython.hypsometric(height, self.press, ref_val, self.temp))
        self.var2 = self.var
        self.varTitle = "500-mb Geopotential Height (m)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "500-mb Geopotential Height (m)"

    def temp_500mb(self):
       
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        height = wrf.unstaggerZ(self.height)
        ref_val = 50000.
        #Switched to Cython
        #var1 = wrf.loglinear_interpolate(self.temp, self.press, ref_val)
        #self.u10 = wrf.loglinear_interpolate(u_corr, self.press, ref_val)
        #self.v10 = wrf.loglinear_interpolate(v_corr, self.press, ref_val)
        var1 = np.array(wrf_cython.loglinear_interpolate(self.temp, self.press, ref_val))
        self.u10 = np.array(wrf_cython.loglinear_interpolate(u_corr, self.press, ref_val))
        self.v10 = np.array(wrf_cython.loglinear_interpolate(v_corr, self.press, ref_val))
        self.var = wrf.convertT_KtoC(var1)
        #self.var2 = wrf.hypsometric(height, self.press, ref_val, self.temp)
        self.var2 = np.array(wrf_cython.hypsometric(height, self.press, ref_val, self.temp))
        self.varTitle = "500-mb Temperature ($^{\circ}$C)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "500-mb Temperature ($^{\circ}$C)"

    def vort_500mb(self):

        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        height = wrf.unstaggerZ(self.height)
        ref_val = 50000.     
        #Switched to Cython
        #self.u10 = wrf.loglinear_interpolate(u_corr, self.press, ref_val)
        #self.v10 = wrf.loglinear_interpolate(v_corr, self.press, ref_val)
        self.u10 = np.array(wrf_cython.loglinear_interpolate(u_corr, self.press, ref_val))
        self.v10 = np.array(wrf_cython.loglinear_interpolate(v_corr, self.press, ref_val))
        self.var = wrf.rel_vort(self.u10, self.v10,
                   self.dataSet.dx[self.dataSet.currentGrid-1],
                   self.dataSet.dy[self.dataSet.currentGrid-1])
        #self.var2 = wrf.hypsometric(height, self.press, ref_val, self.temp)
        self.var2 = np.array(wrf_cython.hypsometric(height, self.press, ref_val, self.temp))
        self.varTitle = "500-mb Relative Vorticity ($10^{-5}$ $s^{-1}$)\n" +\
                         self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "500-mb Relative Vorticity ($10^{-5}$ $s^{-1}$)"

    def rh_700mb(self):

        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        ref_val = 70000.
        #Switched to Cython 
        #self.u10 = wrf.loglinear_interpolate(u_corr, self.press, ref_val)
        #self.v10 = wrf.loglinear_interpolate(v_corr, self.press, ref_val)
        self.u10 = np.array(wrf_cython.loglinear_interpolate(u_corr, self.press, ref_val))
        self.v10 = np.array(wrf_cython.loglinear_interpolate(v_corr, self.press, ref_val))
        var1 = wrf.get_rh(self.temp, self.qvapor, self.press)
#        self.var = wrf.loglinear_interpolate(var1, self.press, ref_val)
        self.var = np.array(wrf_cython.loglinear_interpolate(var1, self.press, ref_val))
        self.var2 = self.var
        self.varTitle = "700-mb Relative Humidity (%)\n" + self.dataSet.getTime() 
        #Set short variable title for time series
        self.sTitle = "700-mb Relative Humidity (%)"

    def temp_850mb(self):

        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        ref_val = 85000.
        #Switched to Cython
        #self.u10 = wrf.loglinear_interpolate(u_corr, self.press, ref_val)
        #self.v10 = wrf.loglinear_interpolate(v_corr, self.press, ref_val)
        self.u10 = np.array(wrf_cython.loglinear_interpolate(u_corr, self.press, ref_val))
        self.v10 = np.array(wrf_cython.loglinear_interpolate(v_corr, self.press, ref_val))
        #var1 = wrf.loglinear_interpolate(self.temp, self.press, ref_val)
        var1 = np.array(wrf_cython.loglinear_interpolate(self.temp, self.press, ref_val))
        self.var = wrf.convertT_KtoC(var1)
        self.var2 = wrf.get_mslp(self.height, self.press, self.temp, self.qvapor)
        self.varTitle = "850-mb Temperature\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "850-mb Temperature"

    def vort_850mb(self):
        
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        ref_val = 85000.
        #Switched to Cython 
        #self.u10 = wrf.loglinear_interpolate(u_corr, self.press, ref_val)
        #self.v10 = wrf.loglinear_interpolate(v_corr, self.press, ref_val)
        self.u10 = np.array(wrf_cython.loglinear_interpolate(u_corr, self.press, ref_val))
        self.v10 = np.array(wrf_cython.loglinear_interpolate(v_corr, self.press, ref_val))
        self.var = wrf.rel_vort(self.u10,self.v10,
                   self.dataSet.dx[self.dataSet.currentGrid-1],
                   self.dataSet.dy[self.dataSet.currentGrid-1])
        self.var2 = wrf.get_mslp(self.height, self.press, self.temp, self.qvapor)
        self.varTitle = "850-mb Relative Voriticity ($10^{-5}$ $s^{-1}$)\n" +\
                        self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "850-mb Relative Voriticity ($10^{-5}$ $s^{-1}$)"

    def acc_pcp(self):

        rainc = self.dataSet.readNCVariable('RAINC')
        rainnc = self.dataSet.readNCVariable('RAINNC')
        rainsh = self.dataSet.readNCVariable('RAINSH')
        self.var = rainc + rainnc + rainsh
        #self.var = wrf.convertP_MMtoIN(var1)
        self.var2 = self.var
        self.varTitle = "Total Accumulated Precipitation (mm)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Total Accumulated Precipitation (mm)"

    def clwp(self):

        qcloud = self.dataSet.readNCVariable('QCLOUD')
        qrain = self.dataSet.readNCVariable('QRAIN')
        qsnow = self.dataSet.readNCVariable('QSNOW')
        qice = self.dataSet.readNCVariable('QICE')
        qgraup = self.dataSet.readNCVariable('QGRAUP')
        self.var = wrf.get_lwp(self.temp, self.height, self.rho, qcloud, qice, qgraup,
                               qsnow, qrain, 'liquid')
        self.var2 = self.var
        self.varTitle = "Cloud Water Path (g m$^{-2}$)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Cloud Water Path (g m$^{-2}$)"

    def iwp(self):

        qcloud = self.dataSet.readNCVariable('QCLOUD')
        qrain = self.dataSet.readNCVariable('QRAIN')
        qsnow = self.dataSet.readNCVariable('QSNOW')
        qice = self.dataSet.readNCVariable('QICE')
        qgraup = self.dataSet.readNCVariable('QGRAUP')
        self.var = wrf.get_lwp(self.temp, self.height, self.rho, qcloud, qice, qgraup,
                               qsnow, qrain, 'ice')
        self.var2 = self.var
        self.varTitle = "Ice Water Path (g m$^{-2}$)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Ice Water Path (g m$^{-2}$)"

    def twp(self):

        qcloud = self.dataSet.readNCVariable('QCLOUD')
        qrain = self.dataSet.readNCVariable('QRAIN')
        qsnow = self.dataSet.readNCVariable('QSNOW')
        qice = self.dataSet.readNCVariable('QICE')
        qgraup = self.dataSet.readNCVariable('QGRAUP')
        self.var = wrf.get_lwp(self.temp, self.height, self.rho, qcloud, qice, qgraup,
                               qsnow, qrain, 'total')
        self.var2 = self.var
        self.varTitle = "Total Water Path (g m$^{-2}$)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Total Water Path (g m$^{-2}$)"

    def divergence(self):
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        #For the WRF-ARW Arakawa C grid we don't need to use the unstagger winds
        # to calculate the gradient at the mass points
        #du = np.gradient(u_corr)[2]
        #dv = np.gradient(v_corr)[1]
        #Get dimensions
        dimsu = u.shape
        dimsv = v.shape
        du = u[:,:,1:dimsu[2]] - u[:,:,0:dimsu[2]-1]
        dv = v[:,1:dimsv[1],:] - v[:,0:dimsv[1]-1,:]
        dx = self.dataSet.dx[self.dataSet.currentGrid-1]
        dy = self.dataSet.dy[self.dataSet.currentGrid-1]
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = (du/dx + dv/dy)*pow(10,5)
        self.var2 = self.var
        self.varTitle = "Divergence (10$^{-5}$ s$^{-1}$)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Divergence (10$^{-5}$ s$^{-1}$)"

    def froude(self):
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        self.u10 = u_corr
        self.v10 = v_corr
        wind = (u_corr*u_corr + v_corr*v_corr)**(0.5)
        height = wrf.unstaggerZ(self.height)
        
        #compute the vertical gradient of theta
        dtheta = np.gradient(self.theta)[0]
        dz = np.gradient(height)[0]
 
        #calculate the Brunt-Vaisala Frequency
        dtheta_dz = dtheta/dz
        dims = dtheta_dz.shape
        N = np.zeros((dims[0],dims[1],dims[2]))
        #Account for statically unstable and neutral conditions
        N[dtheta_dz<=0] = 0.00001
        N[dtheta_dz>0] = ((self.g/self.theta[dtheta_dz>0])*(dtheta_dz[dtheta_dz>0]))**(0.5)
        self.var = wind/(height*N)
        self.var2 = self.var
        self.varTitle = "Froude Number\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Froude Number"

    def lcl_hgt(self):
 
        t2 = self.dataSet.readNCVariable('T2')
        q2 = self.dataSet.readNCVariable('Q2')
        psfc = self.dataSet.readNCVariable('PSFC')
        td2 = wrf.get_td(t2,q2,psfc)
        rh = wrf.get_rh(t2,q2,psfc)
        self.var = wrf.get_lcl(t2,rh)
        self.var2 = self.var
        self.varTitle = "LCL Height (m)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "LCL Height (m)"

    def mslp(self):
        
        self.var = wrf.get_mslp(self.height,self.press,self.temp,self.qvapor)
        self.var2 = self.var
        self.varTitle = "Mean Sea Level Pressure (hPa)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Mean Sea Level Pressure (hPa)"

    def precip_rate(self):
        if self.dataSet.ntimes == 1:
            prevind = self.dataSet.currentFileIndex-1
            currind = self.dataSet.currentFileIndex
        elif len(self.dataSet.fileList) > 1:
            prevind = self.dataSet.currentFileIndex*self.dataSet.ntimes + self.dataSet.currentTimeIndex-1
            currind = self.dataSet.currentFileIndex*self.dataSet.ntimes + self.dataSet.currentTimeIndex
            print(prevind,currind)
        else:
            prevind = self.dataSet.currentTimeIndex-1
            currind = self.dataSet.currentTimeIndex
        self.dataSet.setTimeIndex(currind)
        currtime = self.dataSet.readNCVariable('XTIME')
        tmp1 = self.dataSet.readNCVariable('RAINC')
        tmp2 = self.dataSet.readNCVariable('RAINNC')
        tmp3 = self.dataSet.readNCVariable('RAINSH')

        current = tmp1 + tmp2 + tmp3
        
        self.dataSet.setTimeIndex(prevind)
        if prevind == -1:
            prev = 0
            prevtime = currtime
        else:
            tmp = self.dataSet.getTime()
            prevtime = self.dataSet.readNCVariable('XTIME')
            tmp1 = self.dataSet.readNCVariable('RAINC')
            tmp2 = self.dataSet.readNCVariable('RAINNC')
            tmp3 = self.dataSet.readNCVariable('RAINSH')
            prev = tmp1 + tmp2 + tmp3

        self.dataSet.setTimeIndex(currind)
        dt = currtime - prevtime
        if dt == 0:
            self.var = current - prev
        else:
            self.var = (current - prev)/(dt/60.)
        self.var2 = self.var
        self.varTitle = "Rain Rate (mm hr$^{-1}$)\n" + self.dataSet.getTime() 
        #Set short variable title for time series
        self.sTitle = "Rain Rate (mm hr$^{-1}$)"

    def pwat(self):

        self.var = wrf.get_pwat(self.qvapor, self.rho, self.height)
        self.var2 = wrf.get_mslp(self.height, self.press, self.temp, self.qvapor)
        self.varTitle = "Total Precipitable Water (in)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Total Precipitable Water (in)"

    def refl(self):

        if 'REFL_10CM' not in self.dataSet.variableList:
            qrain = self.dataSet.readNCVariable('QRAIN')
            qsnow = self.dataSet.readNCVariable('QSNOW')
            qgraup = self.dataSet.readNCVariable('QGRAUP')
            height = wrf.unstaggerZ(self.height)
            self.var = wrf.get_refl(qrain, qgraup, qsnow, self.rho, self.temp, self.height)
            self.var2 = self.var
        if 'REFL_10CM' in self.dataSet.variableList:
            tmp = self.dataSet.readNCVariable('REFL_10CM')
            if self.dataSet.dx[self.dataSet.currentGrid-1] < 12000.:
                #get the index 1.5-km AGL
                height = wrf.unstaggerZ(self.height)
                diff = (self.height) - (self.height[0,:,:]+1500.)
                zz = min(np.where(diff > 0)[0])
                self.var = np.amax(tmp[0:zz,:,:], axis=0)
            else:
                 self.var = np.amax(tmp,axis=0)
            self.var2 = self.var
        self.varTitle = "0-1km Simulated Radar Reflectivity (dBZ)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "0-1km Simulated Radar Reflectivity (dBZ)"

    def richardson(self):
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        wind = (u_corr*u_corr + v_corr*v_corr)**(0.5)
        height = wrf.unstaggerZ(self.height)

        #compute the vertical gradient of theta
        dtheta = np.gradient(self.theta)[0]
        du = np.gradient(wind)[0]
        dz = np.gradient(height)[0]
       
        #compute the richardson number
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = ((self.g/self.theta)*(dtheta/dz))/(du/dz)**(2)
        self.var2 = self.var
        self.varTitle = "Richardson Number\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Richardson Number"

    def bulk_rich(self):
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        qv = self.dataSet.readNCVariable('QVAPOR')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        tv = self.temp*(1+0.611*qv)
        thetav = self.theta*(1+0.611*qv)
        height = wrf.unstaggerZ(self.height)

        #compute the vertical gradient of theta
        dtheta = np.gradient(thetav)[0]
        dz = np.gradient(height)[0]
        du = np.gradient(u_corr)[0]
        dv = np.gradient(v_corr)[0]
        print(np.min(dtheta),np.max(dtheta),np.min(dz),np.max(dz))
        print(np.mean(du),np.mean(dv))
        print(np.min(du),np.max(du),np.min(dv),np.max(dv))
        #compute bulk richardson number
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = ((self.g/tv)*(dtheta)*dz)/((du)**2+(dv)**2)
        self.var2 = self.var
        self.varTitle = 'Bulk Richardson number\n' + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = 'Bulk Richardson number'

    def shear0_1km(self):

        #get the height array in reference to ground level
        height = wrf.unstaggerZ(self.height) - wrf.unstaggerZ(self.height)[0,:,:]
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        ref_val = 1000. # meters
        self.u10, self.v10, var1 = wrf.get_shear(u, v, height, ref_val)
        self.var = wrf.convertWind_MStoKT(var1)
        self.var2 = self.var
        self.varTitle = "0-1km Shear(knots)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "0-1km Shear(knots)"

    def shear0_3km(self):

        #get the height array in reference to ground level
        height = wrf.unstaggerZ(self.height) - wrf.unstaggerZ(self.height)[0,:,:]
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        ref_val = 3000. # meters
        self.u10, self.v10, var1 = wrf.get_shear(u, v, height, ref_val)
        self.var = wrf.convertWind_MStoKT(var1)
        self.var2 = self.var
        self.varTitle = "0-3km Shear (knots)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "0-3km Shear (knots)"

    def shear0_6km(self):

        #get the height array in reference to ground level
        height = wrf.unstaggerZ(self.height) - wrf.unstaggerZ(self.height)[0,:,:]
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        ref_val = 6000. # meters
        self.u10, self.v10, var1 = wrf.get_shear(u, v, height, ref_val)
        self.var = wrf.convertWind_MStoKT(var1)
        self.var2 = self.var
        self.varTitle = "0-6km Shear (knots)\n" + self.dataSet.getTime()        
        #Set short variable title for time series
        self.sTitle = "0-6km Shear (knots)"

    def t2m(self):
    
        t2 = self.dataSet.readNCVariable('T2')
        self.var = wrf.convertT_KtoF(t2)
        self.var2 = self.var
        self.varTitle = "2-meter Temperature ($^{\circ}$ F)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "2-meter Temperature ($^{\circ}$ F)"

    def td2m(self):

        t2 = self.dataSet.readNCVariable('T2')
        q2 = self.dataSet.readNCVariable('Q2')
        psfc = self.dataSet.readNCVariable('PSFC')
        var1 = wrf.get_td(t2, q2, psfc)
        self.var = wrf.convertT_KtoF(var1)
        self.var2 = self.var
        self.varTitle = "2-meter Dewpoint Temperature ($^{\circ}$ F)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "2-meter Dewpoint Temperature ($^{\circ}$ F)"

    def temp_3d(self):

        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = self.temp
        self.var2 = self.var
        self.varTitle = "Temperature (K)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Temperature (K)"

    def wind10m(self):

        var1 = wrf.get_bulk_wind(self.u10, self.v10)
        #self.var = wrf.convertWind_MStoKT(var1)
        self.var = var1
        self.var2 = self.var
        #self.varTitle = "10-meter Wind (knots)\n"+ self.dataSet.getTime()
        self.varTitle = "10-meter Wind (m s$^{-1}$)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "10-meter Wind (m s$^{-1}$)"

    def wind_3d(self):

        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        var1 = wrf.get_bulk_wind(u_corr,v_corr)
        #self.var = wrf.convertWind_MStoKT(var1)
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = var1
        self.var2 = self.var
        #self.varTitle = "Total Wind (knots)\n"+ self.dataSet.getTime()
        self.varTitle = "Total Wind (m s$^{-1}$)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Total Wind (m s$^{-1}$)"

    def cape_SB(self):

        #Switched to Cython        
        #self.var, self.var2 = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'surface_based')
        height = wrf.unstaggerZ(self.height)
        self.var, self.var2 = wrf_cython.cape_sb(self.temp,self.qvapor,
                              self.press,height)        
        self.varTitle = "Surface-Based CAPE (J kg$^{-1}$)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Surface-Based CAPE (J kg$^{-1}$)"

    def cape_ML(self):

        #Switched to Cython
        #self.var, self.var2 = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'mixed_layer')
        height = wrf.unstaggerZ(self.height)
        self.var, self.var2 = wrf_cython.cape_ml(self.temp,self.qvapor,
                              self.press,height)
        self.varTitle = "Mixed Layer CAPE (J kg$^{-1}$)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Mixed Layer CAPE (J kg$^{-1}$)"

    def cape_MU(self):
        #Switched to Cython
        #self.var, self.var2 = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'most_unstable')
        height = wrf.unstaggerZ(self.height)
        self.var, self.var2 = wrf_cython.cape_mu(self.temp,self.qvapor,
                              self.press,height)
        self.varTitle = "Most Unstable CAPE (J kg$^{-1}$)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Most Unstable CAPE (J kg$^{-1}$)"

    def cin_SB(self):
        #Switched to Cython
        #self.var2, self.var = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'surface_based')
        height = wrf.unstaggerZ(self.height)
        self.var2, self.var = np.array(wrf_cython.cape_sb(self.temp,self.qvapor,
                              self.press,height))        
        self.var *= -1
        self.varTitle = "Surface-Based CIN (J kg$^{-1}$)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Surface-Based CIN (J kg$^{-1}$)"

    def cin_ML(self):
        #Sitched to Cython
        #self.var2, self.var = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'mixed_layer')
        height = wrf.unstaggerZ(self.height)
        self.var2, self.var = np.array(wrf_cython.cape_ml(self.temp,self.qvapor,
                              self.press,height))
        self.var *= -1
        self.varTitle = "Mixed Layer CIN (J kg$^{-1}$)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Mixed Layer CIN (J kg$^{-1}$)"

    def cin_MU(self):
        #Switched to Cython
        #self.var2, self.var = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'most_unstable')
        height = wrf.unstaggerZ(self.height)
        self.var2, self.var = np.array(wrf_cython.cape_mu(self.temp,self.qvapor,
                              self.press,height))
        self.var *= -1
        self.varTitle = "Most Unstable CIN (J kg$^{-1}$)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Most Unstable CIN (J kg$^{-1}$)"

    def ptheta(self):
        self.var = self.theta
        self.var2 = self.theta
        self.varTitle = "$\mathsf{\Theta}$ (K)\n" + self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "$\mathsf{\Theta}$ (K)"

    def theta_e(self):
        
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = wrf.get_thetae(self.temp,self.qvapor,self.press)
        self.var2 = self.var
        self.varTitle = "$\mathsf{\Theta_{e}\, (K)}$\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "$\mathsf{\Theta_{e}\, (K)}$"

    def PotentialVorticity(self):

        #Read in variables
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        f = self.dataSet.readNCVariable('F')    
 
        #Define grid spacing in meters
        dx = self.dataSet.dx[self.dataSet.currentGrid-1] 
        dy = self.dataSet.dy[self.dataSet.currentGrid-1] 

        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = wrf.pot_vort(u,v,f,dx,dy,self.press,self.theta)
        self.var2 = self.var
        self.varTitle = "Potential Voricity (PVU)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Potential Voricity (PVU)"

    def cloud_water(self):
        
        #3D cloud water --> qcloud, qrain
        qcloud = self.dataSet.readNCVariable('QCLOUD')
        qrain = self.dataSet.readNCVariable('QRAIN')
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = qcloud + qrain
        self.var2 = self.var
        self.varTitle = "Total Cloud Water Mixing Ratio (kg kg$^{-1}$)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Total Cloud Water Mixing Ratio (kg kg$^{-1}$)"        

    def ice_water(self):

        #3D ice water --> qice, qsnow, qgraup
        qice = self.dataSet.readNCVariable('QICE')
        qsnow = self.dataSet.readNCVariable('QSNOW')
        qgraup = self.dataSet.readNCVariable('QGRAUP')
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = qice + qsnow + qgraup
        self.var2 = self.var
        self.varTitle = "Total Ice Mixing Ratio (kg kg$^{-1}$)\n"+ self.dataSet.getTime()  
        #Set short variable title for time series
        self.sTitle = "Total Ice Mixing Ratio (kg kg$^{-1}$)"

    def total_water(self):
 
        #3D water --> qrain, qcloud, qice, qsnow, qgraup
        qcloud = self.dataSet.readNCVariable('QCLOUD')
        qrain = self.dataSet.readNCVariable('QRAIN')
        qice = self.dataSet.readNCVariable('QICE')
        qsnow = self.dataSet.readNCVariable('QSNOW')
        qgraup = self.dataSet.readNCVariable('QGRAUP')
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = qice + qsnow + qgraup + qcloud + qrain
        self.var2 = self.var
        self.varTitle = "Total Water Mixing Ratio (kg kg$^{-1}$)\n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Total Water Mixing Ratio (kg kg$^{-1}$)"

    def cloud_albedo(self):
        
        #Read in radiation variables
        #Instantaneous downwelling flux at the bottom
        swdnb = self.dataSet.readNCVariable('SWDNB')
        #Instantaneous downwelling flux at the bottom  - clear sky
        swdnbc = self.dataSet.readNCVariable('SWDNBC')

        #Instantaneous downwelling flux at the top
        #swdnt = self.dataSet.readNCVariable('SWDNT')
        #Instantaneous upwelling flux at the top - clear sky
        #swuptc = self.dataSet.readNCVariable('SWUPTC')
        #Instantaneous upwelling flux at the top
        #swupt = self.dataSet.readNCVariable('SWUPT')

        #Calculate cloud albedo
        self.var = (1 - (swdnb / swdnbc)) * 100.
        self.var2 = self.var
        self.varTitle = "Cloud Albedo \n"+ self.dataSet.getTime() 
        #Set short variable title for time series
        self.sTitle = "Cloud Albedo"

    def mean_layer_temp(self):
        height = wrf.unstaggerZ(self.height)
        pblh = self.dataSet.readNCVariable('PBLH')
        #Calculate mean layer temperature
        #Switch to Cython
        #self.var = wrf.mean_layer(self.temp,height,0,pblh)
        points = pblh.shape
        ref1 = np.zeros((points[1],points[0]),dtype=np.float32)
        self.var = wrf_cython.mean_layer(self.temp,height,ref1,pblh)
        self.var2 = self.var
        self.varTitle = "Mean Layer Temperature (K) \n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Mean Layer Temperature (K)"

    #Function to calculate the relative humidity
    def relative_humidity(self):

        self.var = wrf.get_rh(self.temp, self.qvapor, self.press)
        self.var2 = self.var
        self.varTitle = "Relative Humidity (%) \n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Relative Humidity (%)" 

    #Function to calculate the saturation mixing ratio
    def sat_qvapor(self):

        #Calculate the saturation vapor pressure - Clausius Clapeyron Equation (Bolton, 1980)
        es = 611.2 * np.exp(17.67*(self.temp-273.15) / (self.temp-273.15+243.5))
        #Calculate saturation mixing ratio
        self.var = 0.62197 * (es/(self.press-es))
        self.var2 = self.var
        self.varTitle = "Saturation Mixing Ratio (kg/kg) \n"+ self.dataSet.getTime()
        #Set short variable title for time series
        self.sTitle = "Saturation Mixing Ratio (kg/kg)"

    #Function to determine the brightness temperature or radiance
    def crtm_wrapper(self,request_var):
         
        #Create full path to crtm coefficient data
        full_path = self.main_path + "utils/crtm_coeff/"        

        #Get all necessary variables - put in fortran order
        qcloud = np.array(self.dataSet.readNCVariable('QCLOUD'),order='F')
        qice = np.array(self.dataSet.readNCVariable('QICE'),order='F')
        qrain = np.array(self.dataSet.readNCVariable('QRAIN'),order='F')
        qsnow = np.array(self.dataSet.readNCVariable('QSNOW'),order='F')
        qgraupel = np.array(self.dataSet.readNCVariable('QGRAUP'),order='F')
        lai = np.array(self.dataSet.readNCVariable('LAI'),order='F')
        seaice = np.array(self.dataSet.readNCVariable('SEAICE'),order='F')
        snowh = np.array(self.dataSet.readNCVariable('SNOWH'),order='F')
        coszen = np.array(self.dataSet.readNCVariable('COSZEN'),order='F')
        vegfrac = np.array(self.dataSet.readNCVariable('VEGFRA'),order='F')
        ptop = np.array(self.dataSet.readNCVariable('P_TOP'),order='F')
        tsk = np.array(self.dataSet.readNCVariable('TSK'),order='F')
        ivegtyp = np.array(self.dataSet.readNCVariable('IVGTYP'),order='F')
        xland = np.array(self.dataSet.readNCVariable('XLAND'),order='F')
        landuse = self.dataSet.readNCGlobalAttr('MMINLU')
        mp_physics = self.dataSet.readNCGlobalAttr('MP_PHYSICS')
        #mp_physics = 2
        lon = np.array(self.dataSet.readNCVariable('XLONG'),order='F')
        lat = np.array(self.dataSet.readNCVariable('XLAT'),order='F')
        dims = self.press.shape
        if (mp_physics == 8):
            nrain = np.array(self.dataSet.readNCVariable('QNRAIN'),order='F')
            nice = np.array(self.dataSet.readNCVariable('QNICE'),order='F')
            nsnow = np.zeros((dims[0],dims[1],dims[2]),order='F',dtype='float32')
            ngraupel = np.zeros((dims[0],dims[1],dims[2]),order='F',dtype='float32')
        elif (mp_physics == 10):
            nrain = np.array(self.dataSet.readNCVariable('QNRAIN'),order='F')
            nice = np.array(self.dataSet.readNCVariable('QNICE'),order='F')
            nsnow = np.array(self.dataSet.readNCVariable('QNSNOW'),order='F')
            ngraupel = np.array(self.dataSet.readNCVariable('QNGRAUPEL'),order='F')
        else:
            nrain = np.zeros((dims[0],dims[1],dims[2]),order='F',dtype='float32')
            nice = np.zeros((dims[0],dims[1],dims[2]),order='F',dtype='float32')
            nsnow = np.zeros((dims[0],dims[1],dims[2]),order='F',dtype='float32')
            ngraupel = np.zeros((dims[0],dims[1],dims[2]),order='F',dtype='float32')
        #Current microphysics implemented does not include these variables 
        ncloud = np.zeros((dims[0],dims[1],dims[2]),order='F',dtype='float32')
        nhail = np.zeros((dims[0],dims[1],dims[2]),order='F',dtype='float32')
        #Use global variables
        #self.press, self.theta, self.qvapor, self.u10, self.v10, 
        #self.sensor, self.channel
        #Create dummy array for qhail
        qhail = np.zeros((dims[0],dims[1],dims[2]),order='F',dtype='float32')
        crtm_out = np.zeros((dims[1],dims[2]),order='F',dtype='float32')
        CRTM.crtm(self.press/100.,self.theta,self.qvapor,qcloud,qice,
                  qrain,qsnow,qgraupel,qhail,ncloud,nice,nrain,nsnow,ngraupel,nhail,
                  lai,self.u10,self.v10,seaice,snowh,coszen,vegfrac,ptop,tsk,
                  ivegtyp,xland,landuse,mp_physics,lat,lon,self.sensor,
                  int(self.channel),full_path,request_var,crtm_out,dims[2],dims[1],dims[0])

        #Define variable for output and plot title
        self.var = crtm_out
        self.var2 = self.var
        if (request_var == "Radiance"):
            self.varTitle = "Simulated Radiance [mW $m^{-2}$ $sr^{-1}$ cm] \n"
            #Set short variable title for time series
            self.sTitle = "Simulated Radiance [mW $m^{-2}$ $sr^{-1}$ cm]"
        elif (request_var == "Brightness Temperature"):
            self.varTitle = "Simulated Brightness Temperature [K]\n"
            #Set short variable title for time series
            self.sTitle = "Simulated Brightness Temperature [K]"
        elif (request_var == "Up Radiance"):
            self.varTitle = "Simulated Atm Upward Radiance [mW $m^{-2}$ $sr^{-1}$ cm]\n"
            #Set short variable title for time series
            self.sTitle = "Simulated Atm Upward Radiance [mW $m^{-2}$ $sr^{-1}$ cm]"
        elif (request_var == "Down Radiance"):
            self.varTitle = "Simulated Atm Downward Radiance [mW $m^{-2}$ $sr^{-1}$ cm]\n"
            #Set short variable title for time series
            self.sTitle = "Simulated Atm Downward Radiance [mW $m^{-2}$ $sr^{-1}$ cm]"
        else:
            self.varTitle = "Simulated Downward Solar Radiance [mW $m^{-2}$ $sr^{-1}$ cm]\n"        
            #Set short variable title for time series
            self.sTitle = "Simulated Downward Solar Radiance [mW $m^{-2}$ $sr^{-1}$ cm]"
        self.varTitle = self.varTitle + "(Sensor : "+self.sensor+"    Channel : "+self.channel+")\n"
        self.varTitle = self.varTitle + self.dataSet.getTime()        
          

