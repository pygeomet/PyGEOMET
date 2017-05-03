import numpy as np
import PyGEOMET.utils.wrf_functions as wrf
import PyGEOMET.utils.wrf_cython as wrf_cython
import datetime

class WRFDerivedVar:

    def __init__(self, dset = None, var = None, ptype = None):
        
        #If input variable is not defined return
        if var == None:
            print("Input varible undefined")
            self.var = None
            self.var2 = None
            self.varTitle = None
        else:
            self.dataSet = dset
            self.ptype = ptype            

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
            else:
                print("Cannot find input variable")
                self.var = None
                self.var2 = None
                self.varTitle = None

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

    def acc_pcp(self):

        rainc = self.dataSet.readNCVariable('RAINC')
        rainnc = self.dataSet.readNCVariable('RAINNC')
        rainsh = self.dataSet.readNCVariable('RAINSH')
        var1 = rainc + rainnc + rainsh
        self.var = wrf.convertP_MMtoIN(var1)
        self.var2 = var1
        self.varTitle = "Total Accumulated Precipitation (in)\n" + self.dataSet.getTime()

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

    def divergence(self):
        u = self.dataSet.readNCVariable('U')
        v = self.dataSet.readNCVariable('V')
        u_corr = wrf.unstaggerX(u)
        v_corr = wrf.unstaggerY(v)
        du = np.gradient(u_corr)[2]
        dv = np.gradient(v_corr)[1]
        dx = self.dataSet.dx[self.dataSet.currentGrid-1]
        dy = self.dataSet.dy[self.dataSet.currentGrid-1]
        self.u10 = u_corr
        self.v10 = v_corr
        self.var = (du/dx + dv/dy)*pow(10,5)
        self.var2 = self.var
        self.varTitle = "Divergence (10$^{-5}$ s$^{-1}$)\n" + self.dataSet.getTime()

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

    def lcl_hgt(self):
 
        t2 = self.dataSet.readNCVariable('T2')
        q2 = self.dataSet.readNCVariable('Q2')
        psfc = self.dataSet.readNCVariable('PSFC')
        td2 = wrf.get_td(t2,q2,psfc)
        rh = wrf.get_rh(t2,q2,psfc)
        self.var = wrf.get_lcl(t2,rh)
        self.var2 = self.var
        self.varTitle = "LCL Height (m)\n" + self.dataSet.getTime()

    def mslp(self):
        
        self.var = wrf.get_mslp(self.height,self.press,self.temp,self.qvapor)
        self.var2 = self.var
        self.varTitle = "Mean Sea Level Pressure (hPa)\n" + self.dataSet.getTime()

    def precip_rate(self):
        if self.dataSet.ntimes == 1:
            prevind = self.dataSet.currentFileIndex-1
            currind = self.dataSet.currentFileIndex

        self.dataSet.setTimeIndex(currind)
        tmp = self.dataSet.getTime()
        currtime = datetime.datetime.strptime(tmp, "%d %b %Y, %H:%M:%S UTC")
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
            prevtime = datetime.datetime.strptime(tmp,"%d %b %Y, %H:%M:%S UTC")
            tmp1 = self.dataSet.readNCVariable('RAINC')
            tmp2 = self.dataSet.readNCVariable('RAINNC')
            tmp3 = self.dataSet.readNCVariable('RAINSH')
            prev = tmp1 + tmp2 + tmp3

        self.dataSet.setTimeIndex(currind)
        dt = currtime - prevtime
        if dt.seconds == 0:
            self.var = current - prev
        else:
            self.var = (current - prev)/(dt.seconds/3600.)
        self.var2 = self.var
        self.varTitle = "Rain Rate (mm hr$^{-1}$)\n" + self.dataSet.getTime() 


    def pwat(self):

        self.var = wrf.get_pwat(self.qvapor, self.rho, self.height)
        self.var2 = wrf.get_mslp(self.height, self.press, self.temp, self.qvapor)
        self.varTitle = "Total Precipitable Water (in)\n" + self.dataSet.getTime()

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

    def t2m(self):
    
        t2 = self.dataSet.readNCVariable('T2')
        self.var = wrf.convertT_KtoF(t2)
        self.var2 = self.var
        self.varTitle = "2-meter Temperature ($^{\circ}$ F)\n"+ self.dataSet.getTime()

    def td2m(self):

        t2 = self.dataSet.readNCVariable('T2')
        q2 = self.dataSet.readNCVariable('Q2')
        psfc = self.dataSet.readNCVariable('PSFC')
        var1 = wrf.get_td(t2, q2, psfc)
        self.var = wrf.convertT_KtoF(var1)
        self.var2 = self.var
        self.varTitle = "2-meter Dewpoint Temperature ($^{\circ}$ F)\n"+ self.dataSet.getTime()

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

    def wind10m(self):

        var1 = wrf.get_bulk_wind(self.u10, self.v10)
        #self.var = wrf.convertWind_MStoKT(var1)
        self.var = var1
        self.var2 = self.var
        #self.varTitle = "10-meter Wind (knots)\n"+ self.dataSet.getTime()
        self.varTitle = "10-meter Wind (m s$^{-1}$)\n"+ self.dataSet.getTime()

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

    def cape_SB(self):

        #Switched to Cython        
        #self.var, self.var2 = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'surface_based')
        height = wrf.unstaggerZ(self.height)
        self.var, self.var2 = wrf_cython.cape_sb(self.temp,self.qvapor,
                              self.press,height)        
        self.varTitle = "Surface-Based CAPE (J kg$^{-1}$)\n"+ self.dataSet.getTime()

    def cape_ML(self):

        #Switched to Cython
        #self.var, self.var2 = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'mixed_layer')
        height = wrf.unstaggerZ(self.height)
        self.var, self.var2 = wrf_cython.cape_ml(self.temp,self.qvapor,
                              self.press,height)
        self.varTitle = "Mixed Layer CAPE (J kg$^{-1}$)\n"+ self.dataSet.getTime()

    def cape_MU(self):
        #Switched to Cython
        #self.var, self.var2 = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'most_unstable')
        height = wrf.unstaggerZ(self.height)
        self.var, self.var2 = wrf_cython.cape_mu(self.temp,self.qvapor,
                              self.press,height)
        self.varTitle = "Most Unstable CAPE (J kg$^{-1}$)\n"+ self.dataSet.getTime()

    def cin_SB(self):
        #Switched to Cython
        #self.var2, self.var = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'surface_based')
        height = wrf.unstaggerZ(self.height)
        self.var2, self.var = np.array(wrf_cython.cape_sb(self.temp,self.qvapor,
                              self.press,height))        
        self.var *= -1
        self.varTitle = "Surface-Based CIN (J kg$^{-1}$)\n"+ self.dataSet.getTime()

    def cin_ML(self):
        #Sitched to Cython
        #self.var2, self.var = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'mixed_layer')
        height = wrf.unstaggerZ(self.height)
        self.var2, self.var = np.array(wrf_cython.cape_ml(self.temp,self.qvapor,
                              self.press,height))
        self.var *= -1
        self.varTitle = "Mixed Layer CIN (J kg$^{-1}$)\n"+ self.dataSet.getTime()

    def cin_MU(self):
        #Switched to Cython
        #self.var2, self.var = wrf.get_cape(self.temp,self.qvapor,
        #                      self.press,self.height,'most_unstable')
        height = wrf.unstaggerZ(self.height)
        self.var2, self.var = np.array(wrf_cython.cape_mu(self.temp,self.qvapor,
                              self.press,height))
        self.var *= -1
        self.varTitle = "Most Unstable CIN (J kg$^{-1}$)\n"+ self.dataSet.getTime()

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

    def cloud_albedo(self):
        
        #Read in radiation variables
        #Instantaneous downwelling flux at the top
        swdnt = self.dataSet.readNCVariable('SWDNT')
        #Instantaneous upwelling flux at the top - clear sky
        swuptc = self.dataSet.readNCVariable('SWUPTC')
        #Instantaneous upwelling flux at the top
        swupt = self.dataSet.readNCVariable('SWUPT')

        #Calculate cloud albedo
        self.var = (swupt - swuptc) / swdnt
        self.var2 = self.var
        self.varTitle = "Cloud Albedo \n"+ self.dataSet.getTime() 

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

