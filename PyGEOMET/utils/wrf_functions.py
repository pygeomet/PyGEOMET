#!/rhome/freitagb/anaconda3/bin/python

#import libraries required for the various wrf functions
import numpy as np
import csv
import matplotlib.path as mplPath
from ftplib import FTP
from scipy.special import gamma
import os

###############################################################################
#                                                                             #
#                            Begin Derived Functions                          #
# 1. GetSFCObs                                                                #
# 2. get_bulk_wind                                                            #
# 3. get_cape    (Written by Andrew White - whiteat@nsstc.uah.edu)            #
# 4. get_cl_albedo (Written by Andrew White - whiteat@nsstc.uah.edu)          #
# 5. get_froude                                                               #
# 6. get_height                                                               #
# 7. get_lcl                                                                  #
# 8. get_lwp                                                                  #
# 9. get_mslp                                                                 #
# 10. get_precip                                                              #
# 11. get_press                                                               #
# 12. get_pwat                                                                #
# 13. get_refl                                                                #
# 14. get_rh                                                                  #
# 15. get_rho                                                                 #
# 16. get_shear                                                               #
# 17. get_temp                                                                #
# 18. get_td                                                                  #
# 19. get_thetae (Written by Andrew White - whiteat@nsstc.uah.edu)            #
# 20. get_wc (Written by Andrew White - whiteat@nsstc.uah.edu)                #
# 21. mean_layer                                                              #
# 22. pot_vort (Written by Andrew White - whiteat@nsstc.uah.edu)              #
# 23. rel_vort                                                                #
#                                                                             #
###############################################################################

#####################  1.Begin function of GetSFCObs()   ######################
## Required libraries: ftplib, csv, os, numpy, matplotlib.path                #
##                                                                            #
## Inputs: xx = list of longitudes for all four grid corners. First and last  #
##         longitude point are the same (i.e. list is 5 elements long).       #
##                                                                            #
##         yy = list of latitudes for all four grid corners. First and last   #
##         latitude point are the same (i.e. list is 5 elements long).        #
##                                                                            #
##         sdate = the start date of the simulation period from the namelist  #
##                                                                            #
##         edate = the end date of the simulation period from the namelist    #
##                                                                            #
###############################################################################

def GetSFCObs(xx, yy, sdate, edate):
    #connect to the NCDC FTP site to get the AWOS station list
    #NOTE: the CSV file was easier to manage than TXT.
    if not os.path.isfile('/rhome/freitagb/wrf_seus/isd-history.csv'):
        ftp = FTP('ftp.ncdc.noaa.gov')
        ftp.login(user ='anonymous', passwd = 'freitagb@nsstc.uah.edu')
        ftp.cwd('pub/data/noaa/')
        ftp.retrbinary('RETR ' + 'isd-history.csv',\
open('~/wrf_seus/isd-history.csv','wb').write)
    
    #set up the arrays of values from the CSV file
    USAF = []; WBAN = []; STATION = []; NAME = []; CTRY = []; ST = []
    CALL = []; LAT = []; LON = []; ELEV = []; BEGIN = []; END = []
    
    #set up the fieldnames for the CSV reader
    fnames = ['USAF', 'WBAN', 'STATION NAME', 'CTRY', 'STATE', 'ICAO',\
'LAT', 'LON', 'ELEV(M)', 'BEGIN', 'END']

    with open('isd-history.csv') as csvfile:
        #set up the dictionary reader with a comma delimiter and 
        #  the given fieldnames
        reader = csv.DictReader(csvfile, fieldnames=fnames, delimiter=',')

        #this is used to skip the header lines. number of lines will 
        #  change with different files
        for i in range(0,18):
            next(reader)
        #for each row beyond the header, read in the data from the field name
        for row in reader:
            if row['LAT'] and row['LON']:
                USAF.append(row['USAF'])
                WBAN.append(row['WBAN'])
                STATION.append(row['STATION NAME'])
                CTRY.append(row['CTRY'])
                ST.append(row['STATE'])
                CALL.append(row['ICAO'])
                LAT.append(float(row['LAT']))
                LON.append(float(row['LON']))
                ELEV.append(row['ELEV(M)'])
                BEGIN.append(row['BEGIN'])
                END.append(row['END'])
    
    #unless you're moved to keep this file, delete the station list
    #os.remove('isd-history.csv')
    loc = np.array((xx,yy),dtype=float)
    loc = loc.T
    #set up a polygon using the input grid corner points in lat/lon coordinates
    gpath = mplPath.Path(loc)
     
    #set up the count and arrays of usaf code, lat, and lon for
    # sites within the grid domain
    count = 0
    usaf_good = []
    lat_good = []
    lon_good = []
    beg_good = []
    end_good = []

    #loop through the lat/lon values in the data,
    #  check if the lat/lon values fall within the grid and the 
    #  begin and end dates of the simulation are contained by 
    #  the observations available at the NCDC site
    for i in range(0,len(LON)):
        location = [LON[i],LAT[i]]
        if gpath.contains_point(location):
            count += 1   
            lat_good.append(LAT[i])
            lon_good.append(LON[i])

    #return the number of observations, and the lon/lat of the observations
    #  so they can be plotted on the map
    return count, lon_good, lat_good

###############################################################################
#########################  End function GetSFCObs()  ##########################
###############################################################################

###################  2.Begin function of get_bulk_wind()   ####################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: u = ndarray of the u-component wind values (Units: m/s)            #
##                                                                            #
##         v = ndarray of the v-component wind values (Units: m/s)            #
##                                                                            #
###############################################################################

def get_bulk_wind(u, v):
    
    #get the bulk wind
    wind = (u*u + v*v)**0.5
    
    #return the bulk wind value (Units: m/s)
    return wind

###############################################################################
#######################  End function get_bulk_wind()  ########################
###############################################################################

#####################  3.Begin function of get_cape()   #######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: u = ndarray of the u-component wind values (Units: m/s)            #
##                                                                            #
##         v = ndarray of the v-component wind values (Units: m/s)            #
##                                                                            #
###############################################################################

def get_cape(temp, qvapor, press, height, cape_type):

    #set constants
    g = 9.81
    rd = 287.04
    cp = 1004.
    Ep = 18.016/29.8

    #Determine array size
    dims = temp.shape
    nz = dims[0]
    ny = dims[1]
    nx = dims[2]
  
    #Unstagger Z
    Z = unstaggerZ(height)

    #convert pressure to hPa
    pressure = press*0.01

    #Calculate virtual temperature
    tempv = temp * (1. + 0.608 * qvapor)

    #Calculate potential temperature
    theta = temp * (1000./pressure)**(rd/cp)

    #Calculate dew point temperature
    vp = (qvapor*pressure) / (qvapor + Ep)
    td = (243.5 * np.log(vp/6.112)) / (17.67 - np.log(vp/6.112)) + 273.15

    if cape_type == 'surface_based': 
        #Calculate Saturation mixing ratio along constant potential temperature
        t_pt = theta[0,:,:] * (pressure/1000.)**(rd/cp) - 273.15
        vpsat = 6.112*np.exp(17.67*t_pt/(243.5+t_pt))
        qsat = Ep*vpsat / (pressure - vpsat)

       #Determine moist parcel path temperature and calculate CAPE and CIN
        t_par = np.zeros((nz,ny,nx))
        t_parv = np.zeros((nz,ny,nx))
        cape = np.zeros((ny,nx))
        cin = np.zeros((ny,nx))
        theta_lcl = theta[0,:,:]
        q_lcl = qvapor[0,:,:]

        for k in range(nz):
              jj_dry, ii_dry = np.where(qsat[k,:,:] > qvapor[0,:,:])
              jj_moist, ii_moist = np.where(qsat[k,:,:] <= qvapor[0,:,:])
              if (len(jj_dry) > 0):
                  t_par[k,jj_dry,ii_dry] = theta[0,jj_dry,ii_dry]\
                  *(pressure[k,jj_dry,ii_dry]/1000.)**(rd/cp)
                  t_parv[k,jj_dry,ii_dry] = t_par[k,jj_dry,ii_dry] \
                  *(1. + 0.608 * qvapor[0,jj_dry,ii_dry])
              if (len(jj_moist) > 0):
                  theta_lcl[jj_moist,ii_moist], q_lcl[jj_moist,ii_moist], t_par[k,jj_moist,ii_moist] \
                              = TParcel(pressure[k,jj_moist,ii_moist],q_lcl[jj_moist,ii_moist],
                                        theta_lcl[jj_moist,ii_moist])
                  t_parv[k,jj_moist,ii_moist] = t_par[k,jj_moist,ii_moist]\
                  * (1. + 0.608 * q_lcl[jj_moist,ii_moist])
              if (k > 0):
                  dz = Z[k,:,:]-Z[k-1,:,:]
                  dtparcel = (t_parv[k,:,:]+t_parv[k-1,:,:])/2.
                  dtenv = (tempv[k,:,:]+tempv[k-1,:,:])/2.
                  dum1 = g * dz * (dtparcel-dtenv) / dtenv
                  dum1[dum1 < 0] = 0.
                  cape += dum1
                  if (np.any(pressure[k-1,:,:] > pressure[0,:,:] - 300.)):
                      dum2 = g * dz * (dtparcel-dtenv) / dtenv
                      dum2[dum2 > 0] = 0.
                      cin += dum2


    if cape_type == 'mixed_layer':
        #Calculate mean layer values
        minP = pressure[0,:,:] - 100.
        diff = pressure - minP
        diff[diff < 0] -= 1000.
        diff2 = diff[1:,:,:] - diff[0:-1,:,:]
        diff2 = MatrixCondition(diff2)
        kk, jj, ii = np.where(diff2 == 1)
        total_temp = np.zeros((ny,nx))
        total_the = np.zeros((ny,nx))
        total_qv = np.zeros((ny,nx))
        count = np.zeros((ny,nx))

        for k in range(max(kk)):
            index = np.where(kk >= k)
            total_temp[jj[index],ii[index]] += np.squeeze(temp[k,jj[index],ii[index]])
            total_qv[jj[index],ii[index]] += np.squeeze(qvapor[k,jj[index],ii[index]])
            total_the[jj[index],ii[index]] += np.squeeze(theta[k,jj[index],ii[index]])
            count[jj[index],ii[index]] += 1

        #Mean Layer temp, mixing ratio and potential temp
        t_ml = total_temp / count
        q_ml = total_qv / count
        theta_ml = total_the / count

        #Calculate Saturation mixing ratio along constant potential temperature
        t_pt = theta_ml * (pressure/1000.)**(rd/cp) - 273.15
        vpsat = 6.112*np.exp(17.67*t_pt/(243.5+t_pt))
        qsat = Ep*vpsat / (pressure - vpsat)

        #Determine moist parcel path temperature and calculate CAPE and CIN
        t_par = np.zeros((nz,ny,nx))
        t_parv = np.zeros((nz,ny,nx))
        cape = np.zeros((ny,nx))
        cin = np.zeros((ny,nx))
        theta_lcl = theta_ml
        q_lcl = q_ml

        for k in range(nz):
              jj_dry, ii_dry = np.where(qsat[k,:,:] > q_ml)
              jj_moist, ii_moist = np.where(qsat[k,:,:] <= q_ml)
              if (len(jj_dry) > 0):
                  t_par[k,jj_dry,ii_dry] = theta_ml[jj_dry,ii_dry] * (pressure[k,jj_dry,ii_dry]/1000.)**(rd/cp)
                  t_parv[k,jj_dry,ii_dry] = t_par[k,jj_dry,ii_dry] * (1. + 0.608 * q_ml[jj_dry,ii_dry])
              if (len(jj_moist) > 0):
                  theta_lcl[jj_moist,ii_moist], q_lcl[jj_moist,ii_moist], t_par[k,jj_moist,ii_moist] \
                              = TParcel(pressure[k,jj_moist,ii_moist],q_lcl[jj_moist,ii_moist],theta_lcl[jj_moist,ii_moist])
                  t_parv[k,jj_moist,ii_moist] = t_par[k,jj_moist,ii_moist] * (1. + 0.608 * q_lcl[jj_moist,ii_moist])
              if (k > 0):
                  dz = Z[k,:,:]-Z[k-1,:,:]
                  dtparcel = (t_parv[k,:,:]+t_parv[k-1,:,:])/2.
                  dtenv = (tempv[k,:,:]+tempv[k-1,:,:])/2.
                  dum1 = g * dz * (dtparcel-dtenv) / dtenv
                  dum1[dum1 < 0] = 0.
                  cape += dum1
                  if (np.any(pressure[k-1,:,:] > pressure[0,:,:] - 300.)):
                      dum2 = g * dz * (dtparcel-dtenv) / dtenv
                      dum2[dum2 > 0] = 0.
                      cin += dum2

    if cape_type == 'most_unstable':
       #Determine the most unstable layer (sfc to sfc-300mb) -- max thetae
        minP = pressure[0,:,:] - 300.
        diff = pressure - minP
        diff[diff < 0] -= 1000.
        diff2 = diff[1:,:,:] - diff[0:-1,:,:]
        diff2 = MatrixCondition(diff2)
        kk, jj, ii = np.where(diff2 == 1)
        thetae = np.zeros((nz,ny,nx))
        thetae_max = np.zeros((ny,nx))
        theta_mu = np.zeros((ny,nx))
        t_mu = np.zeros((ny,nx))
        q_mu = np.zeros((ny,nx))

        for k in range(max(kk)):
            index = np.where(kk >= k)
            vp_sfc = 6.112*np.exp(17.67*(td[k,jj[index],ii[index]]-273.15)\
                   /(243.5+(td[k,jj[index],ii[index]]-273.15)))
            t_lcl = (2840./(3.5*np.log(temp[k,jj[index],ii[index]])-np.log(vp_sfc)-4.805))+55.
            thetae[k,jj[index],ii[index]] = (theta[k,jj[index],ii[index]]\
                                          *np.exp((3.376/t_lcl - 0.00254)*1000.\
                                          *qvapor[k,jj[index],ii[index]]*(1.+0.81*qvapor[k,jj[index],ii[index]])))

            if k > 0:
                i = np.where(thetae_max[jj[index],ii[index]] < thetae[k,jj[index],ii[index]])
                theta_mu[(jj[index])[i],(ii[index])[i]] = theta[k,(jj[index])[i],(ii[index])[i]]
                t_mu[(jj[index])[i],(ii[index])[i]] = temp[k,(jj[index])[i],(ii[index])[i]]
                q_mu[(jj[index])[i],(ii[index])[i]] = qvapor[k,(jj[index])[i],(ii[index])[i]]
                thetae_max[jj[index],ii[index]] = np.maximum(thetae_max[jj[index],ii[index]],thetae[k,jj[index],ii[index]])
            else:
                thetae_max[jj[index],ii[index]] = thetae[k,jj[index],ii[index]]
                theta_mu[jj[index],ii[index]] = theta[k,jj[index],ii[index]]
                t_mu[jj[index],ii[index]] = temp[k,jj[index],ii[index]]
                q_mu[jj[index],ii[index]] = qvapor[k,jj[index],ii[index]]

       #Calculate Saturation mixing ratio along constant potential temperature
        t_pt = theta_mu * (pressure/1000.)**(rd/cp) - 273.15
        vpsat = 6.112*np.exp(17.67*t_pt/(243.5+t_pt))
        qsat = Ep*vpsat / (pressure - vpsat)

        #Determine moist parcel path temperature and calculate CAPE and CIN
        t_par = np.zeros((nz,ny,nx))
        t_parv = np.zeros((nz,ny,nx))
        cape = np.zeros((ny,nx))
        cin = np.zeros((ny,nx))
        theta_lcl = theta_mu
        q_lcl = q_mu

        for k in range(nz):
              jj_dry, ii_dry = np.where(qsat[k,:,:] > q_mu)
              jj_moist, ii_moist = np.where(qsat[k,:,:] <= q_mu)
              if (len(jj_dry) > 0):
                  t_par[k,jj_dry,ii_dry] = theta_mu[jj_dry,ii_dry] \
                  * (pressure[k,jj_dry,ii_dry]/1000.)**(rd/cp)
                  t_parv[k,jj_dry,ii_dry] = t_par[k,jj_dry,ii_dry] \
                  * (1. + 0.608 * q_mu[jj_dry,ii_dry])
              if (len(jj_moist) > 0):
                  theta_lcl[jj_moist,ii_moist], q_lcl[jj_moist,ii_moist], t_par[k,jj_moist,ii_moist] \
                              = TParcel(pressure[k,jj_moist,ii_moist],q_lcl[jj_moist,ii_moist],theta_lcl[jj_moist,ii_moist])
                  t_parv[k,jj_moist,ii_moist] = t_par[k,jj_moist,ii_moist] * (1. + 0.608 * q_lcl[jj_moist,ii_moist])
              if (k > 0):
                  dz = Z[k,:,:]-Z[k-1,:,:]
                  dtparcel = (t_parv[k,:,:]+t_parv[k-1,:,:])/2.
                  dtenv = (tempv[k,:,:]+tempv[k-1,:,:])/2.
                  dum1 = g * dz * (dtparcel-dtenv) / dtenv
                  dum1[dum1 < 0] = 0.
                  cape += dum1
                  if (np.any(pressure[k-1,:,:] > pressure[0,:,:] - 300.)):
                      dum2 = g * dz * (dtparcel-dtenv) / dtenv
                      dum2[dum2 > 0] = 0.
                      cin += dum2

    #return cape/cin (J kg^-1)
    return cape, cin

###############################################################################
##########################  End function get_cape()  ##########################
###############################################################################


####################  4.Begin function of get_cl_albedo()   ###################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: swu = ndarray of upwelling solar radiation at TOA (Units: W/m2)    #
##                                                                            #
##         swucs = ndarray of clear-sky upwelling solar radiation at TOA      #
##                 (Units: W/m2)                                              #
##                                                                            #
##         swd = ndarray of downwelling solar radiation at the surface        #
##                 (Units: W/m2)                                              #
##                                                                            #
###############################################################################

def get_cl_albedo(swu, swucs, swd):

    alb = (swu/swucs) / swd

    #return the height (Units: None)
    return alb

###############################################################################
########################  End function get_cl_albedo()  #######################
###############################################################################


####################  5.Begin function of get_froude()   ######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: height = ndarray of geopotential height (Units: m)                 #
##                                                                            #
##         theta = ndarray of potential temperature (Units: K)                #
##                                                                            #
##         wind = ndarray of bulk wind (Units: m s-1)                         #
##                                                                            #
###############################################################################

def get_froude(height, theta, wind):

    #set the constants
    g = 9.81 # gravity (units: m s^-2)

    # get the difference in potential temperature from the 
    # surface to the mountain top
    dtheta = np.gradient(theta)[-2]
    dz = np.gradient(height)[-2]

    #calculate the Brunt-Vaisala Frequency
    N = ((g/theta)*(dtheta/dz))**(0.5)
    
    #calculate the Froude number
    fr = wind[0,:,:]/(height[0,:,:]*N[0,:,:])  
    #return the height (Units: None)
    return fr

###############################################################################
#########################  End function get_height()  #########################
###############################################################################


####################  6.Begin function of get_height()   ######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: ph = ndarray of the perturbation geopotential field (Units: m2/s2) #
##                                                                            #
##         phb = ndarray of the base geopotential field (Units: m2/s2)        #
##                                                                            #
###############################################################################

def get_height(ph, phb):

    #get the geopotential height
    height = (ph + phb)/9.81
    
    #return the height (Units: m)
    return height

###############################################################################
#########################  End function get_height()  #########################
###############################################################################

######################  7.Begin function of get_lcl()   #######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: t2 = ndarray of the 2-meter temperature values (Units: K)          #
##                                                                            #
##         rh = ndarray of the 2-meter relative humidity values (Units: %)    #
##                                                                            #
###############################################################################

def get_lcl(t2,rh):

    #get the lcl height
    lclhgt = (20 + (t2 - 273.15)/5.)*(100-rh)

    #return the lcl height (Units: m)
    return lclhgt

###############################################################################
##########################  End function get_lcl()  ###########################
###############################################################################

######################  8.Begin function of get_lwp()   #######################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: temp = ndarray of the temperature values (Units: K)                #
##                                                                            #
##         height = ndarray of the height values (Units: m)                   #
##                                                                            #
##         rho = ndarray of the density values (Units: kg/m3)                 #
##                                                                            #
##         qcloud = ndarray of the cloud mixing ratio values (Units: kg/kg)   #
##                                                                            #
##         qice = ndarray of the ice mixing ratio values (Units: kg/kg)       #
##                                                                            #
##         qsnow = ndarray of the snow mixing ratio values (Units: kg/kg)     #
##                                                                            #
###############################################################################
    
def get_lwp(temp, height, rho, qcloud, qice, qgraup, qsnow, qrain, vtype):

    #unstagger the height grid
    height = unstaggerZ(height)

    #get the composite mixing ratio
    if vtype == 'total':
        mr = qcloud + qrain + qice + qsnow + qgraup
    if vtype == 'liquid':
        mr = qcloud + qrain
    if vtype == 'ice':
        mr = qice + qsnow + qgraup
  
    #compute the vertical gradient in the height field
    dz = np.gradient(height)[-3]

    #compute the cloud total water in units of g m^-2
    cloud_LWP = np.sum(mr*rho*dz,axis=0)*1000. #g m^-2

    #return the cloud total water path (Units: g/m2)
    return cloud_LWP

###############################################################################
##########################  End function get_lwp()  ###########################
###############################################################################

######################  9.Begin function of get_mslp()   ######################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: h = ndarray of the height values (Units: m)                        #
##                                                                            #
##         p = ndarray of pressure values (Units: Pa)                         #
##                                                                            #
##         t = ndarray of temperature values (Units: K)                       #
##                                                                            #
##         q = ndarray of water vapor mixing ratio values (Units: kg/kg)      #
##                                                                            #
###############################################################################

def get_mslp(h, p, t, q):
   
    #unstagger the height grid
    h = unstaggerZ(h)
 
    #get the virtual temperature
    vtemp = t*(1+0.611*q)
    Rd = 287.05 # J kg^-1 K^-1 (dry-air gas constant)

    #reduce the temperature moist adiabatically from the surface height to mean sea level
    tmp = vtemp[0,:,:] + (h[0,:,:])*6.5/1000.
    vtemp[0,:,:] = (vtemp[0,:,:]+tmp)/2.0
    
    #compute the mean sea level pressure
    mslp = p[0,:,:]*np.exp((9.81*h[0,:,:])/(Rd*vtemp[0,:,:]))/100. #hPa
    
    #return the mean sea level pressure (Units: hPa)
    return mslp

###############################################################################
#########################  End function get_mslp()  ###########################
###############################################################################

####################  10.Begin function of get_precip()   #####################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: rc = ndarray of sub-grid scale precipitation values (Units: mm)    #
##                                                                            #
##         rnc = ndarray of grid scale precipitation values (Units: mm)       #
##                                                                            #
##         sh = ndarray of shallow cumulus precipitation values (Units: mm)   #
##                                                                            #
###############################################################################

def get_precip(rc, rnc, sh):

    #get the accumulated precipitation
    pcp = rc + rnc + sh
    
    #return the precipitation array (Units: mm)
    return pcp

###############################################################################
#########################  End function get_precip()  #########################
###############################################################################

####################  11.Begin function of get_press()    #####################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: p = ndarray of the perturbation pressure field values (Units: Pa)  #
##                                                                            #
##         pb = ndarray of the base pressure field values (Units: Pa)         #
##                                                                            #
###############################################################################
 
def get_press(p, pb):

    #get the 3D pressure array
    press = pb+p
    
    #return the pressure (Units: Pa)
    return press

###############################################################################
########################  End function get_press()  ###########################
###############################################################################

##################### 12.Begin function of get_pwat()   #######################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: qv = ndarray of water vapor mixing ratio values (Units: kg/kg)     #
##                                                                            #
##         rho = ndarray of atmospherice density values (Units: kg/m3)        #
##                                                                            #
##         height = ndarray of height values (Units: m)                       #
##                                                                            #
###############################################################################

def get_pwat(qv, rho, height):

    #unstagger the height grid
    height = unstaggerZ(height)

    #compute the vertical gradient in the height field
    dz = np.gradient(height)[-3]
    
    #convert m to in
    mtoin = 1000./25.4
    #convert the water vapor to inches (mm/in)    (kg^-1 m3)
    pwat = np.sum(qv*rho*dz,axis=0)*mtoin/1000.
     
     #return precipitable water (Units: inches)
    return pwat

###############################################################################
#########################  End function get_pwat()  ###########################
###############################################################################

###################### 13.Begin function of get_refl()   ######################
## Required libraries: numpy,scipy.special                                    #
##                                                                            #
## Inputs: qrain = ndarray of rain water mixing ratio (Units: kg/kg)          #
##                                                                            #
##         qgraup = ndarray of graupel mixing ratio (Units: kg/kg)            #
##                                                                            #
##         qsnow = ndarray of snow mixing ratio (Units: kg/kg)                #
##                                                                            #
##         rho = ndarray of atmospheric density values (Units: kg/m3)         #
##                                                                            #
##         temp = ndarray of the temperature values (Units: K)                #
##                                                                            #
##         height = ndarray of the height values (Units: m)                   #
##                                                                            #
###############################################################################

def get_refl(qrain,qgraup,qsnow,rho,temp,height):

    #set up the constants
    expnt = 7./4.
    expnt2 = -0.75

    #the following constants are consistent with Thompson et al., 2008
    alpha = 0.224 # accounts for dielectric effects
    rho_r = 1000. # kg m^-3 -- density of rain
    rho_g = 400. # kg m^-3 -- density of graupel
    rho_s = 100. # kg m^-3 -- density of snow
    t0 = 273.15 # K -- triple point of water
    N1 =  1.e10 # m^-4 -- constant 1 for rain intercept parameter
    N2 = 8.e6 # m^-4 -- constant 2 for rain intercept parameter
    qr0 = 1.e-4 # kg kg^-1 -- reference rain water mixing ratio

    #if mixing ratios less than 0, make them 0
    qrain[np.where(qrain < 0)] = 0.0
    qsnow[np.where(qsnow < 0)] = 0.0
    qgraup[np.where(qgraup < 0)] = 0.0

    #set up the dummy arrays for the intercept parameters    
    dims = qrain.shape
    dum1 = np.zeros((dims[0],dims[1],dims[2]), dtype=np.float)
    dum1.fill(-0.001)
    dum2 = np.zeros((dims[0],dims[1],dims[2]), dtype=np.float)
    dum2.fill(2.e8)
    dum3 = np.zeros((dims[0],dims[1],dims[2]), dtype=np.float)
    dum3.fill(1.e4)
    dum4 = np.zeros((dims[0],dims[1],dims[2]), dtype=np.float)
    dum5 = np.zeros((dims[0],dims[1],dims[2]), dtype=np.float)
    dum5.fill(1.e10)

    #compute the minimum temperature difference
    temp_cel = np.minimum(dum1, temp-273.15)
    
    #compute the intercept parameter for snow
    Nos = np.minimum(dum2,2.e6*np.exp(-0.12*temp_cel))
   
    #compute the intercept parameter for graupel
    kk, jj, ii = np.where(qgraup > 1.e-15)
    dum4[kk,jj,ii] = 2.38*(np.pi*rho_g/(rho[kk,jj,ii]*qgraup[kk,jj,ii]))**0.92
    Nog = np.maximum(dum3, np.minimum(dum4, 5.e7))

    #compute the intercept parameter for rain
    kk, jj, ii = np.where(qrain > 1.e-15)
    dum5[kk,jj,ii] = ((N1-N2)/2.0)*np.tanh((qr0 - qrain[kk,jj,ii])/(0.25*qr0)) + (N1+N2)/2.0 # m^-4
    Nor = dum5    

    #Compute the star values for each precip phase
    Norstar = gamma(7)*pow(10,18)*(1/(np.pi*rho_r))**expnt
    Nogstar = gamma(7)*pow(10,18)*(1/(np.pi*rho_g))**expnt*(rho_g/rho_r)**2*alpha
    Nosstar = gamma(7)*pow(10,18)*(1/(np.pi*rho_s))**expnt*(rho_s/rho_r)**2*alpha

    #Compute the reflectivity values in units of mm^6 m^-3
    refl_mm = Norstar*(rho*qrain)**expnt*Nor**expnt2 + Nosstar*(rho*qsnow)**expnt*Nos**expnt2 + Nogstar*(rho*qgraup)**expnt*Nog**expnt2
   
    #set the minimum reflectivity to -30 dBz 
    dum1.fill(0.001)
    refl_mm = np.maximum(dum1, refl_mm)

    #convert reflectivity to dBz
    refl_dbz  = 10.*np.log10(refl_mm)
 
    #get the index 1-km AGL
    diff = abs((height) - (height[0,:,:]+1000.))
    zz = np.where(diff == diff.min())[0]
    #set up the composite reflectivity arrays
    comp_refl = np.amax(refl_dbz[0:zz,:,:], axis=0)   

    #return the composite reflectivity (Units: dBz)
    return comp_refl        

###############################################################################
##########################  End function get_refl()  ##########################
###############################################################################

###################### 14.Begin function of get_rh()   ########################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: t = ndarray of temperature values (Units: K)                       #
##                                                                            #
##         q = ndarray of water vapor mixing ratio values (Units: kg/kg)      #
##                                                                            #
##         p = ndarray of pressure values (Units: Pa)                         #
##                                                                            #
###############################################################################

def get_rh(t, q, p):

    #calculate the saturation vapor pressure - Clausius Clapeyron Equation (Bolton, 1980)
    es = 6.112 * np.exp(17.67*(t-273.15) / (t-273.15+243.5))

    #calculate the environmental vapor pressure - Clausius Clapeyron Equation
    e = (q*p/0.622)/(1+(q/0.622))/100.

    #calculate the relative humidity
    rh = e/es*100.

    #return relative humidity (Units: %)
    return rh

###############################################################################
##########################  End function get_rh()  ############################
###############################################################################

##################### 15.Begin function of get_rho()   ########################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: press = ndarray of pressure values (Units: Pa)                     #
##                                                                            #
##         temp = ndarray of temperature values (Units: K)                    #
##                                                                            #
##         qvapor = ndarray of water vapor mixing ratio values (Units: kg/kg) #
##                                                                            #
###############################################################################

def get_rho(press, temp, qvapor):
    
    #if mixing ratio less than 0, make it 0
    qvapor[np.where(qvapor < 0)] = 0.0   
 
    #get the virtual temperature
    vtemp = temp*(1 + 0.611*qvapor)
    R = 287.05 # J kg^-1 K^-1 -- dry air gas constant
    
    #compute the 3D density field
    rho = press/(R*vtemp)
    
    #return the density array (Units: kg/m3)
    return rho

###############################################################################
##########################  End function get_rho()  ###########################
###############################################################################


#################### 16.Begin function of get_shear()   #######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: u = ndarray of u-wind values (Units: m/s)                          #
##                                                                            #
##         v = ndarray of v-wind values (Units: m/s)     )                    #
##                                                                            #
##         grid = ndarray of the vertical coordinate variable                 #
##                                                                            #
##         ref_val = reference level for computing shear                      #
##                                                                            #
###############################################################################

def get_shear(u, v, grid, ref_val):

    #unstagger the wrf u and v winds  
    ucorr = unstaggerX(u)
    vcorr = unstaggerY(v)

    #interpolate the u and v wind fields to the reference level
    u_lvl = linear_interpolate(ucorr, grid, ref_val)
    v_lvl = linear_interpolate(vcorr, grid, ref_val)
  
    #get the u and v wind fields at the surface
    u_sfc = ucorr[0,:,:]
    v_sfc = vcorr[0,:,:]

    #get the u and v components of the shear
    u_shear = u_lvl - u_sfc
    v_shear = v_lvl - v_sfc

    #get the total shear 
    shear = get_bulk_wind(u_shear,v_shear)
    #return the u and v components of shear and the
    # bulk shear (Units: knots)
    return u_shear, v_shear, shear

###############################################################################
#########################  End function get_shear()  ##########################
###############################################################################

##################### 17.Begin function of get_temp()   #######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: press = ndarray of pressure values (Units: Pa)                     #
##                                                                            #
##         pottemp = ndarray of potential temperature values (Units: K)       #
##                                                                            #
###############################################################################

def get_temp(press, pottemp):
    
    #add the constant to the perturbation potential temperature
    theta = pottemp + 300. # K
    
    #set up the constants
    Rd = 287.05 # J kg^-1 K^-1 -- dry air gas constant
    Cp = 1004.0 # J kg^-1 K^-1 -- specific heat of dry air
    
    #compute the 3D temperature field using Poisson's equation
    temp = theta*((100000./press)**(-1.*Rd/Cp))
    
    #return the temperature array (Units: K)
    return temp

###############################################################################
##########################  End function get_temp() ###########################
###############################################################################

###################### 18.Begin function of get_td()   ########################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: t = ndarray of temperature values (Units: K)                       #
##                                                                            #
##         q = ndarray of water vapor mixing ratio values (Units: kg/kg)      #
##                                                                            #
##         p = ndarray of pressure values (Units: Pa)                         #
##                                                                            #
###############################################################################

def get_td(t, q, p):
    
    #calculate the saturation vapor pressure - Clausius Clapeyron Equation (Bolton, 1980)
    es = 6.112 * np.exp(17.67*(t-273.15) / (t-273.15+243.5))

    #calculate the environmental vapor pressure - Clausius Clapeyron Equation
    e = (q*p/0.622)/(1+(q/0.622))/100.
    
    #calculate the relative humidity
    rh = e/es*100.
    
    #calculate the dew point temperature - Clausius Clapeyron Equation solved for T
    td = 243.04*(np.log(rh/100.) + ((17.625*(t-273.15))/(243.04 + t-273.15)))\
    /(17.625-np.log(rh/100.)-((17.625*(t-273.15))/(243.04+t-273.15)))

    #convert dew point temperature to K
    td +=  273.15
 
    #return the dew point temperature array (Units: K)
    return td

###############################################################################
#########################  End function get_td()  #############################
###############################################################################


###################### 19.Begin function of get_thetae()   ####################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: t = ndarray of temperature values (Units: K)                       #
##                                                                            #
##         q = ndarray of water vapor mixing ratio values (Units: kg/kg)      #
##                                                                            #
##         p = ndarray of pressure values (Units: Pa)                         #
##                                                                            #
###############################################################################

def get_thetae(t, q, p):

    #constants
    rd = 287.05
    cp = 1004.7

    #Calculate theta
    theta = t * (100000./p)**(rd/cp)

    #Surface dew point
    td = get_td(t, q, p)

    #Calculate vapor pressure
    vp = 6.112 * np.exp(17.67*(td-273.15)/(243.5+td-273.15))

    #Calculate lcl temperature
    t_lcl = (2840./(3.5*np.log(t)-np.log(vp)-4.805))+55.

    #Calculate surface theta-e
    theta_e = theta * np.exp((3.376/t_lcl - 0.00254)*1000.*q*(1. + 0.81 * q))

    #return theta-e at the value desired value (Units: K)
    return theta_e

###############################################################################
########################  End function get_thetae()  ##########################
###############################################################################


######################## 20.Begin function of get_wc()   ######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: qcloud = ndarray of cloud water mixing ratio (Units: kg/kg)        #
##                                                                            #
##         qrain = ndarray of rain water mixing ratio (Units: kg/kg)          #
##                                                                            #
##         qice = ndarray of ice water mixing ratio (Units: kg/kg)            #
##                                                                            #
##         qsnow = ndarray of snow mixing ratio (Units: kg/kg)                #
##                                                                            #
##         qgraup = ndarray of graupel mixing ratio (Units: kg/kg)            #
##                                                                            #
##         dtype = type of water content to be calculated (cloud_water,       #
##                 ice_water, total_water)                                    #
##                                                                            #
###############################################################################

def get_wc(qcloud, qrain, qice, qsnow, qgraup, dtype):

    # calculate cloud water content
    if dtype == 'cloud_water':
        wc = qcloud + qrain
    #calculate ice water content
    if dtype == 'ice_water':
        wc = qice + qsnow + qgraup

    #calculate total water content
    if dtype == 'total_water':
        wc = qcloud + qrain + qsnow + qice + qgraup

    #error check
    if dtype != 'cloud_water' or dtype != 'ice_water' or\
    dtype != 'total_water':
        print('This is not a valid option for this function')
    #return the desired water content (Units: kg kg-1)
    else:
        return wc
        
###############################################################################
########################  End function get_thetae()  ##########################
###############################################################################


#####################  4.Begin function of mean_layer()   #####################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: swu = ndarray of upwelling solar radiation at TOA (Units: W/m2)    #
##                                                                            #
##         swucs = ndarray of clear-sky upwelling solar radiation at TOA      #
##                 (Units: W/m2)                                              #
##                                                                            #
##         swd = ndarray of downwelling solar radiation at the surface        #
##                 (Units: W/m2)                                              #
##                                                                            #
###############################################################################

def mean_layer(var,hgt,ref1,ref2):

    #find where the height array matches the lower bound
    dims = var.shape
    var_ml = np.zeros((dims[1],dims[2]),dtype=np.float)
    diff = ref1-hgt
    diff2 = ref2-hgt
 
    for i in range(0,dims[1]):
        for j in range(0,dims[2]):
            lind = min(np.where(diff[:,i,j] < 0)[0])
            if lind == 0:
                lind1 = lind; lind2 = 1
            elif lind == (len(diff[:,i,j])-1):
                lind1 = ind-1; lind2 = lind
            else:
                lind1 = lind-1; lind2 = lind
            hind = min(np.where(diff2[:,i,j] < 0)[0])
            if hind == 0:
                hind1 = hind; hind2 = 1
            elif hind == (len(diff[:,i,j])-1):
                hind1 = hind-1; hind2 = hind
            else:
                hind1 = hind-1; hind2 = hind

            #calculate the mean layer value
            var_ml[i,j] = (((var[lind2,i,j] + var[lind,i,j])/2.) 
                        + ((var[hind2,i,j]+var[hind,i,j])/2.))/2.

    #return the mean layer value
    return var_ml

###############################################################################
########################  End function mean_layer()  ##########################
###############################################################################


##################### 22.Begin function of pot_vort()   #######################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: u = ndarray of u-component wind values (Units: m/s)                #
##                                                                            #
##         v = ndarray of v-component wind values (Units: m/s)                #
##                                                                            #
##         f = ndarray of Coriolis force (Units: s^-1)                        #
##                                                                            #
##         dx = float value representing E/W Grid Spacing (Units: m)          #
##                                                                            #
##         dy = float value representing N/S Grid Spacing (Units: m)          #
##                                                                            #
###############################################################################

def pot_vort(u,v,f,dx,dy,press,theta):
 
    g = 9.81 #m s^-2

    #Unstagger U and V
    u = unstaggerX(u)
    v = unstaggerY(v)

    #np.gradient returns the 3d gradients
    #Gradient in k is 0 index
    #Gradient in j is 1 index
    #Gradient in i is 2 index
    #Change in U
    du = np.gradient(u)
    #Change in V
    dv = np.gradient(v)
    #Change in Potential Temp
    dtheta = np.gradient(theta)
    #Change in Pressure
    dp = np.gradient(press)

    #Change in U in the y-direction
    dudy = du[1]/dy

    #Change in V in the x-direction
    dvdx = dv[2]/dx

    #Vertical gradient of potential temperature
    dthetadp = dtheta[0]/dp[0]
 
    #Calculate absolute vorticity
    av = dvdx - dudy + f

    #Calculate potential vorticity
    pv = -g * av * dthetadp * pow(10,6)

    #return the vorticity (Units: PVU [10^-6 m-2 s-1 K kg-1]) 
    return pv

###############################################################################
##########################  End function pot_vort()  ##########################
###############################################################################


##################### 23.Begin function of rel_vort()   #######################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: u = ndarray of u-component wind values (Units: m/s)                #
##                                                                            #
##         v = ndarray of v-component wind values (Units: m/s)                #
##                                                                            #
##         dx = float value representing E/W Grid Spacing (Units: m)          #
##                                                                            #
##         dy = float value representing N/S Grid Spacing (Units: m)          #
##                                                                            #
###############################################################################

def rel_vort(u,v,dx,dy):
    #get the gradient of the wind in the u- and v-directions
    du = np.gradient(u)
    dv = np.gradient(v)
    
    #compute the relative vorticity (units : 10^-5 s^-1)
    vort = ((dv[-1]/dx) - (du[-2]/dy))*pow(10,5) 

    #return the vorticity (Units: 10^-5 s-1) 
    return vort

###############################################################################
##########################  End function rel_vort()  ##########################
###############################################################################


################################################################################
#                                                                              #
#                            Begin Utility Functions                           #
#                                                                              #
# 1. convertT_KtoC                                                             #
# 2. convertT_KtoF                                                             #
# 3. convertT_CtoF                                                             #
# 4. convertT_FtoC                                                             #
# 5. convertP_MMtoIN                                                           #
# 6. convertWind_MStoKT                                                        #
# 7. convertWind_MStoMPH                                                       #
# 8. get_distance                                                              #
# 9. hypsometric                                                               #
# 10. linear_interpolate                                                       #
# 11. loglinear_interpolate                                                    #
# 12. MatrixCondition  (Written by Andrew White - whiteat@nsstc.uah.edu)       #
# 12. TParcel          (Written by Andrew White - whiteat@nsstc.uah.edu)       #
# 13. unstaggerX                                                               #
# 14. unstaggerY                                                               #
# 15. unstaggerZ                                                               #
#                                                                              #
################################################################################

####################   Begin function of convertT_KtoC()  #####################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: new_t = ndarray of temperature values(Units: K)                    #
##                                                                            #
###############################################################################

def convertT_KtoC(new_t):

    #convert temperature from Kelvin to Celsius
    new_t -= 273.15

    #return the converted temperature (Units: *C)
    return new_t

###############################################################################
######################  End function convertT_KtoC()  #########################
###############################################################################

####################   Begin function of convertT_KtoF  #######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: new_t = ndarray of temperature values (Units: K)                   #
##                                                                            #
###############################################################################

def convertT_KtoF(new_t):

    #convert temperature from Kelvin to Celsius
    new_t -= 273.15
    new_t = (new_t*1.8) + 32.

    #return the converted temperature (Units: *F)
    return new_t

###############################################################################
######################  End function convertT_KtoF()  #########################
###############################################################################

####################   Begin function of convertT_CtoF  #######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: new_t = ndarray of temperature values (Units: *C)                  #
##                                                                            #
###############################################################################

def convertT_CtoF(new_t):

    #convert temperature from Kelvin to Celsius
    new_t = (new_t*1.8) + 32.

    #return the converted temperature (Units *F)
    return new_t

###############################################################################
######################  End function convertT_CtoF()  #########################
###############################################################################

###################   Begin function of convertT_FtoC  ########################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: new_t = ndarray of temperature values (Units: *F)                  #
##                                                                            #
###############################################################################

def convertT_FtoC(new_t):

    #convert temperature from Kelvin to Celsius
    new_t = (new_t - 32.)/1.8

    #return the converted temperature (Units *C)
    return new_t

###############################################################################
######################  End function convertT_FtoC()  #########################
###############################################################################

###################  Begin function of convertP_MMtoIN  #######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: prec = ndarray of precipitation values (Units: mm)                 #
##                                                                            #
###############################################################################

def convertP_MMtoIN(prec):

    #convert precipitation from mm to inches
    prec /= 25.4

    #return converted precipitation (Units: inches)
    return prec

###############################################################################
#####################  End function convertP_MMtoIN()  ########################
###############################################################################

##################  Begin function of convertWind_MStoKT  #####################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: wind = ndarray of wind values(Units: m/s)                          #
##                                                                            #
###############################################################################

def convertWind_MStoKT(wind):

    #convert wind from m s^-1 to Knots
    wind *= 1.943844492574 #1.943844492574 knots per m/s

    #return the wind field (Units: Knots)
    return wind

###############################################################################
####################  End function convertWind_MStoKT()  ######################
###############################################################################

##################  Begin function of convertWind_MStoMPH  ####################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: wind = ndarray of wind values (Units: m/s)                         #
##                                                                            #
###############################################################################

def convertWind_MStoMPH(wind):

    #convert wind from m s^-1 to MPH
    wind *= (3600./1609.334) #1609.334 meters per mile

    #return the wind field (Units: MPH)
    return wind

###############################################################################
###################  End function convertWind_MStoMPH()  ######################
###############################################################################

####################  Begin function of get_distance  #########################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: lat1 = ndarray of latitude points (Units: *)                       #
##                                                                            #
##         lon1 = ndarray of longitude points (Units: *)                      #
##                                                                            #
##         lat2 = ndarray of latitude points (Units: *)                       #
##                                                                            #
##         lon2 = ndarray of longitude points (Units: * )                     #
##                                                                            #
###############################################################################

def get_distance(lat1,lon1,lat2,lon2):

    #set up the constants for the distance calculation
    r = 6371000. #Earth's radius (units: meters)

    #convert the latitudes into radians
    l1 = lat1*(np.pi/180.)    
    l2 = lat2*(np.pi/180.)

    #get the difference in longitudes in radians
    dlon = (abs(lon2) - abs(lon1))*(np.pi/180.)

    #calculate the distance using the spherical law of cosines
    dist = np.arccos((np.sin(l1)*np.sin(l2)) + \
    (np.cos(l1)*np.cos(l2)*np.cos(dlon)))*r

    #return the distance between the two latitudes (Units: km)
    return dist/1000.

###############################################################################
#######################  End function get_distance()  #########################
###############################################################################

#####################  Begin function of hypsometric  #########################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: indata = ndarray of the data that will be interpolated.            #
##                  Typically this will be height (Units: m)                  #
##                                                                            #
##         gdata = ndarray of the grid data.                                  #
##                 Typically this will be pressure (Units: Pa)                #
##                                                                            #
##         ref_val = reference value within the gdata that indata will be     #
##                   interpolated to. Typically this will be pressure         #
##                   (Units: Pa)                                              #
##                                                                            #
##         t = ndarray of temperature values (Units: K)                       #
##                                                                            #
###############################################################################

def hypsometric(indata, gdata, ref_val, t):

    dims = indata.shape
#    diff = np.zeros(dims[0], dtype=np.float)
#    ref = np.zeros((dims[0],dims[1],dims[2]),dtype=np.float)
#    ref.fill(ref_val)
#    t_avg = np.zeros((dims[1],dims[2]),dtype=np.float) 
    newvar = np.zeros((dims[1],dims[2]),dtype=np.float)
    R = 287.05 # dry air gas constant J kg^-1 K^-1
    g = 9.81 # gravity m s^-2
    diff = gdata - ref_val
    for i in range(0, dims[1]):
        for j in range(0, dims[2]):
#    diff = gdata-ref
#    print(diff.shape)
#    kk, jj, ii = np.where(np.sign(diff) == -1)
#    print(kk.shape)
#    kk = kk.min()
#    t_avg[jj,ii] = (t[kk-1,jj,ii] + t[kk,jj,ii])/2.
#    newvar[jj,ii] = indata[kk-1,jj,ii] + ((R*t_avg[jj,ii])/g)*np.log(gdata[kk-1,jj,ii]/ref_val)

#            for sgn in range(0, dims[0]-1):
#                if np.sign(diff[sgn]) != np.sign(diff[sgn+1]):
#                    ind1 = sgn; ind2 = sgn+1
            ind = min(np.where(diff[:,i,j] < 0)[0])
            if ind == 0:
                ind1 = ind; ind2 = 1
            elif ind == (len(diff[:,i,j])-1):
                ind1 = ind-1; ind2 = ind
            else:
                ind1 = ind-1; ind2 = ind

            #get the layer averaged temperature
            t_avg = (t[ind1,i,j] + t[ind2,i,j])/2.
            newvar[i,j] = indata[ind1,i,j] + ((R*t_avg)/g)*np.log(gdata[ind1,i,j]/ref_val)
            if i == 130 and j == 87:
                print(diff[:,j,i],gdata[:,j,i])
                print(ref_val,t_avg,ind)
                print(newvar[j,i])
    # return the height value interpolated to the reference level (Units: m)
    return newvar

###############################################################################
#######################  End function hypsometric()  ##########################
###############################################################################



################  Begin function of linear_interpolate  ####################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: indata = ndarray of the data that will be interpolated.            #
##                                                                            #
##         gdata = ndarray of the grid data.                                  #
##                                                                            #
##         ref_val = reference value within the gdata that indata will be     #
##                   interpolated to. Typically this will be pressure         #
##                   (Units: Pa)                                              #
##                                                                            #
###############################################################################

def linear_interpolate(indata,gdata,ref_val):

    dims = indata.shape
    newvar = np.zeros((dims[1],dims[2]),dtype=np.float)
    diff = gdata - ref_val
    for i in range(0, dims[1]):
        for j in range(0, dims[2]):
            ind = min(np.where(diff[:,i,j] > 0)[0])
            if ind == 0:
                ind1 = ind; ind2 = 1
            elif ind == (len(diff[:,i,j])-1):
                ind1 = ind-1; ind2 = ind
            else:
                ind1 = ind-1; ind2 = ind
#            for sgn in range(0,len(diff)-1):
#                if np.sign(diff[sgn]) != np.sign(diff[sgn+1]):
#                    ind1 = sgn; ind2 = sgn+1
            m = (indata[ind2,i,j]-indata[ind1,i,j])/(gdata[ind2,i,j]
            -gdata[ind1,i,j])
            c = indata[ind2,i,j] - (m*gdata[ind2,i,j])
            newvar[i,j] = m*ref_val + c

    # return the value interpolated to the reference level
    return newvar

###############################################################################
####################  End function linear_interpolate()  ######################
###############################################################################



################  Begin function of loglinear_interpolate  ####################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: indata = ndarray of the data that will be interpolated.            #
##                                                                            #
##         gdata = ndarray of the grid data.                                  #
##                                                                            #
##         ref_val = reference value within the gdata that indata will be     #
##                   interpolated to. Typically this will be pressure         #
##                   (Units: Pa)                                              #
##                                                                            #
###############################################################################

def loglinear_interpolate(indata,gdata,ref_val):

    dims = indata.shape
    newvar = np.zeros((dims[1],dims[2]),dtype=np.float)
    diff = gdata - ref_val
    for i in range(0, dims[1]):
        for j in range(0, dims[2]):
            ind = min(np.where(diff[:,i,j] < 0)[0])
            if ind == 0:
                ind1 = ind; ind2 = 1
            elif ind == (len(diff[:,i,j])-1):
                ind1 = ind-1; ind2 = ind
            else:
                ind1 = ind-1; ind2 = ind
#            for sgn in range(0,len(diff)-1):
#                if np.sign(diff[sgn]) != np.sign(diff[sgn+1]):
#                    ind1 = sgn; ind2 = sgn+1
            m = (indata[ind2,i,j]-indata[ind1,i,j])/(np.log(gdata[ind2,i,j])
            -np.log(gdata[ind1,i,j]))
            c = indata[ind2,i,j] - (m*np.log(gdata[ind2,i,j]))
            newvar[i,j] = m*np.log(ref_val) + c

    # return the value interpolated to the reference level
    return newvar

###############################################################################
###################  End function loglinear_interpolate()  ####################
###############################################################################



####################  Begin function of MatrixCondition() #####################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: indata = input matrix for evaluation                               #
##                                                                            #
###############################################################################

def MatrixCondition(a):

    b = abs(a).max(0)
    condition = abs(a) == b[np.newaxis,...]
    a[condition] = 1
    a[np.logical_not(condition)] = 0
    return a

###############################################################################
######################  End function MatrixCondition()  #######################
###############################################################################



#######################  Begin function of TParcel  ###########################
## Required libraries: numpy                                                  #
##                                                                            #
## Inputs: p = ndarray of the pressure field at a model level. (units: hPa)   #
##                                                                            #
##         w = water vapor mixing ratio at a model level. (units: kg kg^-1)   #
##                                                                            #
##         theta = potential temperature value at a model level (units: K)    #
##                                                                            #
###############################################################################

def TParcel(p, w, theta):

    dims = p.shape
    rd = 287.05
    cp = 1004.7
    kpa = rd/cp
    ep = 18.016/29.87
    L = 2.5e6
    p0 = 1000.
    dwthres = 0.0000001
    dwf = np.zeros(dims[0]) + 1./3.
    tk = theta * (p/p0)**kpa
    vp = 6.112*np.exp(17.67*(tk-273.15)/(243.5+(tk-273.15)))
    ws = ep*vp / (p - vp)
    dw = dwf * (ws - w)

    while(np.any(abs(dw) > dwthres)):

        #Set old values for fail safe
        dwold = dw
        theta_old = theta
        w_old = w
        tk_old = tk
        ws_old = ws
        #Calculate new values
        dw = dwf * (ws - w)
        theta = theta - (dw*L*theta) / (cp*tk)
        w = w + dw
        tk = theta * (p/p0)**kpa
        vp = 6.112*np.exp(17.67*(tk-273.15)/(243.5+(tk-273.15)))
        ws = ep*vp / (p - vp)
        #Check to make sure ws is less than the guess w
        #If so, reduce the iterator and recalculate
        if np.any(ws > w):
            dwf = 0.5 * dwf
            theta = theta_old
            w = w_old
            tk = tk_old
            ws = ws_old

    tv = (theta*(p/p0)**kpa)
    return theta, w, tv

###############################################################################
##########################  End function Tparcel()  ###########################
###############################################################################



########################  Begin function of unstaggerX  #######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: var = ndarray to be unstaggered                                    #
##                                                                            #
###############################################################################

def unstaggerX(var):
    
    #unstagger the variable in the horizontal x-dimension
    dims = var.shape
    newvar = (var[:,:,1:dims[2]] + var[:,:,0:dims[2]-1])/2.
    
    #return the unstaggered variable
    return newvar

###############################################################################
#######################  End function unstaggerX()  ###########################
###############################################################################

########################  Begin function of unstaggerY  #######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: var = ndarray to be unstaggered                                    #
##                                                                            #
###############################################################################

def unstaggerY(var):
    
    #unstagger the variable in the horizontal y-dimension
    dims = var.shape
    newvar = (var[:,1:dims[1],:] + var[:,0:dims[1]-1,:])/2.
    
    #return the unstaggered variable
    return newvar

###############################################################################
#######################  End function unstaggerY()  ###########################
###############################################################################

#####################  Begin function of unstaggerZ()  ########################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: var = ndarray to be unstaggered                                    #
##                                                                            #
###############################################################################

def unstaggerZ(var):
    
    #unstagger the variable in the vertical dimension
    dims = var.shape
    newvar = (var[1:dims[0],:,:] + var[0:dims[0]-1,:,:])/2.
    
    #return the unstaggered variable
    return newvar

###############################################################################
#######################  End function unstaggerZ()  ###########################
###############################################################################

################################################################################
#                                                                              #
#                            Begin Misc Functions                              #
#                                                                              #
# 1. getDvarList                                                               #
#                                                                              #
################################################################################

#####################  Begin function of getDvarList()  #######################
## Required libraries: none                                                   #
##                                                                            #
## Inputs: varlist = list of available variabls                               #
##                                                                            #
###############################################################################

def getDvarList(varlist):
    dvarlist = []
    variables = varlist

    dvarlist.extend(['CAPE_SB','CAPE_ML','CAPE_MU',
                     'CIN_SB','CIN_ML','CIN_MU',
                     'mslp','lcl_hgt','pwat', 't2m','td2m',
                     'temp-3d','temp_ML','theta-e-3d',
                     'wind10m'])

    if ('U' in variables and 'V' in variables):
        dvarlist.extend(['300mb_winds','500mb_hgt','500mb_temp',
                         '500mb_vort','700mb_rh','850mb_vort',
                         '850mb_temp','divergence','froude no',
                         'shear0_1km','shear0_3km','shear0_6km',
                         'wind-3d','pot_vort'])

    if ('RAINNC' in variables and 'RAINSH' in variables
        and 'RAINC' in variables):
        dvarlist.extend(['acc_pcp','rainrate'])

    if ('QCLOUD' in variables and 'QRAIN' in variables
        and 'QSNOW' in variables and 'QICE' in variables
        and 'QGRAUPEL' in variables):
        dvarlist.extend(['clwp','iwp','twp','ice water',
                         'cloud water','total water'])

    if ('QRAIN' in variables and 'QSNOW' in variables
        and 'QGRAUPEL' in variables) or ('REFL_10CM' in variables):
        dvarlist.append('refl')

    if ('SWDNT' in variables and 'SWUPTC' in variables
        and 'SWUPT' in variables):
        dvarlist.append('cloud albedo')
    dvars = np.asarray(dvarlist)
    dvars.flatten()
    return sorted(dvars)

