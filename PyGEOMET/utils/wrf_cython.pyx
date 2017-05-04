import numpy as np
from libc.math cimport exp, log

################  Begin function hypsometric  ###########################
#This subroutine uses the hypsometric 

def hypsometric(float[:,:,:] z1, float[:,:,:] pressure, 
                float ptarget, float[:,:,:] t):

    #Define variables
    cdef int nx = z1.shape[2]
    cdef int ny = z1.shape[1]
    cdef int nz = z1.shape[0]
    #cdef int zz = ptarget.shape[0]
    cdef double[:,:] z2 = np.zeros((ny,nx))
    cdef int i, j, k, kk, ind
    cdef double diff
    cdef double t_avg
    cdef float r = 287.05
    cdef float g = 9.81

    #Start the program
    for i in range(nx):
        for j in range(ny):
            #for k in range(zz):
            #Check if ptarget is below surface before looping
            if pressure[0,j,i] < ptarget:
                z2[j,i] = np.nan
            else:
                 for kk in range(nz):
                     diff = pressure[kk,j,i] - ptarget
                     if (diff <= 0):
                          if (kk == 0):
                              z2[j,i] = np.nan
                          else:
                              ind = kk - 1
                              #Interpolate to level
                              t_avg = (t[ind+1,j,i] + t[ind,j,i])/2.
                              z2[j,i] = (z1[ind,j,i] + ((r*t_avg)/g)* 
                                          log(pressure[ind,j,i]/ptarget))
                          
                              break
                     #Condition if ptarget is less than pressure at model top                      
                     if (kk == nz-1):
                         z2[j,i] = np.nan

    return z2
   

###############################################################################
####################  End function hypsometric()  ##########################
###############################################################################


################  Begin function loglinear_interpolate  ####################
#This function interpolates a variable to a target pressure level

def loglinear_interpolate(float[:,:,:] var, float[:,:,:] pressure, 
                          float ptarget):

    #Define variables
    cdef int nx = var.shape[2]
    cdef int ny = var.shape[1]
    cdef int nz = var.shape[0]
    #cdef int zz = 1 #ptarget.shape[0]
    cdef double[:,:] var_new = np.zeros((ny,nx))
    cdef int i, j, k, kk, ind
    cdef double diff



    #Start the program
    for i in range(nx):
        for j in range(ny):
            #for k in range(zz):
            #Check if ptarget is below surface before looping
            if pressure[0,j,i] < ptarget:
                var_new[j,i] = np.nan
            else:     
                for kk in range(nz):
                    diff = pressure[kk,j,i] - ptarget
                    if (diff <= 0):
                        if (kk == 0):
                            ind = kk
               	            #Interpolate to level
                            var_new[j,i] = ( var[ind+1,j,i] - ((var[ind+1,j,i]-var[ind,j,i])/ 
                                             (log(pressure[ind+1,j,i]/pressure[ind,j,i])))
                                             *log(pressure[ind+1,j,i]/ptarget) )
                        else:
                            ind = kk - 1
                            #Interpolate to level
                            var_new[j,i] = ( var[ind+1,j,i] - ((var[ind+1,j,i]-var[ind,j,i])/ 
                                             (log(pressure[ind+1,j,i]/pressure[ind,j,i])))
                                             *log(pressure[ind+1,j,i]/ptarget) )
                        break
                    #Condition if ptarget is less than pressure at model top
                    if (kk == nz-1):
                        var_new[j,i] = np.nan

    return var_new

###############################################################################
###################  End function loglinear_interpolate()  ####################
###############################################################################

################  Begin function linear_interpolate  ####################
#This function interpolates a variable to a target height level

def linear_interpolate(float[:,:,:] var, float[:,:,:] z,
                       float ztarget):

    #Define variables
    cdef int nx = var.shape[2]
    cdef int ny = var.shape[1]
    cdef int nz = var.shape[0]
    cdef double[:,:] var_new = np.zeros((ny,nx))
    cdef int i, j, k, ind
    cdef double diff

    #Start the program
    for i in range(nx):
        for j in range(ny):
            #Check if ztarget is below surface before looping
            if z[0,j,i] > ztarget:
                var_new[j,i] = np.nan
            else:
                for k in range(nz):
                    diff = z[k,j,i] - ztarget
                    if (diff >= 0):
                        ind = k - 1
                        #Interpolate to level
                        var_new[j,i] = ( var[ind+1,j,i] - ((var[ind+1,j,i]-var[ind,j,i])/
                                         (z[ind+1,j,i]-z[ind,j,i]))*(z[ind+1,j,i]-ztarget))
                        break
                    #Condition if ztarget is less than z at model top
                    if (k == nz-1):
                        var_new[j,i] = np.nan

    return var_new

###############################################################################
###################  End function linear_interpolate()  ####################
###############################################################################

################  Begin function linear_interpolate1d  ####################
#This function interpolates a variable to a target height level
#Input is a 1D array
def linear_interpolate1D(float[:] var, float[:] z,
                       float ztarget):

    #Define variables
    cdef int nz = var.shape[0]
    cdef double var_new 
    cdef int k, ind
    cdef double diff

    #Start the program
    #Check if ztarget is below surface before looping
    if z[0] > ztarget:
        var_new = np.nan
    else:
        for k in range(nz):
            diff = z[k] - ztarget
            if (diff >= 0):
                ind = k - 1
                #Interpolate to level
                var_new = ( var[ind+1] - ((var[ind+1]-var[ind])/
                            (z[ind+1]-z[ind]))*(z[ind+1]-ztarget))
                break
            #Condition if ztarget is less than z at model top
            if (k == nz-1):
                var_new = np.nan

    return var_new

###############################################################################
###################  End function linear_interpolate()  ####################
###############################################################################


################  Begin function mean_layer  ####################
#This function calculates the mean of a layer

def mean_layer(float[:,:,:] var, float[:,:,:] z, 
               float[:,:] ref1, float[:,:] ref2):

    #Define variables
    cdef int nx = var.shape[2]
    cdef int ny = var.shape[1]
    cdef int nz = var.shape[0]
    cdef double[:,:] var_ml = np.zeros((ny,nx))
    cdef double diff, diff2
    cdef int i, j, k, ind1, ind2

    for i in range(nx):
        for j in range(ny):
            #Find the value at reference level 1
            ind1 = 1
            diff = ref1[j,i] - z[0,j,i]
            while (diff > 0):
                diff = ref1[j,i] - z[ind1,j,i]
                #Get out of the loop
                if (ind1 == nz):
                    #To keep minus 2 working below
                    ind1 = ind1 + 1
                    break
                ind1 = ind1 + 1
            #Find layer below condition
            #Minus 2 because we add 1 before conditional check
            ind1 = ind1 - 2
            if (ind1 < 0):
                ind1 = 0
            #Find the value at reference level 2
            ind2 = 1
            diff2 = ref2[j,i] - z[0,j,i]
            while (diff2 > 0):
                diff2 = ref2[j,i] - z[ind2,j,i]
                #Get out of the loop
                if (ind2 == nz):
                    #To keep minus 2 working below
                    ind2 = ind2 + 1
                    break
                ind2 = ind2 + 1
            #Find layer below condition
            #Minus 2 because we add 1 before conditional check
            ind2 = ind2 - 2
            if (ind2 < 0):
                ind2 = 0             
            
            #Calculate the mean layer value
            var_ml[j,i] = ( ((var[ind1+1,j,i] + var[ind1,j,i])/2.) +
                            ((var[ind2+1,j,i] + var[ind2,j,i])/2.))/2.              
            

    return var_ml

###############################################################################
###################  End function ean_layer()  ####################
###############################################################################


################  Begin function cap_sb  ###########################
#This subroutine calculates surface based CAPE

def cape_sb(float[:,:,:] t, float[:,:,:] qv, float[:,:,:] p_in, 
            float[:,:,:] Z):

    #Define variables
    cdef int nx = t.shape[2]
    cdef int ny = t.shape[1]
    cdef int nz = t.shape[0]
    cdef double[:,:] cape = np.zeros((ny,nx))
    cdef double[:,:] cin = np.zeros((ny,nx))
    #cdef double[:] theta = np.zeros(nz)
    cdef double[:] tempv = np.zeros(nz)
    cdef double[:] p = np.zeros(nz)
    #cdef double[:,:] theta_lcl = np.zeros((ny,nx))
    #cdef double[:,:] q_lcl = np.zeros((ny,nx))
    cdef double[:] t_parv = np.zeros(nz)
    cdef double t_pt, vpsat, qsat, t_par
    cdef double dz, dtparcel, dtenv, dum1, dum2
    cdef int i, j, k
    cdef double ep = 18.016/29.87  #Mass vapor over mass of dry air
    cdef double rd = 287.05
    cdef double cp = 1004.0
    cdef double g = 9.81
    cdef double theta_lcl, q_lcl, theta

    #Start the program
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):    
                #Convert to hPa
                p[k] = p_in[k,j,i] * 0.01
                #Calculate virtual temperature
                tempv[k] = t[k,j,i] * (1. + 0.608 * qv[k,j,i])
                #Fill arrays with first guess
                if (k == 0):
                    theta = t[k,j,i] * (1000./p[k])**(rd/cp)
                    theta_lcl = theta
                    q_lcl = qv[0,j,i] 
                #Calculate saturation mixing ratio along constant potential temperature
                t_pt = theta * (p[k]/1000.)**(rd/cp) - 273.15 
                vpsat = 6.112*exp(17.67*t_pt/(243.5+t_pt))
                qsat = ep*vpsat / (p[k] - vpsat)

                #Calculate the parcel temperature
                if (qsat > qv[0,j,i]):
                    #Parcel is not saturated
                    t_par = theta * (p[k]/1000.)**(rd/cp)
                    t_parv[k] = t_par * (1. + 0.608 * qv[0,j,i])
                else:
                    #Parcel is saturated
                    t_par, q_lcl, theta_lcl = TParcel(p[k],q_lcl,theta_lcl)
                    t_parv[k] = t_par * (1. + 0.608 * q_lcl)
                if (k > 0):
                    dz = Z[k,j,i] - Z[k-1,j,i]
                    dtparcel = (t_parv[k]+t_parv[k-1])/2.
                    dtenv = (tempv[k]+tempv[k-1])/2.
                    dum1 = g * dz * (dtparcel-dtenv) / dtenv
                    if (dum1 < 0):
                        dum1 = 0
                    cape[j,i] = cape[j,i] + dum1 
                    if (p[k] > p[0] - 300.):
                        dum2 = g * dz * (dtparcel-dtenv) / dtenv
                        if (dum2 > 0):
                            dum2 = 0
                        cin[j,i] = cin[j,i] + dum2
    
    return cape, cin
 
###############################################################################
####################  End function cape_sb()  ##########################
###############################################################################

################  Begin function cape_ml  ###########################
#This subroutine calculates mixed-layer CAPE/CIN
#Mixed layer is defined as the average of the surface to 100mb above
#  surface layer

def cape_ml(float[:,:,:] t, float[:,:,:] qv, float[:,:,:] p_in,
            float[:,:,:] Z):

    #Define variables
    cdef int nx = t.shape[2]
    cdef int ny = t.shape[1]
    cdef int nz = t.shape[0]
    cdef double[:,:] cape = np.zeros((ny,nx))
    cdef double[:,:] cin = np.zeros((ny,nx))
    cdef double[:] tempv = np.zeros(nz)
    cdef double[:] p = np.zeros(nz)
    cdef double[:] t_parv = np.zeros(nz)
    cdef double t_pt, vpsat, qsat, t_par
    cdef double dz, dtparcel, dtenv, dum1, dum2
    cdef int i, j, k, count
    cdef double ep = 18.016/29.87  #Mass vapor over mass of dry air
    cdef double rd = 287.05
    cdef double cp = 1004.0
    cdef double g = 9.81
    cdef double total_temp, total_qv, total_the, t_ml, q_ml, theta_ml
    cdef double q_lcl, theta_lcl

    #Start the program
    for i in range(nx):
        for j in range(ny):
            #Calculate mixed layer values
            minP = (p_in[0,j,i]*0.01) - 100.
            #Determine where minP occurs in the column
            count = 0
            total_temp = 0
            total_qv = 0 
            total_the = 0
            while (count <= nz):
                ptmp = (p_in[count,j,i]*0.01)
                if (ptmp < minP):
                    break
                else:
                    #Sum up temperature, vapor, theta
                    total_temp = total_temp + t[count,j,i]
                    total_qv = total_qv + qv[count,j,i]
                    total_the = total_the + (t[count,j,i] * (1000./ptmp)**(rd/cp))
                    count = count + 1
                #Fail safe to surface based parcel
                if (count == nz):
                    total_temp = t[0,j,i]
                    total_qv = qv[0,j,i]
                    total_the = t[0,j,i] * (1000./(p_in[0,j,i]*0.01))**(rd/cp)
                    count = 1
                    break
            #Take the average of the layer
            t_ml = total_temp / count
            q_ml = total_qv / count
            theta_ml = total_the / count

            #Set new values for theta and q at the lcl height
            #Will change along k once parcel becomes saturated
            theta_lcl = theta_ml
            q_lcl = q_ml

            for k in range(nz):
                #Convert to hPa
                p[k] = p_in[k,j,i] * 0.01
                #Calculate virtual temperature
                tempv[k] = t[k,j,i] * (1. + 0.608 * qv[k,j,i])
                #Calculate saturation mixing ratio along constant potential temperature
                t_pt = theta_ml * (p[k]/1000.)**(rd/cp) - 273.15
                vpsat = 6.112*exp(17.67*t_pt/(243.5+t_pt))
                qsat = ep*vpsat / (p[k] - vpsat)

                #Calculate the parcel temperature
                if (qsat > q_ml):
                    #Parcel is not saturated
                    t_par = theta_lcl * (p[k]/1000.)**(rd/cp)
                    t_parv[k] = t_par * (1. + 0.608 * q_lcl)
                else:
                    #Parcel is saturated
                    t_par, q_lcl, theta_lcl = TParcel(p[k],q_lcl,theta_lcl)
                    t_parv[k] = t_par * (1. + 0.608 * q_lcl)
                if (k > 0):
                    dz = Z[k,j,i] - Z[k-1,j,i]
                    dtparcel = (t_parv[k]+t_parv[k-1])/2.
                    dtenv = (tempv[k]+tempv[k-1])/2.
                    dum1 = g * dz * (dtparcel-dtenv) / dtenv
                    if (dum1 < 0):
                        dum1 = 0
                    cape[j,i] = cape[j,i] + dum1
                    if (p[k] > p[0] - 300.):
                        dum2 = g * dz * (dtparcel-dtenv) / dtenv
                        if (dum2 > 0):
                            dum2 = 0
                        cin[j,i] = cin[j,i] + dum2

    return cape, cin

###############################################################################
####################  End function cape_ml()  ##########################
###############################################################################

################  Begin function cape_mu  ###########################
#This subroutine calculates most unstable CAPE/CIN
#Most unstable is defined as the max theta e layer (sfc to sfc-300mb)

def cape_mu(float[:,:,:] t, float[:,:,:] qv, float[:,:,:] p_in,
            float[:,:,:] Z):

    #Define variables
    cdef int nx = t.shape[2]
    cdef int ny = t.shape[1]
    cdef int nz = t.shape[0]
    cdef double[:,:] cape = np.zeros((ny,nx))
    cdef double[:,:] cin = np.zeros((ny,nx))
    cdef double[:] tempv = np.zeros(nz)
    cdef double[:] p = np.zeros(nz)
    cdef double[:] t_parv = np.zeros(nz)
    cdef double t_pt, vpsat, qsat, t_par
    cdef double dz, dtparcel, dtenv, dum1, dum2
    cdef int i, j, k, count
    cdef double ep = 18.016/29.87  #Mass vapor over mass of dry air
    cdef double rd = 287.05
    cdef double cp = 1004.0
    cdef double g = 9.81
    cdef double thetae, t_mu, q_mu, theta_mu, theta, thetae_max
    cdef double q_lcl, theta_lcl, t_lcl, vp, td

    #Start the program
    for i in range(nx):
        for j in range(ny):
            #Calculate mixed layer values
            minP = (p_in[0,j,i]*0.01) - 300.
            #Determine the max thetae
            count = 0
            thetae_max = 0.
            while (count <= nz):
                ptmp = (p_in[count,j,i]*0.01)
                if (ptmp < minP):
                    break
                else:
                    #Calculate vapor pressure
                    vp = (qv[count,j,i]*(p_in[count,j,i]*0.01)) / (qv[count,j,i] + ep)
                    #td = (243.5 * log(vp/6.112)) / (17.67 - log(vp/6.112)) + 273.15
                    #Calulate lcl temperature
                    t_lcl = (2840./(3.5*log(t[count,j,i]) - log(vp)-4.805)) + 55.
                    #Calculate theta
                    theta = t[count,j,i] * (1000./(p_in[count,j,i]*0.01))**(rd/cp)
                    #Calculate thetae
                    thetae = (theta*exp((3.376/t_lcl - 0.00254)*1000.*qv[count,j,i]*
                              (1.+0.81*qv[count,j,i])))
                    #Check if greater than max
                    if (thetae > thetae_max):
                        thetae_max = thetae
                        t_mu = t[count,j,i]
                        q_mu = qv[count,j,i]
                        theta_mu =  t[count,j,i] * (1000./(p_in[count,j,i]*0.01))**(rd/cp)
                    count = count + 1
                #Fail safe to surface based parcel
                if (count == nz):
                    t_mu = t[0,j,i]
                    q_mu = qv[0,j,i]
                    theta_mu = t[0,j,i] * (1000./(p_in[0,j,i]*0.01))**(rd/cp)
                    count = 1
                    break

            #Set new values for theta and q at the lcl height
            #Will change along k once parcel becomes saturated
            theta_lcl = theta_mu
            q_lcl = q_mu

            for k in range(nz):
                #Convert to hPa
                p[k] = p_in[k,j,i] * 0.01
                #Calculate virtual temperature
                tempv[k] = t[k,j,i] * (1. + 0.608 * qv[k,j,i])
                #Calculate saturation mixing ratio along constant potential temperature
                t_pt = theta_mu * (p[k]/1000.)**(rd/cp) - 273.15
                vpsat = 6.112*exp(17.67*t_pt/(243.5+t_pt))
                qsat = ep*vpsat / (p[k] - vpsat)

                #Calculate the parcel temperature
                if (qsat > q_mu):
                    #Parcel is not saturated
                    t_par = theta_lcl * (p[k]/1000.)**(rd/cp)
                    t_parv[k] = t_par * (1. + 0.608 * q_lcl)
                else:
                    #Parcel is saturated
                    t_par, q_lcl, theta_lcl = TParcel(p[k],q_lcl,theta_lcl)
                    t_parv[k] = t_par * (1. + 0.608 * q_lcl)
                if (k > 0):
                    dz = Z[k,j,i] - Z[k-1,j,i]
                    dtparcel = (t_parv[k]+t_parv[k-1])/2.
                    dtenv = (tempv[k]+tempv[k-1])/2.
                    dum1 = g * dz * (dtparcel-dtenv) / dtenv
                    if (dum1 < 0):
                        dum1 = 0
                    cape[j,i] = cape[j,i] + dum1
                    if (p[k] > p[0] - 300.):
                        dum2 = g * dz * (dtparcel-dtenv) / dtenv
                        if (dum2 > 0):
                            dum2 = 0
                        cin[j,i] = cin[j,i] + dum2

    return cape, cin

###############################################################################
####################  End function cape_mu()  ##########################
###############################################################################

################# Start Subroutine TParcel######################################
#This subroutine calculates the saturated parcel temperature

cdef TParcel(float p, float qv, float theta):

    cdef float t_par
    cdef float tk, vp, qvsat, dq, dqf 
    cdef float tk_old, dq_old, theta_old, qv_old, qvsat_old
    #Constants
    cdef float rd = 287.05
    cdef float cp = 1004.7
    cdef float kpa = rd/cp
    cdef float ep = 18.016/29.87
    cdef float L = 2.5e6
    cdef float p0 = 1000.0
    cdef float dqthres = 0.00001

    #Calculate temperature
    tk = theta * (p/p0)**kpa
    #Calculate saturation vapor pressure
    vp = 6.112*exp(17.67*(tk-273.15)/(243.5+(tk-273.15)))
    #Calculate saturation mixing ratio
    qvsat = ep*vp / (p-vp)
    #Calculate the fraction difference between the input mixing ratio and saturation
    dqf = 1./2.5
    dq = dqf * (qvsat - qv)

    #Loop until the difference between saturation and actual are negligible
    while (abs(dq) > dqthres):
        #Set old values for fail safe
        dq_old = dq
        theta_old = theta
        qv_old = qv
        tk_old = tk
        qvsat_old = qvsat
        #Calculate new values
        dq = dqf * (qvsat - qv)
        theta = theta - (dq*L*theta) / (cp*tk)
        qv = qv + dq
        tk = theta * (p/p0)**kpa
        vp = 6.112*exp(17.67*(tk-273.15)/(243.5+(tk-273.15)))
        qvsat = ep*vp / (p-vp)
        #Check to make sure qvsat is less that the guess qv
        #If so, reduce the iterator and recalculate
        if (qvsat > qv):
            dqf = 0.5 * dqf
            theta = theta_old
            qv = qv_old
            tk = tk_old
            qvsat = qvsat_old
      
    t_par = (theta*(p/p0)**kpa)

    return t_par, qv, theta
