import numpy as np
from libc.math cimport exp, log
import sys


################  Begin function of hypsometric  ###########################
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
    print("In the Cython program")
    for i in range(nx):
        for j in range(ny):
            #for k in range(zz):
                 for kk in range(nz):
                     diff = pressure[kk,j,i] - ptarget
                     if (diff <= 0):
                          if (kk == 0):
                              z2[j,i] = -9999.
                          else:
                              ind = kk - 1
                              #Interpolate to level
                              t_avg = (t[ind+1,j,i] + t[ind,j,i])/2.
                              z2[j,i] = (z1[ind,j,i] + ((r*t_avg)/g)* 
                                          log(pressure[ind,j,i]/ptarget))
                          
                              break
                          
                     if (kk == nz-1):
                         z2[j,i] = -9999.

    return z2
   

###############################################################################
####################  End function hypsometric()  ##########################
###############################################################################


################  Begin function of loglinear_interpolate  ####################
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
    print("In the Cython program")
    for i in range(nx):
        for j in range(ny):
            #for k in range(zz):
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
                    if (kk == nz-1):
                        var_new[j,i] = -9999.

    return var_new

###############################################################################
###################  End function loglinear_interpolate()  ####################
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
    cdef double[:,:,:] theta = np.zeros((nz,ny,nx))
    cdef double[:,:,:] tempv = np.zeros((nz,ny,nx))
    cdef double[:,:,:] p = np.zeros((nz,ny,nx))
    cdef double[:,:] theta_lcl = np.zeros((ny,nx))
    cdef double[:,:] q_lcl = np.zeros((ny,nx))
    cdef double[:] t_parv = np.zeros(nz)
    cdef double t_pt, vpsat, qsat, t_par
    cdef double dz, dtparcel, dtenv, dum1
    cdef int i, j, k
    cdef double ep = 18.016/29.87  #Mass vapor over mass of dry air
    cdef double rd = 287.05
    cdef double cp = 1004.0
    cdef double g = 9.81

    #Start the program
    print("In the Cython program")
    for i in range(nx):
        for j in range(ny):
            #Fill arrays with first guess
            theta_lcl[j,i] = theta[0,j,i]
            q_lcl[j,i] = qv[0,j,i]            
            for k in range(nz):    
                #Convert to hPa
                p[k,j,i] = p_in[k,j,i] * 0.01
                #Calculate virtual temperature
                tempv[k,j,i] = t[k,j,i] * (1. + 0.608 * qv[k,j,i])
                #Calculate potential temperature
                theta[k,j,i] = t[k,j,i] * (1000./p[k,j,i])**(rd/cp) 
                #Fill arrays with first guess
                if (k == 0):
                    theta_lcl[j,i] = theta[0,j,i]
                    q_lcl[j,i] = qv[0,j,i]                                
                #Calculate saturation mixing ratio along constant potential temperature
                t_pt = theta[0,j,i] * (p[k,j,i]/1000.)**(rd/cp) - 273.15 
                vpsat = 6.112*exp(17.67*t_pt/(243.5+t_pt))
                qsat = ep*vpsat / (p[k,j,i] - vpsat)

                #Calculate the parcel temperature
                if (qsat > qv[0,j,i]):
                    #Parcel is not saturated
                    t_par = theta[0,j,i] * (p[k,j,i]/1000.)**(rd/cp)
                    t_parv[k] = t_par * (1. + 0.608 * qv[0,j,i])
                else:
                    #Parcel is saturated
                    t_par, q_lcl[j,i], theta_lcl[j,i] = TParcel(p[k,j,i],q_lcl[j,i],theta_lcl[j,i])
                    t_parv[k] = t_par * (1. + 0.608 * q_lcl[j,i])
                if (k > 0):
                    dz = Z[k,j,i] - Z[k-1,j,i]
                    dtparcel = (t_parv[k]+t_parv[k-1])/2.
                    dtenv = (tempv[k,j,i]+tempv[k-1,j,i])/2.
                    dum1 = g * dz * (dtparcel-dtenv) / dtenv
                    if (dum1 < 0):
                        dum1 = 0
                    cape[j,i] = cape[j,i] + dum1 
     
    return cape, cin
 
###############################################################################
####################  End function cape_sb()
##########################
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
    dqf = 1./3.
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
