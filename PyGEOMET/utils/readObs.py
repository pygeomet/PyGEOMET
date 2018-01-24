import numpy as np
import csv
import time
import datetime as dt

#This function is designed to read in processed ds472.0 surface data
#Inputs: filename - observation filename
#        begdate  - start date for the loaded WRF dataset
#        enddate  - end date for the loaded WRF dataset
#Written by: Andrew White (andrew.white@nsstc.uah.edu)
#Date: Jan 22, 2018
def readObs(filename,begdate,enddate):
    print(filename)
    #Transform dates to date time objects
    #Necessary because the input is a string
    #Start
    syr = int(begdate[7:11])
    smo = int(time.strptime(begdate[3:6], '%b').tm_mon)
    sdy = int(begdate[0:3])
    shour = int(begdate[13:15])
    bdate = dt.datetime(syr,smo,sdy,shour)
    #End
    eyr = int(enddate[7:11])
    emo = int(time.strptime(enddate[3:6], '%b').tm_mon)
    edy = int(enddate[0:3])
    ehour = int(enddate[13:15])
    edate = dt.datetime(eyr,emo,edy,ehour)

    #Setup list to hold variables
    date = []
    lat = []
    lon = []
    t = []
    q = []
    ws = []
    wd = []
    append_date = date.append
    append_lat = lat.append
    append_lon = lon.append
    append_t = t.append
    append_q = q.append
    append_ws = ws.append
    append_wd = wd.append
    #Read file
    with open(filename,'r') as f:
        reader = csv.reader(f)
        for i in range(6):
            next(reader,None)
        for row in reader:
            yy = int(row[0][0:4])
            mm = int(row[0][5:7])
            dd = int(row[0][8:10])
            hh = int(row[0][12:16])/100
            d_tmp = dt.datetime(yy,mm,dd,hh)
            if (d_tmp > edate):
                break
            if (d_tmp >= bdate and d_tmp <= edate):
                append_date(d_tmp)
                append_lat(row[0][28:35])
                append_lon(row[0][37:45])
                append_t(row[0][54:61])
                append_q(row[0][69:76])
                append_ws(row[0][84:91])
                append_wd(row[0][99:106])


    #Create numpy arrays and convert strings to floats
    date = np.array(date)
    lat = np.array(lat,dtype=np.float)
    lon = np.array(lon,dtype=np.float)
    t = np.array(t,dtype=np.float)
    q = np.array(q,dtype=np.float)
    ws = np.array(ws,dtype=np.float)
    wd = np.array(wd,dtype=np.float)

    #Fill missing data with NaNs (-999 is missing)
    t[t < -900] = np.nan
    q[q < -900] = np.nan
    ws[ws < -900] = np.nan
    wd[wd < -900] = np.nan

    return date, lat, lon, t, q, ws, wd 
