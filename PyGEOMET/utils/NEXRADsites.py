import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import PyGEOMET.utils

WBAN = []
latlon = []
fnames = ['WBAN #', 'STATION_ID', 'STATION_NAME', 
          'LAT/LON', 'ELEV (ft)"', 'TOWER HEIGHT (m)']

def get_sites(): 
    
    path = os.path.dirname(PyGEOMET.utils.__file__)
    with open(os.path.join(path,'radar_sites.csv')) as csvfile:
        reader = csv.DictReader(csvfile,fieldnames=fnames)
        next(reader)
        for row in reader:
            WBAN.append(row['STATION_ID'])
            latlon.append(row['LAT/LON'])

    def get_latlon(ll):
        lat = []
        lon = []
        for l in ll:
            tmp = l.split('/')[0]
            lat.append(int(tmp[0:2])+int(tmp[2:4])/60. + int(tmp[4:6])/3600.)
            tmp = l.split('/')[1]
            if len(tmp) == 8:
                lon.append(-1*(int(tmp[1:4])+int(tmp[4:6])/60. + int(tmp[6:8])/3600.))
            else:
                lon.append(int(tmp[1:4])+int(tmp[4:6])/60. + int(tmp[6:8])/3600.)

        return lon, lat

    longitude, latitude = get_latlon(latlon)

    return sorted(WBAN), longitude, latitude

