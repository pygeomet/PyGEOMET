import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import PyGEOMET.utils

WBAN = []
lat = []
lon = []
fnames = ['NCDCID', 'ICAO', 'WBAN', 'NAME', 'ST', 'COUNTY', 
          'LAT', 'LON', 'ELEV (ft)"', 'TOWER HEIGHT (m)']

def get_sites(): 
    
    path = os.path.dirname(PyGEOMET.utils.__file__)
    with open(os.path.join(path,'radar_sites.csv')) as csvfile:
        reader = csv.DictReader(csvfile,fieldnames=fnames)
        next(reader)
        next(reader)
        for row in reader:
            WBAN.append(row['ICAO'])
            lat.append(np.float(row['LAT']))
            lon.append(np.float(row['LON']))

    return sorted(WBAN), lon, lat

