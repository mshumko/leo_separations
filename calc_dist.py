# This script calculates the cross-spacecraft separation between AC6 and FU3
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as md
import dateutil.parser
from matplotlib.dates import date2num
import sys
import csv

# My libraries
sys.path.append('/home/mike/research/mission-tools/ac6/')
import read_ac_data

Re=6371 # km
FB_ID = '3'
AC_ID = 'A'
DATE_RANGE = [datetime(2018, 4, 11), datetime(2018, 6, 11)]


def greatCircleDist(X1, X2):
    """
    X1 and X2 must be N*3 array of lat, lon, alt. phi = lat, lambda = lon
    """
    X1 = np.asarray(X1)
    X2 = np.asarray(X2)
    R = (Re+(X1[:, 2]+X2[:, 2])/2)
    s = 2*np.arcsin( np.sqrt( np.sin(np.deg2rad(X1[:, 0]-X2[:, 0])/2)**2 + np.cos(np.deg2rad(X1[:, 0]))*np.cos(np.deg2rad(X2[:, 0]))*np.sin(np.deg2rad(X1[:, 1]-X2[:, 1])/2)**2 ))
    return R*s

def km2lag(d):
    """
    Calculates the in-track lag assuming no cross track separation.
    """
    return d/7.5/60

def convert_ax_d_to_lag(ax):
    """
    Update second axis according with first axis.
    """
    y1, y2 = ax.get_ylim()
    ax_t.set_ylim(km2lag(y1), km2lag(y2))
    ax_t.figure.canvas.draw()
    return

def load_ephem(fPath):
    """
    This function loads in a csv formatted ephemeris file.
    """
    ephem = {}
    with open(fPath) as f:
        r = csv.reader(f, quotechar='"')
        keys = next(r)
        # Strip leading whitespce in keys
        keys = [s.lstrip() for s in keys]        
        rawData = list(r)
    for (i, key) in enumerate(keys):
        ephem[key] = np.array(rawData)[:, i]
        if key != 'Time (ISO)':
            ephem[key] = ephem[key].astype(float)
    ephem['dateTime'] = np.array([dateutil.parser.parse(t) for t
                                 in ephem['Time (ISO)']])
    return ephem

def load_daily_ac6_ephem():
    """
    This function will load in the coords data type from the AC6 directory
    and append them all to each other.
    """
    data={}
    data['dateTime'] = np.array([])
    data['lat'] = np.array([])
    data['lon'] = np.array([])
    data['alt'] = np.array([])

    days = [DATE_RANGE[0] + timedelta(t) for t in 
            range((DATE_RANGE[1] - DATE_RANGE[0]).days+1)]
    for d in days:
        # Load AC-6 coordinates
        try:
            rawAc = read_ac_data.read_ac_data_wrapper('A', d, dType='coords') 
        except AssertionError as err:
            if 'None or > 1 AC6 files found' in str(err): 
                continue
            else:
                raise
        data['dateTime'] = np.append(data['dateTime'], rawAc['dateTime'])
        data['lat'] = np.append(data['lat'], rawAc['lat'])
        data['lon'] = np.append(data['lon'], rawAc['lon'])
        data['alt'] = np.append(data['alt'], rawAc['alt'])
    return data

def save_dist(times, dists, fPath):
    """
    given a time array t and distance array d, this function saves that data
    to a csv file given by fPath.
    """

    with open(fPath, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['dateTime', 'distance [km]'])

        for (t, d) in zip(times, dists):
            w.writerow([t, d])


### Load FIREBIRD coordinates ###
fbPath = ('/home/mike/research/mission-tools/orbit/data/'
          'FU{}_{}_{}_LLA_ephemeris_pre.csv'.format(FB_ID, DATE_RANGE[0].date(), 
        (DATE_RANGE[1]).date()))
acPath = ('/home/mike/research/mission-tools/orbit/data/'
          'AEROCUBE_6{}_{}_{}_LLA_ephemeris_pre.csv'.format(AC_ID, DATE_RANGE[0].date(), 
        (DATE_RANGE[1]).date()))

fb = load_ephem(fbPath)
ac = load_ephem(acPath)
#ac = load_daily_ac6_ephem()

### Now find the same indicies across both data sets ###
fbT = date2num(fb['dateTime'])
acT = date2num(ac['dateTime'])

fbInd = np.where(np.in1d(fbT, acT))[0]
acInd = np.where(np.in1d(acT, fbT))[0]

### Calculate separation ###
#X1 = np.array([ac['lat'][acInd], ac['lon'][acInd], ac['alt'][acInd]]).T
X1 = np.array([ac['Lat (deg)'][acInd], ac['Lon (deg)'][acInd], ac['Alt (km)'][acInd]]).T
X2 = np.array([fb['Lat (deg)'][fbInd], fb['Lon (deg)'][fbInd], fb['Alt (km)'][fbInd]]).T
d = greatCircleDist(X1, X2)

### Save distance data ###
save_dist(fb['dateTime'][fbInd], d, 
            '/home/mike/research/leo-lapping-events/data/'
            '{}_{}_FU{}_AC6{}_dist.csv'.format(
            fb['dateTime'][fbInd[0]].date(), 
            fb['dateTime'][fbInd[-1]].date(), 
            FB_ID, AC_ID))

### Make the distance plot ###
fig, ax = plt.subplots(figsize=(10, 8))
ax_t = ax.twinx()

ax.plot(fb['dateTime'][fbInd], d)
ax_t.plot(fb['dateTime'][fbInd], d/7.5) # Assuming a 7.5 km/s orbital velocity
ax_t.set_ylabel('In-track lag [s] (assuming no cross-track)')
ax.set_title('{}-{} | FU{}-AC6{} total separation'.format(
    DATE_RANGE[0].date(), DATE_RANGE[1].date(), FB_ID, AC_ID))
ax.set_xlabel('UTC')
ax.set_ylabel('Separation [km]')
plt.show()