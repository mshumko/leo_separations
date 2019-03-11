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
sys.path.append('/home/mike/research/mission_tools/ac6/')
import read_ac_data

Re=6371 # km

class CalcDist():
    def __init__(self, scA, scB, startDate, endDate, aEphem, bEphem=False):
        """
        This class loads in two ephemeris files that were generated by 
        SGP4, or the daily AC6 coords files, and calculates the total
        distance, in-track, and cross-track separation between the two
        spacecraft.
        """
        self.scA = scA
        self.scB = scB
        self.startDate = startDate
        self.endDate = endDate

        # Load ephem data.
        self.aEphem = self._load_mike_ephem(aEphem)
        # Load AC6 data
        if not bEphem:
            self.bEphem = self._load_ac_ephem()
        else:
            self.bEphem = self._load_mike_ephem(bEphem)
        self._find_common_times() # Filter data by the same time stamps.
        return
        
    def calc_dist(self):
        """
        This is a wrapper function that calculates the total separation using
        greatCircleDist(). This function calculates the psudo in-track lag
        from the difference in latitude. Positive in-track separation implies 
        that you add the lag (or separation) to spacecraft B. 
        """
        # Format inputs for haversine method.
        X1 = np.array([self.aEphem['lat'], self.aEphem['lon'], self.aEphem['alt']]).T
        X2 = np.array([self.bEphem['lat'], self.bEphem['lon'], self.bEphem['alt']]).T
        
        direction = np.convolve([0.5, -0.5], self.bEphem['lat'], mode='same')
        
        self.dTot = self._haversine(X1, X2) # Get total distance
        A = Re+(X1[:, 2]+X2[:, 2])/2 # Mean altitude
        # Find a rough fraction of total distance that is in-track.
        
        self.dInTrack = np.pi/180*A*(X1[:, 0] - X2[:, 0])*np.sign(direction)
        # Use Pathagorean theorem to calculate the cross-track separation.
        self.dCrossTrack = np.sqrt(self.dTot**2 - self.dInTrack**2)
        return self.dTot, self.dInTrack, self.dCrossTrack

    def save_file(self, saveName):
        """
        This method saves the separation data into a csv file.
        """

        with open(saveName, 'w', newline='') as f:
            w = csv.writer(f)
            # Write header
            w.writerow(['dateTime', 'dist_in_track [km]', 
                        'dist_cross_track [km]', 
                        'L_{}'.format(self.scA), 
                        'L_{}'.format(self.scB), 
                        'MLT_{}'.format(self.scA), 
                        'MLT_{}'.format(self.scB)])
            # Save data.
            zz = zip(self.aEphem['dateTime'], self.dInTrack, self.dCrossTrack, 
                    self.aEphem['L'], self.bEphem['L'], self.aEphem['MLT'], 
                    self.bEphem['MLT'])
            for z in zz:
                w.writerow([*z])
            return

    def plot_dist(self):
        """

        """
        fig, ax = plt.subplots(3, figsize=(10, 8), sharex=True)
        ax_t = ax[1].twinx()

        ax[0].plot(self.aEphem['dateTime'], self.dTot)
        ax[1].plot(self.aEphem['dateTime'], self.dInTrack)
        ax_t.plot(self.aEphem['dateTime'], self.dInTrack/7.5) # Assuming a 7.5 km/s orbital velocity
        ax[2].plot(self.aEphem['dateTime'], self.dCrossTrack)
        
        ax[0].set_title('{}-{} | {}-{} separation'.format(
            self.startDate.date(), self.endDate.date(), 
            self.scA, self.scB))
        ax[0].set_ylabel('total separation [km]')
        ax[1].set_ylabel('in-track separation [km]')
        ax[2].set(ylabel='cross-track separation [km]')
        ax_t.set_ylabel('In-track lag [s]')
        ax[-1].set_xlabel('UTC')
        
        return

    def _load_mike_ephem(self, fPath):
        """
        This method loads in the ephemeris (magnetic ephemeris) that was 
        generated by Mike's SGP4 algorithm implementation. 
        """
        ephem = {}
        keys = ['dateTime', 'lat', 'lon', 'alt', 'L', 'MLT']

        with open(fPath) as f:
            r = csv.reader(f, quotechar='"')
            next(r)
            #keys = next(r)
            # Strip leading whitespce in keys
            #keys = [s.lstrip() for s in keys]        
            rawData = list(r)
        for (i, key) in enumerate(keys):
            ephem[key] = np.array(rawData)[:, i]
            if key != 'dateTime':
                ephem[key] = ephem[key].astype(float)
            else:
                ephem['dateTime'] = np.array([dateutil.parser.parse(t) 
                                    for t in ephem['dateTime']])
        return ephem
        
    def _load_ac_ephem(self):
        """
        This function will load in the coords data type from the AC6 directory
        and append them all to each other.
        """
        ephem={}
        ephem['dateTime'] = np.array([])
        ephem['lat'] = np.array([])
        ephem['lon'] = np.array([])
        ephem['alt'] = np.array([])
        ephem['L'] = np.array([])
        ephem['MLT'] = np.array([])

        days = [self.startDate + timedelta(t) for t in 
                range((self.endDate - self.startDate).days+1)]
        for d in days:
            # Load AC-6 coordinates
            try:
                rawAc = read_ac_data.read_ac_data_wrapper('A', d, dType='coords') 
            except AssertionError as err:
                if 'None or > 1 AC6 files found' in str(err): 
                    continue
                else:
                    raise
            ephem['dateTime'] = np.append(ephem['dateTime'], rawAc['dateTime'])
            ephem['lat'] = np.append(ephem['lat'], rawAc['lat'])
            ephem['lon'] = np.append(ephem['lon'], rawAc['lon'])
            ephem['alt'] = np.append(ephem['alt'], rawAc['alt'])
            ephem['L'] = np.append(ephem['L'], rawAc['Lm_OPQ'])
            ephem['MLT'] = np.append(ephem['MLT'], rawAc['MLT_OPQ'])
        return ephem

    def _find_common_times(self):
        """
        This method filters the two ephemeris files to the same time
        stamps.
        """
        # Find common times
        fbT = date2num(self.aEphem['dateTime'])
        acT = date2num(self.bEphem['dateTime'])
        fbInd = np.where(np.in1d(fbT, acT))[0]
        acInd = np.where(np.in1d(acT, fbT))[0]

        # Filter data
        for key in self.aEphem:
            self.aEphem[key] = self.aEphem[key][fbInd]
        for key in self.bEphem:
            self.bEphem[key] = self.bEphem[key][acInd]
        return

        
    def _haversine(self, X1, X2):
        """
        Implementation of the haversine foruma to calculate total distance
        at an average altitude. X1 and X2 must be N*3 array of 
        lat, lon, alt.
        """
        X1 = np.asarray(X1)
        X2 = np.asarray(X2)
        R = (Re+(X1[:, 2]+X2[:, 2])/2)
        s = 2*np.arcsin( np.sqrt( np.sin(np.deg2rad(X1[:, 0]-X2[:, 0])/2)**2 + \
                        np.cos(np.deg2rad(X1[:, 0]))*np.cos(np.deg2rad(X2[:, 0]))*\
                        np.sin(np.deg2rad(X1[:, 1]-X2[:, 1])/2)**2 ))
        return R*s
        
if __name__ == '__main__':
    SC_B = 'REACH'
    for SC_A in ['FU3', 'FU4']:
        DATE_RANGE = [datetime(2019, 1, 21), datetime(2019, 2, 23)]
        pathA = ('./data/magephem/{}_{}_{}_magephem.csv'.format(
                SC_A, DATE_RANGE[0].date(), (DATE_RANGE[1]).date()))

        pathB = ('./data/magephem/reach.20190125.vid-169.magephem.txt')

        saveDir = ('./data/dist/20190125_{}_{}_dist.csv'.format(
                    SC_A, SC_B))
        c = CalcDist(SC_A, SC_B, *DATE_RANGE, pathA, pathB)
        c.calc_dist()
        c.save_file(saveDir)
        c.plot_dist()
        plt.show()
