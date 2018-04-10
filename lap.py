# This class handles the data loading and plotting for lapping 
# events between two LEO spacecraft. 

import csv
import matplotlib.pyplot as plt
from matplotlib.dates import num2date
import dateutil.parser
import numpy as np
from datetime import datetime, timedelta
import sys
import os

import spacepy.datamodel
import spacepy.time

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

class Lap():
    def __init__(self, sepPath, fb_id, ac_id, fbDir=None, acDir=None):
        """
        This class handles the data management and plotting of the 
        FIREBIRD-II - AC6 lapping events. This class needs 
        A) access to the FIREBIRD HiRes data 
        B) access to AC6's 10Hz or survey data
        C) csv separation file.
        """
        self.fb_id = fb_id
        self.ac_id = ac_id
        self.fbDir = fbDir

        if self.fbDir is None:
            self.fbDir = ('/home/mike/research/firebird/Datafiles/'
                        'FU_{}/hires/level2/'.format(self.fb_id))

        self._load_sep(sepPath) # Load separation file.
        return

    def plot_lap_event(self, tRange, acDtype='10Hz'):
        """ This method makes the lapping event plot between FB and AC6 """
        self._load_fb_data(tRange)
        self._load_ac_data(tRange, acDtype)
        self._match_lat_lag(tRange)

        return

    def _load_fb_data(self, tRange):
        """ This method loads in the FIREBIRD-II data """
        hrName = 'FU{}_Hires_{}_L2.txt'.format(self.fb_id, tRange[0].date())
        self.hr = spacepy.datamodel.readJSONheadedASCII(
                            os.path.join(self.fbDir, hrName)).copy()
        self.hr["Time"] = spacepy.time.Ticktock(self.hr["Time"]).UTC
        self.fb_time_shift = np.mean(self.hr['Count_Time_Correction'])
        # self.hr['shifted_time'] = np.array([t + timedelta(seconds=self.fb_time_shift)
        #                             for t in self.hr['Time']])
        return

    def _load_ac_data(self, tRange, dType='10Hz'):
        """ This method loads in the AC-6 data """
        self.acData = read_ac_data.read_ac_data_wrapper(self.ac_id,
            tRange[0], dType=dType, plot=False)
        return

    def _match_lat_lag(self, tRange):
        """ 
        This function will take the AC-6 data, and find the 
        time lag when the latitudes of the two spacecraft line up.
        This code caluclates the start lag at the start and end
        of the interval to check if there were substantial 
        differences in the time lags.
        """
        # Find the first order correction
        sepInd = np.where(self.sep['dateTime'] > tRange[0])[0][0]
        self.ac_time_lag = self.sep['d'][sepInd]/7.5

        # Now find the higher order correction.
        fbSind = np.where(self.hr['Time'] > tRange[0])[0][0]
        print(self.hr['Time'][fbSind])
        return

    def _load_sep(self, fPath):
        """
        This method loads in the separation data file and saves it so self.sep
        """
        with open(fPath) as f:
            r = csv.reader(f)
            keys = next(r)
            rData = np.array(list(r))
        self.sep = {}
        self.sep['dateTime'] = np.array([dateutil.parser.parse(t) for t in rData[:, 0]])
        self.sep['d'] = np.array([float(f) for f in rData[:, 1]])
        return

lapDict = {'0319T0104':{'tRange':[datetime(2018, 3, 19, 1, 4, 30), 
                                datetime(2018, 3, 19, 1, 6, 0)],
                        'latYlim':(55, 65), 'lonYlim':(-43, -30)},
           '0326T0213':{'tRange':[datetime(2018, 3, 26, 2, 13), 
                                datetime(2018, 3, 26, 2, 18)],
                                'latYlim':None, 'lonYlim':None},
           '0326T0348':{'tRange':[datetime(2018, 3, 26, 3, 50, 25), 
                                datetime(2018, 3, 26, 3, 52)],
                                'latYlim':None, 'lonYlim':None},
           '0326T0051':{'tRange':[datetime(2018, 3, 26, 0, 52, 0), 
                                datetime(2018, 3, 26, 0, 53, 30)],
                                'latYlim':None, 'lonYlim':None}                    
          }        

if __name__ == '__main__':
    fb_id = 3
    ac_id = 'A'
    dPath = '2018-02-26_2018-03-29_FU{}_AC6{}_dist.csv'.format(fb_id, ac_id)
    lapTime = '0326T0213'
    le = LappingEvents(dPath, fb_id, ac_id)
    le.make_plot(lapDict[lapTime]['tRange'], reverse=False, 
                    latYlim=lapDict[lapTime]['latYlim'], 
                    lonYlim=lapDict[lapTime]['lonYlim'])