# This class handles the data loading and plotting for lapping 
# events between two LEO spacecraft. 

import csv
import matplotlib.pyplot as plt
import matplotlib.dates
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
    def __init__(self, sepPath, fb_id, ac_id, fbDir=None, acDir=None,
                 startDate=None, endDate=None):
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

        self.startDate = startDate
        self.endDate = endDate

        if self.fbDir is None:
            self.fbDir = ('/home/mike/research/firebird/Datafiles/'
                        'FU_{}/hires/level2/'.format(self.fb_id))

        self._load_sep(sepPath) # Load separation file.
        return

    def plot_lap_event(self, tRange, acDtype='10Hz', lag=None):
        """ This method makes the lapping event plot between FB and AC6 """
        self._load_fb_data(tRange)
        self._load_ac_data(tRange, acDtype)
        
        # Only implement the lag for start of run, and implement 
        # the end time later.
        (self.ac_time_lag, _, start_cross_track, _) = self._get_ac6_lag(tRange) 
        if lag is not None:
            self.ac_time_lag = lag

        fig, ax = plt.subplots(3, figsize=(8, 9))
        self._plot_fb(tRange, ax[0], axPos=ax[2])
        self._plot_ac(tRange, ax[1], axPos=ax[2])

        ### Plot Adjustments ###
        titleStr = ('FU{} - AC6{} Lapping event | {} \n {} s in-track lag ({} km) | cross-track={} km').format(
                                    self.fb_id, self.ac_id, tRange[0].date(), 
                                    round(self.ac_time_lag, 1), 
                                    round(np.abs(self.ac_time_lag)*7.5, 1),
                                    round(start_cross_track, 2))
        ax[0].set_title(titleStr)
        ax[-1].set_xlabel('UTC')

        # Set xlims for all subplots
        ax[0].set_xlim([timedelta(seconds=self.fb_time_shift) + t for t in tRange])
        ax[1].set_xlim([timedelta(seconds=self.ac_time_lag) + t for t in tRange]) # +self.fb_time_shift
        #ax[2].set_xlim([timedelta(seconds=self.fb_time_shift) + t for t in tRange])
        #ax[2].legend()
        # for a in ax[1:]:
        #     a.set_xlim(tRange)
        for a in ax: # Format time stamps for all subplots
            myFmt = matplotlib.dates.DateFormatter('%H:%M:%S')
            a.xaxis.set_major_formatter(myFmt)
        return

    def plot_lap_events(self):
        """
        This method is similar to plot_lap_event, but it automatically
        sifts through any HiRes data avaliable between self.startDate 
        and self.endDate, and plots the AC6-FRIEBIRD data for those
        times, with no regard to their separation at that time. 
        """

        return

    def _load_fb_data(self, tRange):
        """ This method loads in the FIREBIRD-II data """
        hrName = 'FU{}_Hires_{}_L2.txt'.format(self.fb_id, tRange[0].date())
        self.hr = spacepy.datamodel.readJSONheadedASCII(
                            os.path.join(self.fbDir, hrName))
        self.hr['Time'] = spacepy.time.Ticktock(self.hr['Time']).UTC
        self.fb_time_shift = np.mean(self.hr['Count_Time_Correction'])
        return

    def _load_ac_data(self, tRange, dType='10Hz'):
        """ This method loads in the AC-6 data """
        self.acData = read_ac_data.read_ac_data_wrapper(self.ac_id,
            tRange[0], dType=dType, plot=False)
        return

    def _load_sep(self, fPath):
        """
        This method loads in the separation data file and saves it so self.sep
        """
        with open(fPath) as f:
            r = csv.reader(f)
            next(r) # Skip header
            rData = np.array(list(r))
        self.sep = {}
        self.sep['dateTime'] = np.array([dateutil.parser.parse(t) for t in rData[:, 0]])
        self.sep['d_in_track'] = np.array([float(f) for f in rData[:, 1]])
        self.sep['d_cross_track'] = np.array([float(f) for f in rData[:, 2]])
        return

    def _plot_fb(self, tRange, axCounts, axPos=None, axL=True):
        """ This method plots the FIREBIRD col counts data. """
        normTind = np.where((self.hr['Time'] > tRange[0]) & 
                            (self.hr['Time'] < tRange[1]))[0]
        shiftInd = np.where(
            (self.hr['Time'] > tRange[0]-timedelta(seconds=self.fb_time_shift)) & 
            (self.hr['Time'] < tRange[1]-timedelta(seconds=self.fb_time_shift)))[0]
        fbTimes = np.array([t+timedelta(seconds=self.fb_time_shift) 
                            for t in self.hr['Time'][normTind]])
        for E, Elabel in enumerate(self.hr['Col_counts'].attrs['ELEMENT_LABELS']):
            axCounts.plot(fbTimes, self.hr['Col_counts'][normTind, E],
                    label=Elabel)
        axCounts.set(ylabel='FU{} counts/bin'.format(self.fb_id), yscale='log')
        axCounts.legend()

        if axPos is not None:
            axPos.plot(self.hr['Time'][normTind], self.hr['Lat'][normTind], 
                        label='FU{} lat'.format(self.fb_id))
            axPos.set_ylabel('latitude [deg]')
        if axL:
            axL = axCounts.twinx()
            axL.plot(self.hr['Time'], np.abs(self.hr['McIlwainL']), 'k')
            axL.set_ylabel('McIlwain L (T89) (black curve)')
            axL.set_ylim(3, 10)
        return

    def _plot_ac(self, tRange, axCounts, axPos=None, axL=True):
        """
        This method plots the AC6 dosimiter count rates and position 
        in a similar way to _plot_fb()
        """
        # Plot dosimiter counts
        for key in ['dos1rate', 'dos2rate', 'dos3rate']:
            validCounts = np.where((self.acData[key] != -1E31) & 
                            (self.acData['dateTime'] > tRange[0]+timedelta(seconds=self.ac_time_lag)) & 
                            (self.acData['dateTime'] < tRange[1]+timedelta(seconds=self.ac_time_lag)) )[0]
            axCounts.plot(self.acData['dateTime'][validCounts], 
                        self.acData[key][validCounts],
                        label=key)
        axCounts.set_yscale('log')
        axCounts.set_ylabel('Dos rate [counts/s]')
        axCounts.legend()
        # Plot position
        if axPos is not None:
            shiftedTimes = np.array([timedelta(seconds=-self.ac_time_lag) + t
                            for t in self.acData['dateTime']])
            axPos.plot(shiftedTimes[validCounts], self.acData['lat'][validCounts], 
                        label='AC6-{} lat'.format(self.ac_id))
        if axL:
            axL = axCounts.twinx()
            validL = np.where(self.acData['Lm_OPQ'] != -1E31)[0]
            axL.plot(self.acData['dateTime'][validL], self.acData['Lm_OPQ'][validL], 'k')
            axL.set_ylabel('McIlwain L (OPQ) (black curve)')
            axL.set_ylim(3, 10)
        return

    def _get_ac6_lag(self, tRange):
        """
        This method returns the in-track at the start and end time
        specified in tRange.
        """
        idt = np.where((self.sep['dateTime'] > tRange[0]) & 
                        (self.sep['dateTime'] < tRange[1]))[0]
        start_in_track = self.sep['d_in_track'][idt[0]]/7.5
        end_in_track = self.sep['d_in_track'][idt[-1]]/7.5
        start_cross_track = self.sep['d_cross_track'][idt[0]]
        end_cross_track = self.sep['d_cross_track'][idt[-1]]
        return start_in_track, end_in_track, start_cross_track, end_cross_track

if __name__ == '__main__':
    acDtype = 'survey'
    fb_id = 4
    ac_id = 'A'
<<<<<<< Updated upstream
    dPath = './data/dist/2018-04-11_2018-06-11_FU{}_AC6{}_dist_v2.csv'.format(fb_id, ac_id)
    tRange = [datetime(2018, 4, 21, 11, 23), datetime(2018, 4, 21, 11, 29)]

    l = Lap(dPath, fb_id, ac_id)
    l.plot_lap_event(tRange, acDtype=acDtype, lag=244+40+11*60)
=======
    for fb_id in [3, 4]:
        START_DATE = datetime(2018, 7, 27)
        END_DATE = datetime(2018, 8, 31)
        #dPath = './data/dist/2018-04-11_2018-06-11_FU{}_AC6{}_dist_v2.csv'.format(fb_id, ac_id)
        #START_DATE = datetime(2018, 4, 19)
        #END_DATE = datetime.now()
>>>>>>> Stashed changes

    plt.tight_layout()
    plt.show()
    #plt.savefig('./plots/{}/{}_FU{}-AC6{}_{}_lap_event.png'.format(acDtype, lapTime, fb_id, ac_id, acDtype))