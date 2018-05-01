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
        if not hasattr(self, 'hr'):
            self._load_fb_data(tRange)
            print('Loading single day HiRes')
        self._load_ac_data(tRange, acDtype)
        
        # Only implement the lag for start of run, and implement 
        # the end time later.
        self._get_bounds(tRange) 
        if lag is not None:
            self.ac_time_lag = lag
        in_track_lag = (self.fbBounds[0] - self.ac6Bounds[0]).total_seconds()

        fig, ax = plt.subplots(2, figsize=(8, 9))
        self._plot_fb(tRange, ax[0])
        self._plot_ac(tRange, ax[1])
        
        ### Plot Adjustments ###
        titleStr = ('FU{} - AC6{} Lapping event | {} \n {} s in-track lag ({} km)').format(
                                    self.fb_id, self.ac_id, tRange[0].date(), 
                                    round(in_track_lag, 1), 
                                    round(np.abs(in_track_lag)*7.5, 1)
                                    )
        ax[0].set_title(titleStr)
        ax[-1].set_xlabel('UTC [hh:mm]')

        # Set xlims for all subplots
        ax[0].set_xlim(*self.fbBounds)
        ax[1].set_xlim(*self.ac6Bounds)
        
        for a in ax: # Format time stamps for all subplots
            myFmt = matplotlib.dates.DateFormatter('%H:%M')
            a.xaxis.set_major_formatter(myFmt)
        return

    def plot_lap_events(self, saveDir=None, acDtype='10Hz'):
        """
        This method is similar to plot_lap_event, but it automatically
        sifts through any HiRes data avaliable between self.startDate 
        and self.endDate, and plots the AC6-FRIEBIRD data for those
        times, with no regard to their separation at that time. 
        """
        self._load_all_fb_data()

        if saveDir is None:
            saveDir = '/home/mike/research/leo-lapping-events/plots/{}/{}/'.format(
                            datetime.now().date(), acDtype)
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
            print('Made directory at', saveDir)

        # Find all of the start/end times.
        dt = (self.hr['Time'][1:] - self.hr['Time'][:-1])
        tJump = np.where(dt > timedelta(minutes=1))[0]
        idt = sorted(np.concatenate((tJump, tJump+1, [0], [len(self.hr['Time'])-1] )))

        # Now loop over the HiRes times and call the plot_lap_event() function.
        for (i, j) in zip(idt[::2], idt[1::2]):
            self.fb_time_shift = self.hr['Count_Time_Correction'][i]
            self.plot_lap_event([self.hr['Time'][i], self.hr['Time'][j]], acDtype=acDtype)
            saveDate = self.hr['Time'][i].isoformat().replace(':', '').replace('-', '')
            saveName = '{}_FU{}_AC6{}_lap.png'.format(saveDate, self.fb_id, self.ac_id)
            plt.savefig(os.path.join(saveDir, saveName))
        return

    def _load_fb_data(self, tRange):
        """ This method loads in the FIREBIRD-II data """
        hrName = 'FU{}_Hires_{}_L2.txt'.format(self.fb_id, tRange[0].date())
        self.hr = spacepy.datamodel.readJSONheadedASCII(
                            os.path.join(self.fbDir, hrName))
        self.hr['Time'] = spacepy.time.Ticktock(self.hr['Time']).UTC
        self.fb_time_shift = np.mean(self.hr['Count_Time_Correction'])
        return

    def _load_all_fb_data(self):
        """
        This method will load in the HiRes data betwee startDate
        """
        # Load in the FIREBIRD HiRes data between specified time range, and find all HiRes times.
        days = [self.startDate + timedelta(days=i) for i in range((self.endDate-self.startDate).days)]
        self.hr = {'Time':np.array([]), 'Count_Time_Correction':np.array([]), 
                    'Col_counts':np.nan*np.ones((0, 6), dtype=int), 'Lat':np.array([]), 
                    'Lon':np.array([]), 'McIlwainL':np.array([]), 'MLT':np.array([])}
        for day in days:
            hrName = 'FU{}_Hires_{}_L2.txt'.format(self.fb_id, day.date())
            try:
                hrTemp = spacepy.datamodel.readJSONheadedASCII(
                                    os.path.join(self.fbDir, hrName))
            except FileNotFoundError: # If no file found, move on.
                continue

            hrTemp['Time'] = spacepy.time.Ticktock(hrTemp['Time']).UTC

            for key in filter(lambda x: x != 'Col_counts', self.hr):
                self.hr[key] = np.append(self.hr[key], hrTemp[key])      
            self.hr['Col_counts'] = np.concatenate((self.hr['Col_counts'], hrTemp['Col_counts']))
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

    def _plot_fb(self, tRange, axCounts, axL=True):
        """ This method plots the FIREBIRD col counts data. """
        normTind = np.where((self.hr['Time'] > tRange[0]) & 
                            (self.hr['Time'] < tRange[1]))[0]
        shiftInd = np.where(
            (self.hr['Time'] > tRange[0]-timedelta(seconds=self.fb_time_shift)) & 
            (self.hr['Time'] < tRange[1]-timedelta(seconds=self.fb_time_shift)))[0]
        fbTimes = np.array([t+timedelta(seconds=self.fb_time_shift) 
                            for t in self.hr['Time'][normTind]])
        for E in range(6):
            axCounts.plot(fbTimes, self.hr['Col_counts'][normTind, E],
                    label='ch{}'.format(E))
        axCounts.set(ylabel='FU{} counts/bin'.format(self.fb_id), yscale='log')
        axCounts.legend()
        if axL:
            axL = axCounts.twinx()
            axL.plot(self.hr['Time'], np.abs(self.hr['McIlwainL']), 'k')
            axL.set_ylabel('McIlwain L (T89) (black curve)')
            axL.set_ylim(3, 17)
        return

    def _plot_ac(self, tRange, axCounts, axL=True):
        """
        This method plots the AC6 dosimiter count rates and position 
        in a similar way to _plot_fb()
        """
        # Plot dosimiter counts
        for key in ['dos1rate', 'dos2rate', 'dos3rate']:
            validCounts = np.where((self.acData[key] != -1E31) & 
                            (self.acData['dateTime'] > self.ac6Bounds[0]) & 
                            (self.acData['dateTime'] < self.ac6Bounds[1]) )[0]
            axCounts.plot(self.acData['dateTime'][validCounts], 
                        self.acData[key][validCounts],
                        label=key)
        axCounts.set_yscale('log')
        axCounts.set_ylabel('Dos rate [counts/s]')
        axCounts.legend()
        # Plot position 
        if axL:
            axL = axCounts.twinx()
            validL = np.where(self.acData['Lm_OPQ'] != -1E31)[0]
            axL.plot(self.acData['dateTime'][validL], self.acData['Lm_OPQ'][validL], 'k')
            axL.set_ylabel('McIlwain L (OPQ) (black curve)')
            axL.set_ylim(3, 17)
        return

    def _get_bounds(self, tRange, thresh=180):
        """
        This method uses the in-track lag value along with the MagEphem
        file to match up L shells and return times when AC6 crossed the
        same L shells as FIREBIRD did (defined by tRange). I should soon
        implement the OPQ model for FIREBIRD for a direct model comparison. 
        """
        # Calculate FIREBIRD start/end L shells
        fbIdt = np.where((self.hr['Time'] > tRange[0]) & 
                        (self.hr['Time'] < tRange[1]))[0]
        for i in fbIdt:
            if np.abs(self.hr['McIlwainL'][i]) != 1E31:
                fbStartL = self.hr['McIlwainL'][i]
                fbStartI = i
                break
        for i in reversed(fbIdt):
            if np.abs(self.hr['McIlwainL'][i]) != 1E31:
                fbEndL = self.hr['McIlwainL'][i]
                fbEndI = i
                break
        self.fbBounds = [self.hr['Time'][fbStartI], self.hr['Time'][fbEndI]]
                
        # Calculate the in-track lag as a first guess for the AC6 times.
        idt = np.where((self.sep['dateTime'] > tRange[0]) & 
                        (self.sep['dateTime'] < tRange[1]))[0]
        tLag = self.sep['d_in_track'][idt[0]]/7.5
        # Get AC6 L shells around this time with a window.
        id6t = np.where((self.acData['dateTime'] > tRange[0] + timedelta(seconds=tLag) - 
                        timedelta(seconds=thresh)) & 
                        (self.acData['dateTime'] < tRange[1] + timedelta(seconds=tLag) +
                        timedelta(seconds=thresh)))[0]
        ac6L = self.acData['Lm_OPQ'][id6t]
        ac6L[ac6L == -1E31] = np.nan # Replace IRBEM error values with nans.
        # Calculate where AC6 L crosses FIREBIRD's L shells.
        acStartL = np.nanargmin(np.abs(ac6L - fbStartL))
        acEndL = np.nanargmin(np.abs(ac6L - fbEndL))
        self.ac6Bounds = [self.acData['dateTime'][id6t[0] + acStartL], 
                          self.acData['dateTime'][id6t[0] + acEndL]]
        
        # If the difference in the bounds is > 1 (no AC6 data to that high of L 
        # shell, then recalculate the FIREBIRD bounds
        if np.abs(self.acData['Lm_OPQ'][id6t[0]+acStartL] - fbStartL) > 1:
            fbStartI = np.nanargmin(np.abs(self.hr['McIlwainL'][fbIdt] - 
                                self.acData['Lm_OPQ'][id6t[0]+acStartL]))
            self.fbBounds[0] = self.hr['Time'][fbStartI+fbIdt[0]] 
            
        if np.abs(self.acData['Lm_OPQ'][id6t[0]+acEndL] - fbEndL) > 1:
            fbEndI = np.nanargmin(np.abs(self.hr['McIlwainL'][fbIdt] - 
                                self.acData['Lm_OPQ'][id6t[0]+acEndL]))
            self.fbBounds[1] = self.hr['Time'][fbEndI+fbIdt[0]]
               
        return [self.acData['dateTime'][id6t[0] + acStartL], self.acData['dateTime'][id6t[0] + acEndL]]

if __name__ == '__main__':
    acDtype = 'survey'
    fb_id = 3
    ac_id = 'A'
    dPath = './data/dist/2018-04-11_2018-06-11_FU{}_AC6{}_dist_v2.csv'.format(fb_id, ac_id)
    tRange = [datetime(2018, 4, 21, 11, 23), datetime(2018, 4, 21, 11, 29)]
    START_DATE = datetime(2018, 4, 19)
    END_DATE = datetime.now()

    l = Lap(dPath, fb_id, ac_id, startDate=START_DATE, endDate=END_DATE)
    l.plot_lap_events(acDtype='survey')
    #l.plot_lap_events(acDtype='10Hz')
    #l.plot_lap_event(tRange, acDtype=acDtype, lag=244+40+11*60)

    #plt.tight_layout()
    #plt.show()
    #plt.savefig('./plots/{}/{}_FU{}-AC6{}_{}_lap_event.png'.format(acDtype, lapTime, fb_id, ac_id, acDtype))
