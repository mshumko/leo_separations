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

        fig, ax = plt.subplots(3, figsize=(8, 9))
        self._plot_fb(tRange, ax[0], axPos=ax[2])
        self._plot_ac(tRange, ax[1], axPos=ax[2])

        ### Plot Adjustments ###
        titleStr = ('FU{} - AC6{} Lapping event | {} | {} s lag ({} km)').format(
                                    self.fb_id, self.ac_id, tRange[0].date(), 
                                    round(np.abs(self.ac_time_lag), 1), 
                                    round(np.abs(self.ac_time_lag)*7.5, 1))
        ax[0].set_title(titleStr)
        ax[-1].set_xlabel('UTC')

        # Set xlims for all subplots
        ax[0].set_xlim([timedelta(seconds=self.fb_time_shift) + t for t in tRange])
        ax[1].set_xlim([timedelta(seconds=self.ac_time_lag+self.fb_time_shift) + t for t in tRange])
        ax[2].set_xlim([timedelta(seconds=self.fb_time_shift) + t for t in tRange])
        ax[2].legend()
        # for a in ax[1:]:
        #     a.set_xlim(tRange)
        for a in ax: # Format time stamps for all subplots
            myFmt = matplotlib.dates.DateFormatter('%H:%M:%S')
            a.xaxis.set_major_formatter(myFmt)
        return

    def _load_fb_data(self, tRange):
        """ This method loads in the FIREBIRD-II data """
        hrName = 'FU{}_Hires_{}_L2.txt'.format(self.fb_id, tRange[0].date())
        self.hr = spacepy.datamodel.readJSONheadedASCII(
                            os.path.join(self.fbDir, hrName))#.copy()
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
        # Find the first order correction (unsigned)
        sepInd = np.where(self.sep['dateTime'] > tRange[0])[0][0]
        self.ac_time_lag = self.sep['d'][sepInd]/7.5
        #print('Zeroth order time lag', self.ac_time_lag)
        # Determine the sign of self.ac_time_lag
        dLatPast = np.where(self.acData['dateTime'] > 
                    tRange[0] - timedelta(seconds=self.ac_time_lag))[0][0]
        dLatFuture = np.where(self.acData['dateTime'] > 
                    tRange[0] + timedelta(seconds=self.ac_time_lag))[0][0]
        fbSind = np.where(self.hr['Time'] > tRange[0])[0][0]
        pastDiff = np.abs(self.hr['Lat'][fbSind]-self.acData['lat'][dLatPast])
        futureDiff = np.abs(self.hr['Lat'][fbSind]-self.acData['lat'][dLatFuture])
        if pastDiff > futureDiff:
            self.ac_time_lag *= -1
            acInd = dLatPast
        else:
            acInd = dLatFuture
        
        # Now find the higher order correction (signed)
        self.ac_time_lag = ( self.hr['Time'][fbSind] - 
            self.acData['dateTime'][acInd] ).total_seconds()
        #print('First order time lag', self.ac_time_lag)
        #print(self.hr['Lat'][fbSind], self.acData['lat'][dLatPast], 
        #            self.acData['lat'][dLatFuture])

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
            axPos.plot(self.hr['Time'][normTind], self.hr['Lon'][normTind], 
                        label='FU{} Lon'.format(self.fb_id))
            axPos.set_ylabel('Longitude [deg]')
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
        axCounts.legend()
        # Plot position
        if axPos is not None:
            shiftedTimes = np.array([timedelta(seconds=-self.ac_time_lag) + t
                            for t in self.acData['dateTime']])
            axPos.plot(shiftedTimes[validCounts], self.acData['lon'][validCounts], 
                        label='AC6-{} lon'.format(self.ac_id))
        if axL:
            axL = axCounts.twinx()
            validL = np.where(self.acData['Lm_OPQ'] != -1E31)[0]
            axL.plot(self.acData['dateTime'][validL], self.acData['Lm_OPQ'][validL], 'k')
            axL.set_ylabel('McIlwain L (OPQ) (black curve)')
            axL.set_ylim(3, 10)
        
        return

def getLapTimes(fb_id, ac_id):
    if fb_id == 3 and ac_id == 'A':
        lapTimes = {
                # March 19th
                '20180319T0048':{'tRange':[datetime(2018, 3, 19, 0, 48, 41), 
                                        datetime(2018, 3, 19, 0, 54, 0)]
                                        },
                '20180319T0104':{'tRange':[datetime(2018, 3, 19, 1, 3, 30), 
                                        datetime(2018, 3, 19, 1, 6, 0)]
                                        },
                # March 25th
                '20180325T0105':{'tRange':[datetime(2018, 3, 25, 1, 5, 51), 
                                        datetime(2018, 3, 25, 1, 8, 30)]
                                        },
                # March 26th
                '20180326T0038':{'tRange':[datetime(2018, 3, 26, 0, 38, 0), 
                                        datetime(2018, 3, 26, 0, 44, 0)]
                                        },
                '20180326T0051':{'tRange':[datetime(2018, 3, 26, 0, 52, 0), 
                                        datetime(2018, 3, 26, 0, 53, 30)]
                                        },
                '20180326T0213':{'tRange':[datetime(2018, 3, 26, 2, 13), 
                                        datetime(2018, 3, 26, 2, 17)]
                                        },
                '20180326T0348':{'tRange':[datetime(2018, 3, 26, 3, 48, 25), 
                                        datetime(2018, 3, 26, 3, 52)]
                                        },
                # # March 28th
                # '0328T1257':{'tRange':[datetime(2018, 3, 28, 12, 57, 40), 
                #                         datetime(2018, 3, 28, 13, 1, 40)]
                #                         },          
                }         
    elif fb_id == 4 and ac_id == 'A':
        lapTimes= {
                # February 27th
                '20180227T0942':{'tRange':[datetime(2018, 2, 27, 9, 42),
                                        datetime(2018, 2, 27, 9, 38, 27)]
                                },
                '20180227T1033':{'tRange':[datetime(2018, 2, 27, 10, 33),
                                        datetime(2018, 2, 27, 10, 35, 27)]
                                },
                # March 1st
                '20180301T1919':{'tRange':[datetime(2018, 3, 1, 19, 19),
                                        datetime(2018, 3, 1, 19, 24, 0)]
                                }

                }
    else:
        raise NotImplementedError('No times dict for this spacecraft selection!')
    return lapTimes

if __name__ == '__main__':
    acDtype = 'survey'
    fb_id = 4
    ac_id = 'A'
    dPath = './data/2018-02-26_2018-03-29_FU{}_AC6{}_dist.csv'.format(fb_id, ac_id)

    for lapTime, value in getLapTimes(fb_id, ac_id).items():
        print('Making plot for', lapTime)
        l = Lap(dPath, fb_id, ac_id)
        l.plot_lap_event(value['tRange'], acDtype=acDtype)

        plt.tight_layout()
        plt.savefig('./plots/{}/{}_FU{}-AC6{}_{}_lap_event.png'.format(acDtype, lapTime, fb_id, ac_id, acDtype))
        #plt.show()