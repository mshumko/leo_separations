# This script calculates the good time intervals to target.
from matplotlib.dates import date2num
from datetime import timedelta
import dateutil.parser
import numpy as np
import csv

class LapTimes:
    def __init__(self, sc_a, sc_b, sepPath):
        self.sc_a = sc_a
        self.sc_b = sc_b
        self.sepData = self._load_sep(sepPath)
        return     

    def calcLapTimes(self, thresh=500):
        """ 
        This method calculates the start and end times when the two
        spacecraft where within thresh km separation.
        """
        sepInd = np.where(np.abs(self.sepData['d']) < thresh)[0]
        sepInd = np.append(sepInd, -9999)
        conv = np.convolve([1, -1], sepInd, mode = 'valid') - 1
        consecutiveFlag = np.where(conv != 0)[0] + 1
        startInd = np.insert(consecutiveFlag, 0, 0)[:-1]
        endInd = np.insert(consecutiveFlag, len(consecutiveFlag), len(sepInd)-1)[:-1]
        self.startTime = self.sepData['dateTime'][sepInd[startInd]]
        self.endTime = self.sepData['dateTime'][sepInd[endInd-1]]
        
        # This is a python indexing thing, but when start and end 
        # indicies are off by 1, the start/end time is the same.
        # So here I am arbitarily adding a minute.
        nStart = date2num(self.startTime)
        nEnd = date2num(self.endTime)
        i_same = np.where(np.isin(nEnd, nStart))[0]
        for i in i_same:
            self.endTime[i] += timedelta(minutes=1)
        # Calc lapping event duration (in minutes)    
        self.duration = np.array([(t_f - t_i).seconds/60 for t_i, t_f 
                        in zip(self.startTime, self.endTime)])
        
        # Calculate min separation
        self._calc_min_sep(sepInd[startInd], sepInd[endInd-1])
        return

    def saveData(self, fPath):
        with open(fPath, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['lapStartTime', 'lapEndTime', 
                        'lapDuration [min]', 'minDist [km]', 
                        '{}_L_at_min'.format(self.sc_a), 
                        '{}_L_at_min'.format(self.sc_b)])
            zTuple = zip(self.startTime, self.endTime, self.duration,
                         self.dmin, self.scALmin, self.scBLmin)
            for z in zTuple:
                w.writerow([*z])
            return

    def _calc_min_sep(self, startInd, endInd):
        """ 
        For each lapping event, this method calculates the closest separation.
        """
        self.dmin = np.nan*np.zeros(len(startInd))
        self.scALmin = np.nan*np.zeros(len(startInd))
        self.scBLmin = np.nan*np.zeros(len(startInd))
        # To avoid a crash in case there is only one data point below threshold.
        # Should not effect the results.
        endInd += 1 
        for i, (sI, eI) in enumerate(zip(startInd, endInd)):
            self.dmin[i] = np.min(self.sepData['d'][sI:eI])
            iMin = np.argmin(self.sepData['d'][sI:eI])+sI

            self.scALmin[i] = self.sepData['L_{}'.format(self.sc_a)][iMin]
            self.scBLmin[i] = self.sepData['L_{}'.format(self.sc_b)][iMin]
        return

    def _load_sep(self, path):
        """ This method loads in the separation file """
        with open(path) as f:
            r = csv.reader(f)
            keys = next(r)
            rawData = np.array(list(r))
        sepData = {}
        for (i, key) in enumerate(keys):
            if key == 'dateTime':
                sepData['dateTime'] = np.array([dateutil.parser.parse(t)
                                for t in rawData[:, i]])
            #elif 'dist_in_track' in key:
            #    sepData['d'] = np.array([float(d) for d in rawData[:, i]])  
            else:
                sepData[key] = np.array([float(d) for d in rawData[:, i]]) 
        sepData['d'] = np.sqrt(sepData['dist_in_track [km]']**2 + 
                        sepData['dist_cross_track [km]']**2)
        return sepData

if __name__ == '__main__':
    from datetime import date
    sc = ['FU3', 'ELFIN_A']
    dates = [date(2018, 10, 31), date(2019, 2, 1)]
    L = LapTimes(*sc, '/home/mike/research/leo-lapping-events/data/dist/'
                '{}_{}_{}_{}_dist_v2.csv'.format(*dates, *sc))
    L.calcLapTimes()
    L.saveData('./data/lap_times/{}_{}_{}_{}_lap_times.csv'.format(
                *dates, *sc))
    
