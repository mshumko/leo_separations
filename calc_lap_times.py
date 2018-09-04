# This script calculates the good time intervals to target.
import dateutil.parser
import numpy as np
import csv

class LapTimes:
    def __init__(self, fb_id, sepPath):
        self.fb_id = fb_id
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
        # Calculate min separation
        self._calc_min_sep(sepInd[startInd], sepInd[endInd-1])
        return

    def saveData(self, fPath):
        with open(fPath, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['lapStartTime', 'lapEndTime', 'minDist [km]', 
                        'FU{}_L_at_min'.format(self.fb_id), 'AC6A_L_at_min'])

            for z in zip(self.startTime, self.endTime, self.dmin, self.FBLmin, self.AC6Lmin):
                w.writerow([*z])
            return

    def _calc_min_sep(self, startInd, endInd):
        """ 
        For each lapping event, this method calculates the closest separation.
        """
        self.dmin = np.nan*np.zeros(len(startInd))
        self.AC6Lmin = np.nan*np.zeros(len(startInd))
        self.FBLmin = np.nan*np.zeros(len(startInd))
        # To avoid a crash in case there is only one data point below threshold.
        # Should not effect the results.
        endInd += 1 
        for i, (sI, eI) in enumerate(zip(startInd, endInd)):
            self.dmin[i] = np.min(self.sepData['d'][sI:eI])
            iMin = np.argmin(self.sepData['d'][sI:eI])+sI

            self.AC6Lmin[i] = self.sepData['L_AC6A'][iMin]
            self.FBLmin[i] = self.sepData['L_FU{}'.format(self.fb_id)][iMin]
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
    from datetime import datetime
    fb_id = 4
    dates = [datetime(2018, 7, 27).date(), datetime(2018, 8, 31).date()]
    L = LapTimes(fb_id, '/home/mike/research/leo-lapping-events/data/dist/'
                '{}_{}_FU{}_AC6A_dist_v2.csv'.format(*dates, fb_id))
    L.calcLapTimes()
    L.saveData('./data/{}_{}_FU{}_AC6A_lap_times.csv'.format(*dates, fb_id))
    