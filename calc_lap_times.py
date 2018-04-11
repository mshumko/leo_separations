# This script calculates the good time intervals to target.
import dateutil.parser
import numpy as np
import csv

class LapTimes:
    def __init__(self, sepPath):
        self.sepData = self._load_sep(sepPath)
        return     

    def calcLapTimes(self, thresh=500):
        """ 
        This method calculates the start and end times when the two
        spacecraft where within thresh km separation.
        """
        sepInd = np.where(self.sepData['d'] < thresh)[0]
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
            w.writerow(['lapStartTime', 'lapEndTime', 'minDist [km]'])

            for (sT, eT, dMin) in zip(self.startTime, self.endTime, self.dmin):
                w.writerow([sT, eT, dMin])
            return

    def _calc_min_sep(self, startInd, endInd):
        """ 
        For each lapping event, this method calculates the closest separation.
        """
        self.dmin = np.nan*np.zeros(len(startInd))
        # To avoid a crash in case there is only one data point below threshold.
        # Should not effect the results.
        endInd += 1 
        for i, (sI, eI) in enumerate(zip(startInd, endInd)):
            self.dmin[i] = np.min(self.sepData['d'][sI:eI])
        return

    def _load_sep(self, path):
        """ This method loads in the separation file """
        with open(path) as f:
            r = csv.reader(f)
            keys = next(r)
            rawData = np.array(list(r))
        sepData = {}
        sepData['dateTime'] = np.array([dateutil.parser.parse(t)
                             for t in rawData[:, 0]])
        sepData['d'] = np.array([float(d) for d in rawData[:, 1]])  
        return sepData

if __name__ == '__main__':
    fb_id = 4
    L = LapTimes('/home/mike/research/leo-lapping-events/data/'
                '2018-04-11_2018-06-10_FU{}_AC6A_dist.csv'.format(fb_id))
    L.calcLapTimes()
    L.saveData('./data/FU{}_AC6A_lap_times.csv'.format(fb_id))
    