### This code makes prelim magnetic ephemeris ###
import csv
import numpy as np
import os
import dateutil.parser

import IRBEM

class AppendMagEphem(IRBEM.MagFields):
    def __init__(self, ephemPath, kext='T89'):
        IRBEM.MagFields.__init__(self, kext=kext)
        self._load_ephem(ephemPath)
        return

    def _load_ephem(self, path):
        """
        This code reads in the magnetic ephemeris file.
        """
        with open(path) as f:
            r = csv.reader(f)
            keys = next(r)
            rawData = np.array(list(r))

        self.eph = {}
        for (i, key) in enumerate(keys):
            if 'Time' in key:
                self.eph['dateTime'] = np.array([dateutil.parser.parse(t) 
                                        for t in rawData[:, i]])
            elif 'Lat' in key:
                self.eph['Lat'] = np.array([float(i) for i in rawData[:, i]])
            elif 'Lon' in key:
                self.eph['Lon'] = np.array([float(i) for i in rawData[:, i]])
            elif 'Alt' in key:
                self.eph['Alt'] = np.array([float(i) for i in rawData[:, i]])
        return

    def calc_magephem(self, maginput=None):
        """
        This method loops over the epehem times and calculates L and MLT
        """

        self.L = np.nan*np.ones(len(self.eph['Alt']))
        self.MLT = np.nan*np.ones(len(self.eph['Alt']))

        z = zip(self.eph['dateTime'], self.eph['Alt'],
                self.eph['Lat'], self.eph['Lon'])
        
        for i, (time, alt, lat, lon) in enumerate(z):
            X = {'dateTime':time, 'x1':alt, 'x2':lat, 'x3':lon}
            self.make_lstar(X, maginput)

            self.L[i] = self.make_lstar_output['Lm'][0]
            self.MLT[i] = self.make_lstar_output['MLT'][0]
        return

if __name__ == '__main__':
    ephemDir = '/home/mike/research/mission-tools/orbit/data/'
    ephemName = 'AEROCUBE_6A_2018-04-11_2018-06-11_LLA_ephemeris_pre.csv'
    a = AppendMagEphem(os.path.join(ephemDir, ephemName))
    a.calc_magephem(maginput={'Kp':20})