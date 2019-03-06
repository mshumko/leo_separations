### This code makes prelim magnetic ephemeris ###
import csv
import numpy as np
import os
import dateutil.parser

import IRBEM

class AppendMagEphem(IRBEM.MagFields):
    def __init__(self, ephemPath, kext='T89'):
        IRBEM.MagFields.__init__(self, kext=kext)
        self.extModel = kext
        if 'reach' not in ephemPath:
            self._load_ephem(ephemPath)
        else:
            self._load_reach(ephemPath)
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

    def save_magephem(self, path):
        """ 
        This method appends the L and MLT to the ephemeris and saves it to path.
        """
        z = zip(self.eph['dateTime'], self.eph['Lat'], self.eph['Lon'],
                self.eph['Alt'], self.L, self.MLT)
        with open(path, 'w', newline='') as f:
            w = csv.writer(f)
            keys = ['Time (ISO)','Lat (deg)','Lon (deg)','Alt (km)',
                     'Lm_{}'.format(self.extModel), 
                     'MLT_{}'.format(self.extModel)]
            w.writerow(keys)

            for line in z:
                w.writerow([*line])
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
        
    def _load_reach(self, path):
        self.eph = pd.read_csv(path)
        date = path.split('.')[1]
        self.eph['dateTime'] = pd.to_datetime(date[0:4]) +
                               pd.to_timedelta(self.eph['DoY']-1, unit='D')
        self.eph = self.eph.set_index('dateTime')
        return

if __name__ == '__main__':
#    from datetime import datetime
#    for sc_id in ['FU3', 'FU4', 'ELFIN_A']:
#        START_DATE = datetime(2018, 12, 10)
#        END_DATE = datetime(2019, 1, 30)
#        ephemDir = './data/ephem'
#        ephemName = '{}_{}_{}_LLA_ephemeris.csv'.format(
#                       sc_id,
#                       START_DATE.date(),
#                       END_DATE.date())
#        magephemName = ephemName.split('_')
#        magephemName[-1] = 'magephem.csv'
#        magephemName.pop(-2) 
#        magephemName = '_'.join(magephemName)
#        a = AppendMagEphem(os.path.join(ephemDir, ephemName))
#        a.calc_magephem(maginput={'Kp':20})
#        a.save_magephem(os.path.join('./data/magephem/', magephemName))

    ephemDir = '/home/mike/research/reach/data'
    ephemName = 'reach.20190125.vid-169.txt'
    a = AppendMagEphem(os.path.join(ephemDir, ephemName))
    a.calc_magephem(maginput={'Kp':20})
    a.save_magephem(os.path.join('./data/magephem/', magephemName))
