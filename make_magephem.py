### This code makes prelim magnetic ephemeris ###
import csv
import numpy as np
import os
import dateutil.parser

import pandas as pd

import IRBEM

class AppendMagEphem(IRBEM.MagFields):
    def __init__(self, ephemPath, kext='T89'):
        IRBEM.MagFields.__init__(self, kext=kext)
        self.extModel = kext
        self.ephemPath = ephemPath
        # if 'reach' not in ephemPath:
        #     self._load_ephem(ephemPath)
        # else:
        #     self._load_reach(ephemPath)
        self.load_ephem(ephemPath)
        return

    def calc_magephem(self, maginput=None):
        """
        This method loops over the epehem times and calculates L and MLT
        """

        self.L = np.nan*np.ones(len(self.eph['Alt']))
        self.MLT = np.nan*np.ones(len(self.eph['Alt']))

        z = zip(self.eph['dateTime'].dt.to_pydatetime(),
                self.eph['Alt'], self.eph['Lat'], 
                self.eph['Lon'])
        
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
            keys = ['dateTime','Lat','Lon','Alt',
                     'Lm_{}'.format(self.extModel), 
                     'MLT_{}'.format(self.extModel)]
            w.writerow(keys)

            for line in z:
                w.writerow([*line])
        return

    def load_ephem(self, path):
        """
        This method reads in the ephemeris file.
        """
        self.eph = pd.read_csv(path)
        # Convert times
        self._convert_times(self.eph)
        # Convert lat/lon/alt
        self._convert_lla(self.eph)
        return

    def _convert_times(self, data):
        """
        This helper method intelegently converts the time column(s)
        to datetime objects. If a column with the word 'time' is found,
        ity will be converted to datetime objects. Otherwise if 
        no such column is found, it will attempt to combine other 
        keys such as 'year', 'month', 'day' etc. and convert that.
        """
        keys = data.keys()
        time_keys = [key for key in keys if 'time' in key.lower()]

        if len(time_keys) == 1:
            # Easy case
            data['dateTime'] = pd.to_datetime(data[time_keys[0]])

        elif len(time_keys) == 0 and ('reach' in self.ephemPath.lower()):
            # Harder case with reach.
            date = os.path.basename(self.ephemPath).split('.')[1]
            data['dateTime'] = pd.to_datetime(date[0:4]) + \
                               pd.to_timedelta(data['DoY']-1, unit='D')
        else:
            raise ValueError(f'More than one time key found!'
                            f'\ntime_keys={time_keys}')
        return data

    def _convert_lla(self, data):
        """ 
        This method attempts to convert lat/lon/alt values in a conistant manner.
        """
        # Find all instances of LLA keys in the data keys.
        keys = data.keys()
        alt_key = [key for key in keys if ('alt' in key.lower())]
        lat_key = [key for key in keys if (('lat' in key.lower()) 
                    and ('inv' not in key.lower()))]
        lon_key = [key for key in keys if ('lon' in key.lower())]

        # Check that it correctly picked out unique keys.
        if len(alt_key) != 1 or len(lat_key) != 1 or len(lon_key) != 1:
            raise ValueError(f'Incorrect number of LLA keys found!\n'
                            f'alt_key={alt_key}, lat_key={lat_key}, '
                            f'lon_key={lon_key}')

        data['Alt'] = data[alt_key[0]]
        data['Lat'] = data[lat_key[0]]
        data['Lon'] = data[lon_key[0]]
        return data


if __name__ == '__main__':
#   from datetime import datetime
#    for sc_id in ['FU3', 'FU4']:
#        START_DATE = datetime(2019, 1, 21)
#        END_DATE = datetime(2019, 2, 23)
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
    magephemName = ephemName.split('.')
    magephemName.insert(-1, 'magephem')
    magephemName = '.'.join(magephemName)
    print(magephemName)
    a = AppendMagEphem(os.path.join(ephemDir, ephemName))
    a.calc_magephem(maginput={'Kp':20})
    a.save_magephem(os.path.join('./data/magephem/', magephemName))
