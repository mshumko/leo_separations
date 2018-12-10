# This is a wrapper to process the lap times completely. 
from datetime import datetime
import os

import calc_dist
import calc_lap_times
import make_magephem

START_DATE = datetime(2018, 12, 10)
END_DATE = datetime(2019, 1, 30)

sc_a_arr = ['FU3', 'FU4']
sc_b_arr = ['ELFIN_A']

print('Making magnetic ephemeris')
for sc_id in sc_a_arr + sc_b_arr:
    ephemDir = './data/ephem'
    ephemName = '{}_{}_{}_LLA_ephemeris.csv'.format(
                   sc_id,
                   START_DATE.date(),
                   END_DATE.date())
    magephemName = ephemName.split('_')
    magephemName[-1] = 'magephem.csv'
    magephemName.pop(-2) 
    magephemName = '_'.join(magephemName)
    a = make_magephem.AppendMagEphem(os.path.join(ephemDir, ephemName))
    a.calc_magephem(maginput={'Kp':20})
    a.save_magephem(os.path.join('./data/magephem/', magephemName))
    
print('Calculating separations')
for a_id in sc_a_arr:
    for b_id in sc_b_arr:
        pathA = ('./data/magephem/{}_{}_{}_magephem.csv'.format(
                a_id, START_DATE.date(), END_DATE.date()))
        pathB = ('./data/magephem/{}_{}_{}_magephem.csv'.format(
                b_id, START_DATE.date(), END_DATE.date()))
        saveDir = ('./data/dist/{}_{}_{}_{}_dist_v2.csv'.format(
                    START_DATE.date(), END_DATE.date(), a_id, b_id))
        c = calc_dist.CalcDist(a_id, b_id, START_DATE, END_DATE, pathA, pathB)
        c.calc_dist()
        c.save_file(saveDir)
        #c.plot_dist()
        #plt.show()
        
print('Calculating lap times')
for a_id in sc_a_arr:
    for b_id in sc_b_arr:
        L = calc_lap_times.LapTimes(a_id, b_id, '/home/mike/research/leo-lapping-events/data/dist/'
                    '{}_{}_{}_{}_dist_v2.csv'.format(START_DATE.date(), END_DATE.date(), a_id, b_id))
        L.calcLapTimes()
        L.saveData('./data/lap_times/{}_{}_{}_{}_lap_times.csv'.format(
                    START_DATE, END_DATE, a_id, b_id))
