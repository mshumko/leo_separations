# This script calculates how often close lapping events occur.

import numpy as np
import csv
import dateutil.parser
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

fname = './data/2018-04-11_2018-06-11_FU4_AC6A_lap_times_500km_thresh_v2.csv'

with open(fname) as f:
    r = csv.reader(f)
    keys = next(r)
    rawData = np.array(list(r))

lapData = {}
for i, key in enumerate(keys):
    if 'Time' in key:
        lapData[key] = np.array([dateutil.parser.parse(t) 
                                for t in rawData[:, i]])
    else:
        lapData[key] = np.array(rawData[:, i], dtype=float)
z = zip(lapData['lapStartTime'][1:], lapData['lapStartTime'][:-1])
dt = [(tii - ti).total_seconds() for (tii, ti) in z]

plt.hist(dt, bins=np.linspace(28000, 34000))
plt.show()