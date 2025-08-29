

from pyitm.modeldata import satellite
from pyitm.fileio import gitmio
import datetime
import numpy as np

dat = gitmio.read_gitm_all_files('data/3DALL*')
dat['time'] = np.array(dat['times'])

dt_sat = 30 #seconds
t0sat = 5 # min before first gitm output
t1sat = 5 # min after last gitm output
# ^ set either of these negative to not have the satellite span outside the model times


sats = {'time': [],
        'lat':np.linspace(-90,90, 500),
        'lon':np.zeros(500)+120,
        'alt':np.zeros(500)+350}

sattime = dat['time'][0] - datetime.timedelta(minutes=t0sat)
# end time:
satendtime = dat['time'][-1] + datetime.timedelta(minutes=t0sat)

while sattime < satendtime:
    sats['time'].append(sattime)
    sattime += datetime.timedelta(seconds=dt_sat)

a, b = satellite.extract_1d(sats, dat, False, True, interpVar=34)
a, b = satellite.extract_1d(sats, dat, False, True)

print(a, b)
