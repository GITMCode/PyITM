

from bin import plot_satellite
from pyitm.modeldata import satellite
from pyitm.fileio import gitmio, satelliteio
import datetime
import glob
import numpy as np

dat = gitmio.read_gitm_all_files('data/3DALL_t0212*')
dat['times'] = np.array(dat['times'])

dt_sat = 30 # seconds
t0sat = 5 # min before first gitm output
t1sat = 5 # min after last gitm output
# ^ set either of these negative to not have the satellite span outside the model times


sats = {'times': [],
        'lats':np.linspace(-90,90, 500),
        'lons':np.zeros(500)+120,
        'alts':np.zeros(500)+350}

sattime = dat['times'][0] - datetime.timedelta(minutes=t0sat)
# end time:
satendtime = dat['times'][-1] + datetime.timedelta(minutes=t0sat)

while sattime < satendtime:
    sats['times'].append(sattime)
    sattime += datetime.timedelta(seconds=dt_sat)

a1 = satellite.extract_1d(sats, dat, False, True, interpVar=34)
a2 = satellite.extract_1d(sats, dat, False, True)

print(a1.keys())
print(a2.keys())

sat1 = satelliteio.read_sat_file('data/satfiles/CH_DNS_ACC_2002_12_v01.txt', verbose=True)
print(sat1.keys())

a3 = satellite.extract_1d(sat1, dat, verbose=False)
for k in a3.keys():
    print(k, len(a3[k]))




plot_satellite.main(satfiles=['data/satfiles/CH_DNS_ACC_2002_12_v01.txt'],
                    modeldatapath='data/3DALL_t0212*',
                    vars2plot='rho',
                    verbose=True)


plot_satellite.main(satfiles=['data/satfiles/CH_DNS_ACC_2002_12_v01.txt'],
                    modeldatapath='data/3DALL_t0212*',
                    vars2plot=3,
                    verbose=True)

newfiles = glob.glob('champ*')

print("-> Looks like I made some new files. You can `rm champ_*` to delete:")

print('>>>>>   ', newfiles)
