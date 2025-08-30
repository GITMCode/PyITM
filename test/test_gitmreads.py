#!/usr/bin/env python3

# This file is meant to test all the read routines we have available

# Nothing is actually done other than reading data in
import sys
sys.path.insert(0,'/home/ridley/Software/PyITM/')


from pyitm.fileio import util

def display_info(data):
    print('  keys : ', data.keys())
    print('  variables (read) : ', data['vars'])
    if ('data' in data.keys()):
        print('  shape of data (read) : ', data['data'].shape)
    print('   -> ntimes : ', data['ntimes'])
    if (data['nblocks'] > 0):
        print('   -> nblocks : ', data['nblocks'])
    if (data['nvars'] > 1):
        print('   -> nvars : ', data['nvars'])
    print('   -> nlons : ', data['nlons'])
    print('   -> nlats : ', data['nlats'])
    print('   -> nalts : ', data['nalts'])

    return


dat0 = util.read_all_files("data/2DANC_t021221_002000.bin", verbose=True)
print('Info on first file (dat0) : ')
display_info(dat0)


dat1 = util.read_all_files("data/3DALL_t021221_002000.bin", verbose=True, varsToRead = [3,2,1,6,17])
print('Info on second file (dat1) : ')
display_info(dat1)

header2 = util.read_all_headers("data/2DANC*.bin", verbose=True)
print('Info on third header files (header2) : ')
display_info(header2)

dat3 = util.read_all_files("data/2DG*", verbose=True)
print('Info on forth file (dat3) : ')
display_info(dat3)

dat3 = util.read_all_files("data/2DG*", verbose=True)
print('Info on forth file (dat3) : ')
display_info(dat3)

nc1 = util.read_all_files("data/3DALG_20110320_001000.nc")
print('Info on first nc file (nc1) : ')
display_info(nc1)

nc2 = util.read_all_files("data/3DALG_20110320_001000.nc", varsToRead = [15])
print('Info on first nc file (nc2, reading var = 15) : ')
display_info(nc2)

nc3 = util.read_all_files("data/3DALG_20110320_001000.nc", varsToRead = ['Tn'])
print('Info on first nc file (nc3, reading var = Tn) : ')
display_info(nc3)


# dat = gitmio.read_one_to_xr("data/2DGEL_t021221_002000.bin", verbose=True)
# dat = gitmio.read_multiple_to_xr("data/2DG*", verbose=True)
# for dat in [dat0, dat1, dat3]:
#     print([d for d in dat['data']])


print("good work chump")