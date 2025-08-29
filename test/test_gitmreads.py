#!/usr/bin/env python3

# This file is meant to test all the read routines we have available

# Nothing is actually done other than reading data in

from pyitm.fileio import gitmio


dat0 = gitmio.read_gitm_one_file("data/2DANC_t021221_002000.bin", verbose=True)

dat1 = gitmio.read_gitm_one_file("data/3DALL_t021221_002000.bin", verbose=True, varlist=[3,2,1,6,17])

dat = gitmio.read_gitm_headers("data/2DANC*.bin", verbose=True)

dat3 = gitmio.read_gitm_all_files("data/2DG*", verbose=True)

# dat = gitmio.read_one_to_xr("data/2DGEL_t021221_002000.bin", verbose=True)
# dat = gitmio.read_multiple_to_xr("data/2DG*", verbose=True)
print(dat1['data'][3].shape)
# for dat in [dat0, dat1, dat3]:
#     print([d for d in dat['data']])


print("good work chump")