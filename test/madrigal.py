from pyitm.fileio.madrigalio import read_madrigal_file

data = read_madrigal_file('data/madrigal/gps021221g.002.nc', verbose=True)
for k in data.keys():
    print(k, data[k].shape)


data = read_madrigal_file('data/madrigal/dms_20150716_18s1.001.hdf5', verbose=True)
for k in data.keys():
    print(k, data[k].shape)
print(data['times'][0:10])

data = read_madrigal_file('data/cul130621.17201.hdf5', verbose=True)
for k in data.keys():
    print(k, data[k].shape)
print(data['times'][0:10])
