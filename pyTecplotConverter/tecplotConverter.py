import os, sys
sys.path.append(os.getcwd())

import numpy as np
import struct
from tecio_szl import open_file, zone_write_double_values,\
    zone_write_float_values, close_file, create_FE_zone, \
        zone_write_nodemap, FD_DOUBLE, FD_FLOAT, FEBRICK

#%% Function for displaying progress bar
def progressbar(progress, barLength=15):
    percent = float(progress) * 100 / barLength
    done = '#' * progress
    rem = '.' * (barLength - len(done))

    print('[%s%s] %d %%' % (done, rem, percent), end='\r')

#%% Variables
print('Setting up environment for datafile conversion: ')
progress = 0
bar = 100
progressbar(progress, barLength=bar)
nn = 2965732
ne = 9732564
curr = 0
tot = nn + ne
target = tot / bar 
ndof = 5
nvar = 7
gamma = 1.4
cv = 1/(gamma*(gamma-1))
x = np.empty(nn, dtype='float64')
y = np.empty(nn, dtype='float64')
z = np.empty(nn, dtype='float64')
ien = np.empty((ne,8), dtype='int32')
d = np.empty((nn,ndof), dtype='float64')
d_trans = np.empty((nn,nvar), dtype='float64')

#%% Filenames
start = 17
end = 20
dataFile = ['data.'+ str(i).zfill(2) for i in range(start,end+1)]
outFile = ['Re10k_' + str(i).zfill(2) + '.szplt' for i in range(start,end+1)]
#dataFile = ['data.28']
#outFile = ['Re10k_bp1p4_28000.szplt']          
xyzFile = '../mesh_info/mxyz'
ienFile = '../mesh_info/mien'

#%% Read xyz and ien data
xyzfile = open(xyzFile,'rb')
ienfile = open(ienFile,'rb')

for i in range(nn):
    xyz = xyzfile.read(24)
    x[i], y[i], z[i] = struct.unpack('ddd', xyz)
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar
        progressbar(progress, barLength=bar)
curr = nn

for i in range(ne):
    for j in range(8):
        iendata = ienfile.read(4)
        ien[i,j] = struct.unpack('i', iendata)[0]
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar
        progressbar(progress, barLength=bar)
curr += ne    
xyzfile.close()
ienfile.close()

progressbar(progress + 1, barLength=bar)
print('\n')

#%% Read and write data
for iFile, dataFilename in enumerate(dataFile):
    print('Converting ' + dataFilename + ' to ' + outFile[iFile] + ': ')
    tot = nn + 1
    curr = 0
    progress = 0
    target = tot / bar
    progressbar(progress, barLength=bar)
    datafile = open(dataFilename,'rb')
    for i in range(nn):
        for j in range(ndof):
            data = datafile.read(8)
            d[i,j] = struct.unpack('d', data)[0]
        d_trans[i,0] = d[i,0]
        d_trans[i,1] = d[i,1]/d[i,0]
        d_trans[i,2] = d[i,2]/d[i,0]
        d_trans[i,3] = d[i,3]/d[i,0]
        d_trans[i,4] = (1/cv)*(d[i,4]/d[i,0] - 0.5*(d_trans[i,1]**2 \
                     + d_trans[i,2]**2 + d_trans[i,3]**2))
        d_trans[i,5] = (d_trans[i,1]**2 + d_trans[i,2]**2 \
                     + d_trans[i,3]**2)**0.5/d_trans[i,4]**0.5
        d_trans[i,6] = d_trans[i,0]*d_trans[i,4]/gamma
        if i >= target - curr:
            progress += 1
            target = (progress + 1) * tot / bar
            progressbar(progress, barLength=bar)
    curr = nn
    datafile.close()
    varNames = ['x','y','z','rho','uvel','vvel','wvel','temp','mach','pressure']
    ndata = 3 + nvar
    file_handle = open_file(outFile[iFile], "Y-duct", varNames)
    #print(outFile[iFile] + ' opened.')
    zone = create_FE_zone(file_handle, "Zone", nn, ne, FEBRICK, None, \
                               [FD_FLOAT for i in range(len(varNames))])
    zone_write_float_values(file_handle, zone, 1, x)
    zone_write_float_values(file_handle, zone, 2, y)
    zone_write_float_values(file_handle, zone, 3, z)
    for i in range(4, len(varNames) + 1):
        zone_write_float_values(file_handle, zone, i, d_trans[:, i - 4])
    zone_write_nodemap(file_handle, zone, ien.reshape(-1))
    close_file(file_handle)
    progressbar(progress + 1, barLength=bar)
    print('\n')
    #print(outFile[iFile] + ' closed.')


