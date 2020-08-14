import numpy as np
import struct
from tecio_szl import open_file, zone_write_double_values,\
    zone_write_float_values, close_file, create_FE_zone, \
        zone_write_nodemap, FD_DOUBLE, FD_FLOAT, FEBRICK

#%% Variables
nn = 1601441
ne = 4659888
ndof = 5
nvar = 5
gamma = 1.4
cv = 1/(gamma*(gamma-1))
x = np.empty(nn, dtype='float64')
y = np.empty(nn, dtype='float64')
z = np.empty(nn, dtype='float64')
ien = np.empty((ne,8), dtype='int32')
d = np.empty((nn,ndof), dtype='float64')
d_trans = np.empty((nn,nvar), dtype='float64')

#%% Filenames
dataFile = 'data.0012'
outFile = 'Re10k_bp1p6_8000float.szplt'
xyzFile = 'xyz.bin'
ienFile = 'ien.bin'

#%% Read xyz and ien data
xyzfile = open(xyzFile,'rb')
ienfile = open(ienFile,'rb')

for i in range(nn):
    xyz = xyzfile.read(24)
    x[i], y[i], z[i] = struct.unpack('ddd', xyz)

for i in range(ne):
    for j in range(8):
        iendata = ienfile.read(4)
        ien[i,j] = struct.unpack('i', iendata)[0]
xyzfile.close()
ienfile.close()

#%% Read data
datafile = open(dataFile,'rb')
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
#    d_trans[i,5] = (d_trans[i,1]**2 + d_trans[i,2]**2 \
#                    + d_trans[i,3]**2)**0.5/d_trans[i,4]**0.5
#    d_trans[i,6] = d_trans[i,0]*d_trans[i,4]/gamma
datafile.close()

#%% Write data
varNames = ['x','y','z','rho','uvel','vvel','wvel','temp']
file_handle = open_file(outFile, "Y-duct", varNames)
zone = create_FE_zone(file_handle, "Zone", nn, ne, FEBRICK, None, \
                               [FD_FLOAT for i in range(len(varNames))])
zone_write_float_values(file_handle, zone, 1, x)
zone_write_float_values(file_handle, zone, 2, y)
zone_write_float_values(file_handle, zone, 3, z)
for i in range(4, len(varNames) + 1):
    zone_write_float_values(file_handle, zone, i, d_trans[:, i - 4])
zone_write_nodemap(file_handle, zone, ien.reshape(-1))
close_file(file_handle)


