import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import struct

nn = 1601441
bnd1nn = 3086
bnd1ne = 5990
bnd2nn = 3973
bnd2ne = 7784
bnd3nn = 195
bnd3ne = 344
ndof = 5
nvar = 2    # rho and u are only required
gamma = 1.4
cv = 1/(gamma*(gamma-1))

dataFile = 'data.0012'
xyzFile = 'xyz.bin'
bnd1File = 'inletUpper.bin'
bnd2File = 'merger.bin'
bnd3File = 'outlet.bin'

x = np.empty(nn, dtype='float64')
y = np.empty(nn, dtype='float64')
z = np.empty(nn, dtype='float64')

xyzfile = open(xyzFile,'rb')

for i in range(nn):
    xyz = xyzfile.read(24)
    x[i], y[i], z[i] = struct.unpack('ddd', xyz)

xyzfile.close()


bnd1nodes = np.empty(bnd1nn, dtype='int64')
bnd1ien = np.empty((bnd1ne,3), dtype='int64')
bnd1file = open(bnd1File,'rb')
for i in range(bnd1nn):
    inode = bnd1file.read(8)
    bnd1nodes[i] = struct.unpack('q',inode)[0]
for i in range(bnd1ne):
    for j in range(3):
        iendata = bnd1file.read(8)
        bnd1ien[i][j] = struct.unpack('q',iendata)[0]
bnd1file.close()

bnd2nodes = np.empty(bnd2nn, dtype='int64')
bnd2ien = np.empty((bnd2ne,3), dtype='int64')
bnd2file = open(bnd2File,'rb')
for i in range(bnd2nn):
    inode = bnd2file.read(8)
    bnd2nodes[i] = struct.unpack('q',inode)[0]
for i in range(bnd2ne):
    for j in range(3):
        iendata = bnd2file.read(8)
        bnd2ien[i][j] = struct.unpack('q',iendata)[0]
bnd2file.close()

bnd3nodes = np.empty(bnd3nn, dtype='int64')
bnd3ien = np.empty((bnd3ne,3), dtype='int64')
bnd3file = open(bnd3File,'rb')
for i in range(bnd3nn):
    inode = bnd3file.read(8)
    bnd3nodes[i] = struct.unpack('q',inode)[0]
for i in range(bnd3ne):
    for j in range(3):
        iendata = bnd3file.read(8)
        bnd3ien[i][j] = struct.unpack('q',iendata)[0]
bnd3file.close()

d = np.empty((nn,ndof), dtype='float64')
d_trans = np.empty((nn,nvar), dtype='float64')
datafile = open(dataFile,'rb')
for i in range(nn):
    for j in range(ndof):
        data = datafile.read(8)
        d[i,j] = struct.unpack('d', data)[0]
    d_trans[i,0] = d[i,0]
    d_trans[i,1] = d[i,1]/d[i,0]
datafile.close()

rg = np.array([2/3, 1/6, 1/6], dtype='float64')
sg = np.array([1/6, 2/3, 1/6], dtype='float64')
wg = np.array([1/3, 1/3, 1/3], dtype='float64')
ze = np.empty(3, dtype='float64')
ye = np.empty(3, dtype='float64')
ue = np.empty(3, dtype='float64')
rhoe = np.empty(3, dtype='float64')
sh = np.empty(3, dtype='float64')
flowRate1 = 0.0
c = 0.5
for ie, elem in enumerate(bnd1ien):
    ze = z[elem].copy()
    ye = y[elem].copy()
    ue = d_trans[elem,1].copy()
    rhoe = d_trans[elem,0].copy()
    j = np.abs((ze[0]-ze[2])*(ye[1]-ye[2]) - (ze[1]-ze[2])*(ye[0]-ye[2]))
    for iquad in range(3):
        sh[0] = rg[iquad]
        sh[1] = sg[iquad]
        sh[2] = 1 - rg[iquad] - sg[iquad]
        ug = sh.dot(ue)
        rhog = sh.dot(rhoe)
        flowRate1 += c*ug*rhog*j*wg[iquad]

rg = np.array([2/3, 1/6, 1/6], dtype='float64')
sg = np.array([1/6, 2/3, 1/6], dtype='float64')
wg = np.array([1/3, 1/3, 1/3], dtype='float64')
ze = np.empty(3, dtype='float64')
ye = np.empty(3, dtype='float64')
ue = np.empty(3, dtype='float64')
rhoe = np.empty(3, dtype='float64')
sh = np.empty(3, dtype='float64')
flowRate2 = 0.0
c = 0.5
for ie, elem in enumerate(bnd2ien):
    ze = z[elem].copy()
    ye = y[elem].copy()
    ue = d_trans[elem,1].copy()
    rhoe = d_trans[elem,0].copy()
    j = np.abs((ze[0]-ze[2])*(ye[1]-ye[2]) - (ze[1]-ze[2])*(ye[0]-ye[2]))
    for iquad in range(3):
        sh[0] = rg[iquad]
        sh[1] = sg[iquad]
        sh[2] = 1 - rg[iquad] - sg[iquad]
        ug = sh.dot(ue)
        rhog = sh.dot(rhoe)
        flowRate2 += c*ug*rhog*j*wg[iquad]

rg = np.array([2/3, 1/6, 1/6], dtype='float64')
sg = np.array([1/6, 2/3, 1/6], dtype='float64')
wg = np.array([1/3, 1/3, 1/3], dtype='float64')
ze = np.empty(3, dtype='float64')
ye = np.empty(3, dtype='float64')
ue = np.empty(3, dtype='float64')
rhoe = np.empty(3, dtype='float64')
sh = np.empty(3, dtype='float64')
flowRate3 = 0.0
c = 0.5
for ie, elem in enumerate(bnd3ien):
    ze = z[elem].copy()
    ye = y[elem].copy()
    ue = d_trans[elem,1].copy()
    rhoe = d_trans[elem,0].copy()
    j = np.abs((ze[0]-ze[2])*(ye[1]-ye[2]) - (ze[1]-ze[2])*(ye[0]-ye[2]))
    for iquad in range(3):
        sh[0] = rg[iquad]
        sh[1] = sg[iquad]
        sh[2] = 1 - rg[iquad] - sg[iquad]
        ug = sh.dot(ue)
        rhog = sh.dot(rhoe)
        flowRate3 += c*ug*rhog*j*wg[iquad]

print('mdot at x = 0.0572: ' + str(flowRate1*2))
print('mdot at x = 1.0: ' + str(flowRate2))
print('mdot at x = 2.8: ' + str(flowRate3))