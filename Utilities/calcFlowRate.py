import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import struct
import timeit

#%%
nn = 1601441
ne = 4659888

#%% Read xyz and ien data
x = np.empty(nn, dtype='float64')
y = np.empty(nn, dtype='float64')
z = np.empty(nn, dtype='float64')
ien = np.empty((ne,8), dtype='int32')

xyzfile = open('/home/debartha/MTProject/1p6m_mesh/xyz.bin','rb')
ienfile = open('/home/debartha/MTProject/1p6m_mesh/ien.bin','rb')

for i in range(nn):
    xyz = xyzfile.read(24)
    x[i], y[i], z[i] = struct.unpack('ddd', xyz)

for i in range(ne):
    for j in range(8):
        iendata = ienfile.read(4)
        ien[i,j] = struct.unpack('i', iendata)[0]
xyzfile.close()
ienfile.close()

#%% inlet
bnd1nn = 3086
bnd1ne = 5990
bnd1nodes = np.empty(bnd1nn, dtype='int64')
bnd1ien = np.empty((bnd1ne,3), dtype='int64')
bnd1file = open('/home/debartha/MTProject/1p6m_mesh/inletUpper.bin','rb')
for i in range(bnd1nn):
    inode = bnd1file.read(8)
    bnd1nodes[i] = struct.unpack('q',inode)[0]
for i in range(bnd1ne):
    for j in range(3):
        iendata = bnd1file.read(8)
        bnd1ien[i][j] = struct.unpack('q',iendata)[0]
bnd1file.close()

#%% merger
bnd2nn = 3973
bnd2ne = 7784
bnd2nodes = np.empty(bnd2nn, dtype='int64')
bnd2ien = np.empty((bnd2ne,3), dtype='int64')
bnd2file = open('/home/debartha/MTProject/1p6m_mesh/merger.bin','rb')
for i in range(bnd2nn):
    inode = bnd2file.read(8)
    bnd2nodes[i] = struct.unpack('q',inode)[0]
for i in range(bnd2ne):
    for j in range(3):
        iendata = bnd2file.read(8)
        bnd2ien[i][j] = struct.unpack('q',iendata)[0]
bnd2file.close()

#%% outlet
bnd3nn = 195
bnd3ne = 344
bnd3nodes = np.empty(bnd3nn, dtype='int64')
bnd3ien = np.empty((bnd3ne,3), dtype='int64')
bnd3file = open('/home/debartha/MTProject/1p6m_mesh/outlet.bin','rb')
for i in range(bnd3nn):
    inode = bnd3file.read(8)
    bnd3nodes[i] = struct.unpack('q',inode)[0]
for i in range(bnd3ne):
    for j in range(3):
        iendata = bnd3file.read(8)
        bnd3ien[i][j] = struct.unpack('q',iendata)[0]
bnd3file.close()

#%%
ndof = 5
nvar = 7
gamma = 1.4
cv = 1/(gamma*(gamma-1))
d = np.empty((nn,ndof), dtype='float64')
d_trans = np.empty((nn,nvar), dtype='float64')
datafile = open('/home/debartha/MTProject/1p6m_mesh/data.0012','rb')
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
datafile.close()
'''
#%%
bnd1ienl = np.empty((bnd1ne,3), dtype='int64')
for i in range(bnd1ne):
    for j in range(3):
        bnd1ienl[i,j] = np.where(bnd1nodes == bnd1ien[i,j])[0]
        
#%%
triang1 = tri.Triangulation(z[bnd1nodes], y[bnd1nodes], bnd1ienl)
levelVal = np.linspace(0,1.6,17)
plt.figure(1)
plt.tricontourf(triang1, d_trans[bnd1nodes,5], levels=levelVal, cmap='jet')
plt.colorbar()
plt.axis('equal')
plt.show()

#%%
bnd2ienl = np.empty((bnd2ne,3), dtype='int64')
for i in range(bnd2ne):
    for j in range(3):
        bnd2ienl[i,j] = np.where(bnd2nodes == bnd2ien[i,j])[0]
        
#%%
triang2 = tri.Triangulation(z[bnd2nodes], y[bnd2nodes], bnd2ienl)
levelVal = np.linspace(0,1.6,17)
plt.figure(2)
plt.tricontourf(triang2, d_trans[bnd2nodes,5], levels=levelVal, cmap='jet')
plt.colorbar()
plt.axis('equal')
plt.show()

#%%
bnd3ienl = np.empty((bnd3ne,3), dtype='int64')
for i in range(bnd3ne):
    for j in range(3):
        bnd3ienl[i,j] = np.where(bnd3nodes == bnd3ien[i,j])[0]
        
#%%
triang3 = tri.Triangulation(z[bnd3nodes], y[bnd3nodes], bnd3ienl)
levelVal = np.linspace(0,1.6,17)
plt.figure(3)
plt.tricontourf(triang3, d_trans[bnd3nodes,5], levels=levelVal, cmap='jet')
plt.colorbar()
plt.axis('equal')
plt.show()
'''
#%%
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

#%%
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

#%%
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
        
#%%
print(flowRate1*2)
print(flowRate2)
print(flowRate3)