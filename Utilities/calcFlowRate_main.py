import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import struct

def progressbar(progress, barLength=15):
    percent = float(progress) * 100 / barLength
    done = '#' * progress
    rem = '.' * (barLength - len(done))

    print('[%s%s] %d %%' % (done, rem, percent), end='\r')

# Input

start = 1
end = 15
filename = 'ssa4_bp1p6_flowRate.dat'

# Code main
print('Setting up environment for mass flow rate calculations: ')
progress = 0
bar = 100
progressbar(progress, barLength=bar)
nn = 1601441
bnd1nn = 3086
bnd1ne = 5990
bnd2nn = 3973
bnd2ne = 7784
bnd3nn = 195
bnd3ne = 344
bnd4nn = 3086
bnd4ne = 5990
tot = nn + bnd1nn + bnd1ne + bnd2nn + bnd2ne + bnd3nn + bnd3ne + bnd4nn + bnd4ne
curr = 0
target = tot / bar 
ndof = 5
nvar = 2    # rho and u are only required
gamma = 1.4
cv = 1/(gamma*(gamma-1))

dataFiles = ['data.' + str(i).zfill(2) for i in range(start, end+1)]
xyzFile = 'mxyz'
bnd1File = 'inletUpper.bin'
bnd2File = 'merger.bin'
bnd3File = 'outlet.bin'
bnd4File = 'inletLower.bin'

x = np.empty(nn, dtype='float64')
y = np.empty(nn, dtype='float64')
z = np.empty(nn, dtype='float64')

xyzfile = open(xyzFile,'rb')

for i in range(nn):
    xyz = xyzfile.read(24)
    x[i], y[i], z[i] = struct.unpack('ddd', xyz)
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar
        progressbar(progress, barLength=bar)
curr = nn
xyzfile.close()

bnd1nodes = np.empty(bnd1nn, dtype='int64')
bnd1ien = np.empty((bnd1ne,3), dtype='int64')
bnd1file = open(bnd1File,'rb')
for i in range(bnd1nn):
    inode = bnd1file.read(8)
    bnd1nodes[i] = struct.unpack('q',inode)[0]
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar 
        progressbar(progress, barLength=bar)
curr += bnd1nn
for i in range(bnd1ne):
    for j in range(3):
        iendata = bnd1file.read(8)
        bnd1ien[i][j] = struct.unpack('q',iendata)[0]
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar
        progressbar(progress, barLength=bar)
curr += bnd1ne
bnd1file.close()

bnd2nodes = np.empty(bnd2nn, dtype='int64')
bnd2ien = np.empty((bnd2ne,3), dtype='int64')
bnd2file = open(bnd2File,'rb')
for i in range(bnd2nn):
    inode = bnd2file.read(8)
    bnd2nodes[i] = struct.unpack('q',inode)[0]
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar
        progressbar(progress, barLength=bar)
curr += bnd2nn
for i in range(bnd2ne):
    for j in range(3):
        iendata = bnd2file.read(8)
        bnd2ien[i][j] = struct.unpack('q',iendata)[0]
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar
        progressbar(progress, barLength=bar)
curr += bnd2ne
bnd2file.close()

bnd3nodes = np.empty(bnd3nn, dtype='int64')
bnd3ien = np.empty((bnd3ne,3), dtype='int64')
bnd3file = open(bnd3File,'rb')
for i in range(bnd3nn):
    inode = bnd3file.read(8)
    bnd3nodes[i] = struct.unpack('q',inode)[0]
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar
        progressbar(progress, barLength=bar)
curr += bnd3nn
for i in range(bnd3ne):
    for j in range(3):
        iendata = bnd3file.read(8)
        bnd3ien[i][j] = struct.unpack('q',iendata)[0]
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar
        progressbar(progress, barLength=bar)
curr += bnd3ne
bnd3file.close()

bnd4nodes = np.empty(bnd4nn, dtype='int64')
bnd4ien = np.empty((bnd4ne,3), dtype='int64')
bnd4file = open(bnd4File,'rb')
for i in range(bnd4nn):
    inode = bnd4file.read(8)
    bnd4nodes[i] = struct.unpack('q',inode)[0]
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar 
        progressbar(progress, barLength=bar)
curr += bnd4nn
for i in range(bnd4ne):
    for j in range(3):
        iendata = bnd4file.read(8)
        bnd4ien[i][j] = struct.unpack('q',iendata)[0]
    if i >= target - curr:
        progress += 1
        target = (progress + 1) * tot / bar
        progressbar(progress, barLength=bar)
curr += bnd4ne
bnd4file.close()

d = np.empty((nn,ndof), dtype='float64')
d_trans = np.empty((nn,nvar), dtype='float64')

progressbar(progress + 1, barLength=bar)
print('\n')

for dataFile in dataFiles:
    print('Calculating mass flow rates for ' + dataFile + ': ')
    tot = nn + bnd1ne + bnd2ne + bnd3ne + bnd4ne
    curr = 0
    progress = 0
    target = tot / bar
    progressbar(progress, barLength=bar)
    datafile = open(dataFile,'rb')
    for i in range(nn):
        for j in range(ndof):
            data = datafile.read(8)
            d[i,j] = struct.unpack('d', data)[0]
        d_trans[i,0] = d[i,0]
        d_trans[i,1] = d[i,1]/d[i,0]
        if i >= target - curr:
            progress += 1
            target = (progress + 1) * tot / bar
            progressbar(progress, barLength=bar)
    curr = nn
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
        if ie >= target - curr:
            progress += 1
            target = (progress + 1) * tot / bar
            progressbar(progress, barLength=bar)
        for iquad in range(3):
            sh[0] = rg[iquad]
            sh[1] = sg[iquad]
            sh[2] = 1 - rg[iquad] - sg[iquad]
            ug = sh.dot(ue)
            rhog = sh.dot(rhoe)
            flowRate1 += c*ug*rhog*j*wg[iquad]
    curr += bnd1ne

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
        if ie >= target - curr:
            progress += 1
            target = (progress + 1) * tot / bar
            progressbar(progress, barLength=bar)
        for iquad in range(3):
            sh[0] = rg[iquad]
            sh[1] = sg[iquad]
            sh[2] = 1 - rg[iquad] - sg[iquad]
            ug = sh.dot(ue)
            rhog = sh.dot(rhoe)
            flowRate2 += c*ug*rhog*j*wg[iquad]
    curr += bnd2ne

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
        if ie >= target - curr:
            progress += 1
            target = (progress + 1) * tot / bar
            progressbar(progress, barLength=bar)
        for iquad in range(3):
            sh[0] = rg[iquad]
            sh[1] = sg[iquad]
            sh[2] = 1 - rg[iquad] - sg[iquad]
            ug = sh.dot(ue)
            rhog = sh.dot(rhoe)
            flowRate3 += c*ug*rhog*j*wg[iquad]
    curr += bnd3ne
    
    ze = np.empty(3, dtype='float64')
    ye = np.empty(3, dtype='float64')
    ue = np.empty(3, dtype='float64')
    rhoe = np.empty(3, dtype='float64')
    sh = np.empty(3, dtype='float64')
    flowRate4 = 0.0
    c = 0.5
    for ie, elem in enumerate(bnd4ien):
        ze = z[elem].copy()
        ye = y[elem].copy()
        ue = d_trans[elem,1].copy()
        rhoe = d_trans[elem,0].copy()
        j = np.abs((ze[0]-ze[2])*(ye[1]-ye[2]) - (ze[1]-ze[2])*(ye[0]-ye[2]))
        if ie >= target - curr:
            progress += 1
            target = (progress + 1) * tot / bar
            progressbar(progress, barLength=bar)
        for iquad in range(3):
            sh[0] = rg[iquad]
            sh[1] = sg[iquad]
            sh[2] = 1 - rg[iquad] - sg[iquad]
            ug = sh.dot(ue)
            rhog = sh.dot(rhoe)
            flowRate4 += c*ug*rhog*j*wg[iquad]
    curr += bnd4ne
    
    #print('mdot at x = 0.0572: ' + str(flowRate1*2))
    #print('mdot at x = 1.0: ' + str(flowRate2))
    #print('mdot at x = 2.8: ' + str(flowRate3))

    writefile = open(filename, 'a')
    writefile.write(str(flowRate1) + '\t' + str(flowRate4) + '\t' + str(flowRate2) + '\t' + str(flowRate3) + '\n')
    writefile.close()
    
    progressbar(progress+1, barLength=bar)
    print('\n')
