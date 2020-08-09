import numpy as np
import matplotlib.pyplot as plt
import struct
import timeit

#%%
def area(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1)) -np.dot(y,np.roll(x,1)))

#%%
nn = 1601441
ne = 4659888

#%% ASCII Read data
#time0 = timeit.default_timer()
file = open('/home/debartha/MTProject/1p6m_mesh/yduct_change_main.dat', 'r')

x, y, z = np.loadtxt(file, dtype='float64', unpack=True, max_rows=nn)
ien = np.loadtxt(file, dtype='int32', max_rows=ne)
file.close()
#time1 = timeit.default_timer()
#print('ASCII Read time: '+ str(time1-time0))
 
#%% Binary Read data
#time0 = timeit.default_timer()
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
#time1 = timeit.default_timer()
#print('Binary Read time: '+ str(time1-time0))

#%%
x1 = 0.0572
index1 = (np.abs(x-x1)).argmin()
nodelist1 = []
for i in range(nn):
    if x[i] == x[index1] and abs(y[i]-0.25) < 0.15 and abs(z[i]) < 0.15:
        r = ((y[i] - 0.20309)**2 + z[i]**2)**0.5
        if r < 0.1166:
            nodelist1.append(i)
inode1 = np.array(nodelist1)
#plt.plot(z[inode1], y[inode1], '.')

#%%
x2 = 1.0
index2 = (np.abs(x-x2)).argmin()
inode2 = np.where(x == x[index2])[0]
#plt.plot(z[inode2], y[inode2], '.')

#%%
x3 = 2.8
index3 = (np.abs(x-x3)).argmin()
inode3 = np.where(x == x[index3])[0]
#plt.plot(z[inode3], y[inode3], '.')

#%%
elemList1 = []
for inode in inode1:
    adjElems = np.where(ien == inode+1)[0]
    for ielem in adjElems:
        elem1 = set([inode])
        for inen in ien[ielem,:]:
            if inen-1 in inode1:
                elem1.add(inen-1)
        if elem1 not in elemList1 and len(elem1) >= 3:
            elemList1.append(elem1)

for ipoint, points in enumerate(elemList1):
    ptlist = list(points)
    elemList1[ipoint] = ptlist
    
#%%
elemList2 = []
for inode in inode2:
    adjElems = np.where(ien == inode+1)[0]
    for ielem in adjElems:
        elem2 = set([inode])
        for inen in ien[ielem,:]:
            if inen-1 in inode2:
                elem2.add(inen-1)
        if elem2 not in elemList2 and len(elem2) >= 3:
            elemList2.append(elem2)

for ipoint, points in enumerate(elemList2):
    ptlist = list(points)
    elemList2[ipoint] = ptlist

#%%
elemList3 = []
for inode in inode3:
    adjElems = np.where(ien == inode+1)[0]
    for ielem in adjElems:
        elem3 = set([inode])
        for inen in ien[ielem,:]:
            if inen-1 in inode3:
                elem3.add(inen-1)
        if elem3 not in elemList3 and len(elem3) >= 3:
            elemList3.append(elem3)
            
for ipoint, points in enumerate(elemList3):
    ptlist = list(points)
    elemList3[ipoint] = ptlist

#%%
for ielem, elem in enumerate(elemList1):
    if len(elem) == 4:
        nodeOrder = [[elem[0], elem[1], elem[2], elem[3]],
                     [elem[1], elem[0], elem[2], elem[3]],
                     [elem[2], elem[1], elem[0], elem[3]],
                     [elem[3], elem[1], elem[2], elem[0]]]
        maxArea = 0.0
        imaxArea = 0
        for i in range(4):
            elemArea = area(z[nodeOrder[i]], y[nodeOrder[i]])
            if elemArea > maxArea:
                maxArea = elemArea
                imaxArea = i
        elemList1[ielem] = nodeOrder[imaxArea].copy()           
        
#%%
triList1 = []
for elem in elemList1:
    if len(elem) == 4:
        triList1.append([elem[0], elem[1], elem[2]])
        triList1.append([elem[2], elem[3], elem[0]])
    else:
        triList1.append(elem)

#%%
for ielem, elem in enumerate(elemList2):
    if len(elem) == 4:
        nodeOrder = [[elem[0], elem[1], elem[2], elem[3]],
                     [elem[1], elem[0], elem[2], elem[3]],
                     [elem[2], elem[1], elem[0], elem[3]],
                     [elem[3], elem[1], elem[2], elem[0]]]
        maxArea = 0.0
        imaxArea = 0
        for i in range(4):
            elemArea = area(z[nodeOrder[i]], y[nodeOrder[i]])
            if elemArea > maxArea:
                maxArea = elemArea
                imaxArea = i
        elemList2[ielem] = nodeOrder[imaxArea].copy()           
        
#%%
triList2 = []
for elem in elemList2:
    if len(elem) == 4:
        triList2.append([elem[0], elem[1], elem[2]])
        triList2.append([elem[2], elem[3], elem[0]])
    else:
        triList2.append(elem)

#%%
triList3 = []
for elem in elemList3:
    if len(elem) == 4:
        triList3.append([elem[0], elem[1], elem[2]])
        triList3.append([elem[2], elem[3], elem[0]])
    else:
        triList3.append(elem)
        
#%%
quad1 = 0
tri1 = 0
for elem in elemList1:
    if len(elem) == 4:
        quad1 += 1
    else:
        tri1 += 1

#%%
quad2 = 0
tri2 = 0
for elem in elemList2:
    if len(elem) == 4:
        quad2 += 1
    else:
        tri2 += 1
        
#%%
for ptlist in triList2:
    if len(ptlist) == 4:
        zc = np.array([z[ptlist[0]], z[ptlist[1]], z[ptlist[2]], z[ptlist[3]], z[ptlist[0]]])
        yc = np.array([y[ptlist[0]], y[ptlist[1]], y[ptlist[2]], y[ptlist[3]], y[ptlist[0]]])
        plt.plot(zc, yc, '-k')
    else:
        zc = np.array([z[ptlist[0]], z[ptlist[1]], z[ptlist[2]], z[ptlist[0]]])
        yc = np.array([y[ptlist[0]], y[ptlist[1]], y[ptlist[2]], y[ptlist[0]]])
        plt.plot(zc, yc, '-k')
        
#%%
xfile = open('/home/debartha/MTProject/1p6m_mesh/xyz.bin','wb')
for i in range(nn):
    xfile.write(x[i])
    xfile.write(y[i])
    xfile.write(z[i])
xfile.close()

#%%
xtest = np.empty(nn, dtype='float64')
ytest = np.empty(nn, dtype='float64')
ztest = np.empty(nn, dtype='float64')
xyzfile = open('/home/debartha/MTProject/1p6m_mesh/xyz.bin','rb')
for i in range(nn):
    xyz = xyzfile.read(24)
    xtest[i], ytest[i], ztest[i] = struct.unpack('ddd', xyz)
xyzfile.close()

#%%
ienfile = open('/home/debartha/MTProject/1p6m_mesh/ien.bin','wb')
for i in range(ne):
    for j in range(8):
        ienfile.write(ien[i,j])
ienfile.close()

#%%
ientest = np.empty((ne,8), dtype='int32')
ienfile = open('/home/debartha/MTProject/1p6m_mesh/ien.bin','rb')
for i in range(ne):
    for j in range(8):
        iendata = ienfile.read(4)
        ientest[i,j] = struct.unpack('i', iendata)[0]
ienfile.close()

#%%
bnd1file = open('/home/debartha/MTProject/1p6m_mesh/inletUpper.bin','wb')
bnd1nn = len(inode1)
bnd1ne = len(triList1)
for i in range(bnd1nn):
    bnd1file.write(inode1[i])
for i in range(bnd1ne):
    for j in range(3):
        bnd1file.write(triList1[i][j])
bnd1file.close()

#%%
bnd2file = open('/home/debartha/MTProject/1p6m_mesh/merger.bin','wb')
bnd2nn = len(inode2)
bnd2ne = len(triList2)
for i in range(bnd2nn):
    bnd2file.write(inode2[i])
for i in range(bnd2ne):
    for j in range(3):
        bnd2file.write(triList2[i][j])
bnd2file.close()

#%%
bnd3file = open('/home/debartha/MTProject/1p6m_mesh/outlet.bin','wb')
bnd3nn = len(inode3)
bnd3ne = len(triList3)
for i in range(bnd3nn):
    bnd3file.write(inode3[i])
for i in range(bnd3ne):
    for j in range(3):
        bnd3file.write(triList3[i][j])
bnd3file.close()

#%%
