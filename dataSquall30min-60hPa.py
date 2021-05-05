import numpy as np
from netCDF4 import Dataset

data = Dataset(
    'NetCDF/rceiso_squall2_U10_H1000_m60hPad_rceiso01_3D_dt30min_jour40.nc')
data2 = Dataset(
    'NetCDF/rceiso_squall2_U10_H1000_m60hPad_rceiso01_3D_dt30min_jour41.nc')
data3 = Dataset(
    'NetCDF/rceiso_squall2_U10_H1000_m60hPad_rceiso01_3D_dt30min_jour42.nc')
data4 = Dataset(
    'NetCDF/rceiso_squall2_U10_H1000_m60hPad_rceiso01_3D_dt30min_jour43.nc')
data5 = Dataset(
    'NetCDF/rceiso_squall2_U10_H1000_m60hPad_rceiso01_3D_dt30min_jour44.nc')
data6 = Dataset(
    'NetCDF/rceiso_squall2_U10_H1000_m60hPad_rceiso01_3D_dt30min_jour45.nc')
data7 = Dataset(
    'NetCDF/rceiso_squall2_U10_H1000_m60hPad_rceiso01_3D_dt30min_jour46.nc')
data8 = Dataset(
    'NetCDF/rceiso_squall2_U10_H1000_m60hPad_rceiso01_3D_dt30min_jour47.nc')
data9 = Dataset(
    'NetCDF/rceiso_squall2_U10_H1000_m60hPad_rceiso01_3D_dt30min_jour48.nc')
#data10 = Dataset('NetCDF/rceiso_squall2_U10_H1000_m60hPad_rceiso01_3D_dt30min_jour49.nc')

print("data:")
for d in data.variables:
    print(d, data.variables[d].dimensions, data.variables[d].shape)

print('PROGRESS')

print(0/13*100, '%')

qc = data.variables['QC'][:]
qi = data.variables['QI'][:]
qsat = data.variables['QSAT'][:]
qv = data.variables['QV'][:]
x = np.asarray(data.variables['x'][:])
y = np.asarray(data.variables['y'][:])
z = np.asarray(data.variables['z'][:])
t = data.variables['time'][:]

print(1/13*100, '%')

qc2 = data2.variables['QC'][:]
qi2 = data2.variables['QI'][:]
qsat2 = data2.variables['QSAT'][:]
qv2 = data2.variables['QV'][:]
# x2 = np.asarray(data2.variables['x'][:])
# y2 = np.asarray(data2.variables['y'][:])
# z2 = np.asarray(data2.variables['z'][:])
t2 = data2.variables['time'][:]

print(2/13*100, '%')

qc3 = data3.variables['QC'][:]
qi3 = data3.variables['QI'][:]
qsat3 = data3.variables['QSAT'][:]
qv3 = data3.variables['QV'][:]
# x3 = np.asarray(data3.variables['x'][:])
# y3 = np.asarray(data3.variables['y'][:])
# z3 = np.asarray(data3.variables['z'][:])
t3 = data3.variables['time'][:]

print(3/13*100, '%')

qc4 = data4.variables['QC'][:]
qi4 = data4.variables['QI'][:]
qsat4 = data4.variables['QSAT'][:]
qv4 = data4.variables['QV'][:]
# x4 = np.asarray(data4.variables['x'][:])
# y4 = np.asarray(data4.variables['y'][:])
# z4 = np.asarray(data4.variables['z'][:])
t4 = data4.variables['time'][:]

print(4/13*100, '%')

qc5 = data5.variables['QC'][:]
qi5 = data5.variables['QI'][:]
qsat5 = data5.variables['QSAT'][:]
qv5 = data5.variables['QV'][:]
# x5 = np.asarray(data5.variables['x'][:])
# y5 = np.asarray(data5.variables['y'][:])
# z5 = np.asarray(data5.variables['z'][:])
t5 = data5.variables['time'][:]

print(5/13*100, '%')

qc6 = data6.variables['QC'][:]
qi6 = data6.variables['QI'][:]
qsat6 = data6.variables['QSAT'][:]
qv6 = data6.variables['QV'][:]
# x6 = np.asarray(data6.variables['x'][:])
# y6 = np.asarray(data6.variables['y'][:])
# z6 = np.asarray(data6.variables['z'][:])
t6 = data6.variables['time'][:]

print(6/13*100, '%')

qc7 = data7.variables['QC'][:]
qi7 = data7.variables['QI'][:]
qsat7 = data7.variables['QSAT'][:]
qv7 = data7.variables['QV'][:]
# x7 = np.asarray(data7.variables['x'][:])
# y7 = np.asarray(data7.variables['y'][:])
# z7 = np.asarray(data7.variables['z'][:])
t7 = data7.variables['time'][:]

print(8/13*100, '%')

qc8 = data8.variables['QC'][:]
qi8 = data8.variables['QI'][:]
qsat8 = data8.variables['QSAT'][:]
qv8 = data8.variables['QV'][:]
# x8 = np.asarray(data8.variables['x'][:])
# y8 = np.asarray(data8.variables['y'][:])
# z8 = np.asarray(data8.variables['z'][:])
t8 = data8.variables['time'][:]

print(9/13*100, '%')

qc9 = data9.variables['QC'][:]
qi9 = data9.variables['QI'][:]
qsat9 = data9.variables['QSAT'][:]
qv9 = data9.variables['QV'][:]
# x9 = np.asarray(data9.variables['x'][:])
# y9 = np.asarray(data9.variables['y'][:])
# z9 = np.asarray(data9.variables['z'][:])
t9 = data9.variables['time'][:]

print(10/13*100, '%')

# qc10 = data10.variables['QC'][:]
# qi10 = data10.variables['QI'][:]
# qsat10 = data10.variables['QSAT'][:]
# qv10 = data10.variables['QV'][:]
# x10 = np.asarray(data10.variables['x'][:])
# y10 = np.asarray(data10.variables['y'][:])
# z10 = np.asarray(data10.variables['z'][:])
# t10 = data10.variables['time'][:]

nsamples = 9

t = t-23
for i in np.arange(2, nsamples+1):
    locals()['t'+str(i)] = locals()['t'+str(i)] - 23

w = data.variables['W'][:]
w2 = data2.variables['W'][:]
w3 = data3.variables['W'][:]
w4 = data4.variables['W'][:]
w5 = data5.variables['W'][:]
w6 = data6.variables['W'][:]
w7 = data7.variables['W'][:]
w8 = data8.variables['W'][:]
w9 = data9.variables['W'][:]
# w10 = data10.variables['W'][:]

print(11/13*100, '%')

p = data.variables['p'][:]
p2 = data2.variables['p'][:]
p3 = data3.variables['p'][:]
p4 = data4.variables['p'][:]
p5 = data5.variables['p'][:]
p6 = data6.variables['p'][:]
p7 = data7.variables['p'][:]
p8 = data8.variables['p'][:]
p9 = data9.variables['p'][:]
# p10 = data10.variables['p'][:]

print(12/13*100, '%')

plist = [p, p2, p3, p4, p5, p6, p7, p8, p9]#, p10]
nzlist = []
nz = len(p)
nzlist.append(nz)
for i in np.arange(2, nsamples+1):
    locals()["nz"+str(i)] = len(plist[i-1])
    nzlist.append(locals()["nz"+str(i)])

omegamaxlist = []
omegamax = -60
omegamaxlist.append(omegamax)
# for i in np.arange(2,nsamples+1):
#     locals()["omegamax" + str(i)] = 0
#     omegamaxlist.append(locals()["omegamax" + str(i)])

for i in np.arange(2, nsamples+1):
    locals()["omegamax" + str(i)] = -60
    omegamaxlist.append(locals()["omegamax" + str(i)])

print(13/13*100, '%')
