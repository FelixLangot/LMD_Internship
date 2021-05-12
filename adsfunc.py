import numpy as np
from numba import njit, jit

@njit
def Zparcel(Z):
    i = 0
    while Z[i] <= max(Z):
        if Z[i] >= 5000:
            Zparcel = Z[i]
            iZparcel = i
            break
        i += 1
    return iZparcel, Zparcel


@njit
def Zloc(Z,loc):
    i = 0
    while Z[i] <= max(Z):
        if Z[i] >= loc:
            Zloc = Z[i]
            iZloc = i
            break
        i += 1
    return iZloc, Zloc

@njit
def LSA(T, X, Y, Z, QI, QC, QSAT):
    zclouds = np.zeros((len(T), len(X), len(Y)))
    izclouds = np.zeros((len(T), len(X), len(Y)))
    qsatzclouds = np.zeros((len(T), len(X), len(Y)))
    thresh = 1e-6
    for time in np.arange(len(T)):
        for i in np.arange(len(X)):
            for j in np.arange(len(Y)):
                cond = QC[time, Zparcel(Z)[0]:, j, i] + QI[time, Zparcel(Z)[0]:, j, i]
                if np.max(cond) - thresh >= 0:
                    toto = np.where(cond-thresh >= 0)
                    (matrix,) = toto
                    zclouds[time, i, j] = Z[Zparcel(Z)[0]+np.min(matrix)]
                    izclouds[time, i, j] = int(Zparcel(Z)[0]+np.min(matrix))
                    qsatzclouds[time, i, j] = QSAT[time,Zparcel(Z)[0]+np.min(matrix), j, i]
                    if qsatzclouds[time, i, j] > 100:
                        qsatzclouds[time, i, j] = qsatzclouds[time, i, j]/1000
                else:
                    zclouds[time, i, j] = np.max(Z)
                    izclouds[time, i, j] = int(np.argmax(Z))
                    qsatzclouds[time, i, j] = QSAT[time, np.argmax(Z), j, i]
                    if qsatzclouds[time, i, j] > 100:
                        qsatzclouds[time, i, j] = qsatzclouds[time, i, j]/1000
    return zclouds, izclouds, qsatzclouds


@njit
def RELHUM(T, X, Y, Z, QV, QSAT):
    RH = np.zeros((len(T), len(X), len(Y)))
    for time in np.arange(len(T)):
        for i in np.arange(len(X)):
            for j in np.arange(len(Y)):
                if 100 < QV[time, Zparcel(Z)[0], j, i]/QSAT[time, Zparcel(Z)[0], j, i]:
                    RH[time, i, j] = (
                        QV[time, Zparcel(Z)[0], j, i]/QSAT[time, Zparcel(Z)[0], j, i])/1000
                else:
                    RH[time, i, j] = (
                        QV[time, Zparcel(Z)[0], j, i]/QSAT[time, Zparcel(Z)[0], j, i])

    return RH


@njit
def RELHUMP(T, X, Y, Z, QSAT, iZclouds):
    RHP = np.empty((len(X), len(Y)))
    RHP[:, :] = np.nan
    for time in np.arange(len(T)):
        for i in np.arange(len(X)):
            for j in np.arange(len(Y)):
                RHP[i, j] = QSAT[time, iZclouds[i, j], j, i] / \
                    QSAT[time, Zparcel(Z)[0], j, i]
    return RHP


@njit
def MWENV(T, X, Y, Z, W, QC, QI):
    thresh = 1e-6
    Wenv = np.zeros(W.shape)
    for time in np.arange(len(T)):
        for i in np.arange(len(X)):
            for j in np.arange(len(Y)):
                for k in np.arange(len(Z)):
                    cond2 = QC[time, k, j, i] + QI[time, k, j, i]
                    if cond2 < thresh:
                        Wenv[time, k, j, i] = W[time, k, j, i]
                    else:
                        Wenv[time, k, j, i] = np.nan
    return Wenv


from const_thermo import *

SST = 30
pp = 1000
pvzero = 100
H = 10
pmax = 500

def WLS(NZ, P, PMAX, OMEGAMAX):
    vmaxu = OMEGAMAX/86400*100.0
    pvmaxu = PMAX*100.0
    ppu = pp*100.0
    pvzerou = pvzero*100.0
    omeLS = np.zeros((NZ))
    wLS = np.zeros((NZ))
    for i in np.arange(NZ):
        plevu = P[i]*100
        xm = ppu - pvmaxu
        xtw = ppu - pvzerou
        Asb = -vmaxu*(3.*xm-2.*xtw)/(xm*(xtw-xm)**2)
        Bsb = -vmaxu*(xtw-2.*xm)/(xm*xm*(xtw-xm)**2)
        xsb = ppu-plevu
        omeu = (xtw-xsb)*(Asb*xsb+Bsb*xsb*xsb)

        if plevu < pvzerou:
            omeu = 0
        if omeu*OMEGAMAX < 0:
            omeu = 0

        ratio = (plevu+1e-10)/ppu
        zloc = -H*np.log(ratio)
        T = SST+degK-7.0*zloc/1000.0
        rho = plevu/Rd/T
        omeLS[i] = omeu*86400.0/100.0
        wLS[i] = -omeu/rho/g
    return wLS, omeLS


def WTOT(MWENV, WLS):
    wtot = np.zeros(len(WLS))
    for i in np.arange(len(wtot)):
        wtot[i] = MWENV[i]+WLS[i]
    return wtot


def ZTRAJ(TIME, W, NOBS, Z, Dt):
    Ztrajlist = []
    Ttot = len(TIME)
    nobs = NOBS
    for N in np.arange(Ttot-nobs, Ttot, 1):
        ztraj = np.zeros(N+1)
        ztraj[N] = Zparcel(Z)[1]
        maxtraj = np.max(ztraj)
        i = N-1
        for j in np.arange(Zparcel(Z)[0], np.argmax(Z)):
            while maxtraj < Z[j+1] and i >= 0:
                ztraj[i] = ztraj[i+1]-W[j]*Dt
                maxtraj = np.max(ztraj)
                i = i-1
        Ztrajlist.append(ztraj)
    return Ztrajlist


def ZTRAJBINHIST(Ztraj,Z, ZbinL, ZbinH):
    Ztrajbinlist = np.zeros(len(Ztraj), dtype=np.ndarray)
    Ztrajhistlist = np.zeros(len(Ztraj), dtype=np.ndarray)
    for i in np.arange(len(Ztrajbinlist)):
        Ztrajbinlist[i] = np.zeros(len(Ztraj[i]))
        Ztrajhistlist[i] = np.histogram(Ztraj[i], bins=Z[ZbinL:ZbinH])
        a = 0
        b = 0
        for j in np.arange(len(Ztrajhistlist[i][0])):
            b += Ztrajhistlist[i][0][j]
            Ztrajbinlist[i][a:b] = Ztrajhistlist[i][1][j]
            a = b
        Ztrajbinlist[i] = np.flip(Ztrajbinlist[i])
    return Ztrajbinlist, Ztrajhistlist


@njit
def LSADYN(X, Y, Z, QCC, QIC, Ztrajbin, Ztrajhist, Step):
    maxtrajloc = Zloc(Z, np.max(Ztrajbin))[0]
    zparceltraj = Zparcel(Z)[0]
    zclouds = np.zeros((len(X), len(Y)))
    izclouds = np.zeros((len(X), len(Y)))
    thresh = 1e-6
    for i in np.arange(len(X)):
        for j in np.arange(len(Y)):
            a = 0
            b = 0
            cond = np.zeros(len(Ztrajbin))
            for k, iz in zip(np.arange(len(Ztrajhist[0])), np.arange(zparceltraj, maxtrajloc+1)):
                b += Ztrajhist[0][k]
                for ttime in np.arange(a, b, Step):
                    cond[ttime] = QCC[ttime, iz, j, i] + QIC[ttime, iz, j, i]
                    a = b
            (matrix,) = np.where(cond > thresh)
            if matrix.shape == (0,):
                zclouds[i, j] = np.max(Ztrajbin)
                izclouds[i, j] = np.argmax(Ztrajbin)
            else:
                zclouds[i, j] = Ztrajbin[-np.min(matrix)]
                izclouds[i, j] = -np.min(matrix)
    return izclouds, zclouds

# def ZCLOUDS(X,Y,Z,QC,QI,Ztrajbin,Ztrajhist, Nsamples):
#     step = 1
#     Zcloudslist = np.zeros((10, 128*128), dtype=np.ndarray)
#     Zcloudslist[0] = LSADYN(X, Y, Z, QC, QI, Ztrajbin[0], Ztrajhist[0], step)[1].flatten()
#     iZclouds = LSADYN(X, Y, Z, QC, QI, Ztrajbin[0], Ztrajhist[0], step)[0].astype(int)
#     iZZclouds = np.zeros(iZclouds.shape)
#     for i in tqdmn(np.arange(len(X)), leave=False):
#         for j in np.arange(len(Y)):
#             (index,) = np.where(Ztrajbin[-1][iZclouds[i, j]] == Z)
#             iZZclouds[i, j] = index
#     iZZclouds = iZZclouds.astype(int)
#     for j in tqdmn(np.arange(2, Nsamples+1), leave=False):
#         Zcloudslist[j-1] = LSADYN(X, Y, Z, QC, QI, Ztrajbin[j-1],Ztrajhist[j-1], step)[1].flatten()
#         locals()['iZclouds_'+str(j)] = LSADYN(X, Y, Z, QC, QI,Ztrajbin[j-1], Ztrajhist[j-1], step)[0].astype(int)
#         locals()['iZZclouds_'+str(j)] = np.zeros(locals()['iZclouds_'+str(j)].shape)
#         for k in np.arange(len(X)):
#             for l in np.arange(len(Y)):
#                 (index,) = np.where(
#                     Ztrajbin[j-1][locals()['iZclouds_'+str(j)][k, l]] == z)
#                 locals()['iZZclouds_'+str(j)][k, l] = index
#         locals()['iZZclouds_'+str(j)
#                 ] = locals()['iZZclouds_'+str(j)].astype(int)
#     Zcloudslist = Zcloudslist.flatten()
#     iZZcloudslist = np.asaray([iZZclouds, iZZclouds_2, iZZclouds_3, iZZclouds_4, iZZclouds_5, iZZclouds_6, iZZclouds_7, iZZclouds_8, iZZclouds_9, iZZclouds_10])
#     return iZZclouds_list
