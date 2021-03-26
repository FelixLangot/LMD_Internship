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
    RHP = np.empty((len(T), len(X), len(Y)))
    RHP[:, :, :] = np.nan
    for time in np.arange(len(T)):
        for i in np.arange(len(X)):
            for j in np.arange(len(Y)):
                RHP[time, i, j] = QSAT[time, iZclouds[time, i, j], j, i] / \
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
