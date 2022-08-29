########################################################################
#                                                                      #
#     OMNI-STRAIN FAILURE ENVELOPE PROGRAM developed using the         #
#     3D INVARIANT-BASED FAILURE CRITERIA for FRP composites           #
#                                                                      #
########################################################################
#                                                                      #
#     Implemented by                                                   #
#     Giuseppe Corrado, PhD candidate at FEUP (Porto)                  #
#                     and ITN-EID Marie Skłodowska-Curie Fellow        #
#                     (https://optimacs.net/)                          #
#     email: gcorrado@fe.up.pt                                         #
#                                                                      #
#     Created on 03-Mar-2019                                           #
########################################################################
#
#   Main references:
#
#   - P.P. Camanho et al. /Int. Journal of Solids and Structures 55
#       (2015) 92–107
#   - S.W. Tsai, J.D.D. Melo /Composites Science and Techn. 100 (2014) 
#       237–243
#   - S.W. Tsai, J.D.D. Melo, S. Sihn, A. Arteiro, R. Rainsberger.
#       Composite Laminates: Theory and practice of analysis, design and 
#       automated layup. Stanford Aeronautics & Astronautics, 2017.
#
# -------------------------------------------------------------------- #

# -------------------------- Imports --------------------------------- #
from numpy import *
import numpy as np
#
from _includes.Inv3DFCFunctions import FailureMatrix
from _includes.Inv3DFCFunctions import FailureFiber

# -------------------------------------------------------------------- #

# ------------------------ local functions --------------------------- # 
def get_F_max(s11, s22, s33, s12, s13, s23, nu12, TOLER, DEG2RAD, 
              a1, a2, a32T, a3T, a32C, a3C, PHIC, XT, XC,
              YC, SL, ST, beta, G12):
    PHI_D = PSI_D = 0
    FM_MAX, FLAG_M, THETADEG_M = FailureMatrix(
        s22, s33, s12, s13, s23, TOLER, DEG2RAD, YC, SL, ST,
        a1, a2, a32T, a3T, a32C, a3C, PHIC)

    FF_MAX, FLAG_F, THETADEG_F, PHI_D, PSI_D = FailureFiber(
        s11, s22, s33, s12, s13, s23, nu12,
        TOLER, DEG2RAD, a1, a2, a32T, a3T, a32C, a3C, PHIC, XT, XC,
        YC, SL, ST, beta, G12)

    if max(FF_MAX,0) >= max(FM_MAX,0):
        F_MAX = FF_MAX
        FLAG = FLAG_F
        THETADEG = THETADEG_F
    else:
        F_MAX = FM_MAX
        FLAG = FLAG_M
        THETADEG = THETADEG_M
    if FLAG == 3 and abs(PHI_D)< TOLER:
        FLAG = 1

    return F_MAX, FLAG, THETADEG, PHI_D, PSI_D

def get_sij(SI, SII, e11, e22, e33, e12, e23, e13):
    if SI == 'e11':
        e11 = SII
    elif SI == 'e22':
        e22 = SII
    elif SI == 'e33':
        e33 = SII
    elif SI == 'e12':
        e12 = SII
    elif SI == 'e23':
        e23 = SII
    elif SI == 'e13':
        e13 = SII

    return e11, e22, e33, e12, e23, e13

def get_SXX_SYY(R, GAMM, SX, SY, T2, Kstif):
    e11 = e22 = e33 = e12 = e23 = e13 = 0
    SXX = R * cos(GAMM)
    SYY = R * sin(GAMM)

    e11, e22, e33, e12, e23, e13 = get_sij(SX, SXX, e11, e22, e33, e12, 
                                            e23, e13)
    e11, e22, e33, e12, e23, e13 = get_sij(SY, SYY, e11, e22, e33, e12, 
                                            e23, e13)

    Eps = np.array([e11, e22, e33, 2 * e23, 2 * e13, 2 * e12])
    EpsR = np.matmul(T2, Eps)
    SkR = np.matmul(Kstif, EpsR)
    s11 = SkR[0]
    s22 = SkR[1]
    s33 = SkR[2]
    s23 = SkR[3]
    s13 = SkR[4]
    s12 = SkR[5]
    return SXX, SYY, s11, s22, s33, s12, s23, s13

def get_delta_F(R, GAMM, SX, SY, TOLER, DEG2RAD, a1, a2, a32T, a3T, 
                            a32C, a3C, PHIC, XT, XC, YC, SL, ST, 
                            beta, G12, nu12, T2, Kstif):
                
    SXX, SYY, s11, s22, s33, s12, s23, s13 = get_SXX_SYY(R, GAMM, SX, 
                                                        SY, T2, Kstif)

    F_MAX, FLAG, THETADEG, PHI_D, PSI_D = get_F_max(s11, s22, s33, s12, 
                    s13, s23, nu12, TOLER, DEG2RAD, a1, a2, a32T, a3T,
                    a32C, a3C, PHIC, XT, XC, YC, SL, ST, beta, G12)
                                
    DELTA_F = 1.0 - F_MAX

    return DELTA_F, FLAG, THETADEG, PHI_D, PSI_D, SXX, SYY

# --------------------------------------------------------------------------- #