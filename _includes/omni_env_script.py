################################################################################
#                                                                              #
#      OMNI ENVELOPE PROGRAM applied to the                                    #
#      3D INVARIANT-BASED FAILURE CRITERIA FOR LAMINATED COMPOSITES            #
#                                                                              #
################################################################################
#                                                                              #
#   Implemented by                                                             #
#   Giuseppe Corrado, PhD candidate FEUP (Porto) - gcorrado@fe.up.pt           #
#                     and ITN-EID Marie Skłodowska-Curie Fellow                #
#                     (https://optimacs.net/)                                  #
#                                                                              #
#   Created on 03-Mar-2019                                                     #
################################################################################
#
#   Main references:
#
#   - P.P. Camanho et al. /Int. Journal of Solids and Structures 55(2015) 92–107
#   - S.W. Tsai, J.D.D. Melo /Composites Science and Techn. 100 (2014) 237–243
#
# ---------------------------------------------------------------------------- #
# BRIEF USERGUIDE:
# After the definition of the Material Properties in "INPUT", the number of 
# failure envelopes for the generation of omni strain FPF envelopes can be 
# controlled by the angular increment "deltaTheta".
# If you wish to increase the accuracy of the failure envelope, you can increase
# the number of failure loci, which is controlled by the parameter "ENNE", i.e.
# the number of failure points computed for each angular increment in strain space
# for each failure envelope.

# -------------------------------- Imports ------------------------------------ #
from numpy import *
import math
import time
import numpy as np
from matplotlib import pyplot as plt
import os
import xlsxwriter
#
import _include.material_data
from _include.Inv3DFCFunctions import FailureMatrix
from _include.Inv3DFCFunctions import FailureFiber
from _include.DegradedPropEvaluation import get_DegradedProp

# ------------------------------ User INPUT ----------------------------------- #

OmniLPF = 'NO'          #'YES' or 'NO': 
#                         Omni Last Ply Failure envelopes are generated 
#                         using degraded properties
ExcelOutput = 'YES'     #'YES' or 'NO'
PlotFigure = 'NO'       #'YES' or 'NO'
SAVEFIG = 'NO'          #'YES' or 'NO' 
#                         Note: This is a suboption of PlotFigure! If you want to save the figure,
#                               you need to set PlotFigure = 'YES'
#
## Material properties
#
# Choose the material from the database -
# materials = IM78552, custom ...
material = material_data.IM78552
mat_label = 'IM78552'

## Laminate design (just if failure envelope needs to be represented in STRESS SPACE)
if mat_label == 'AS43501_6':
    Layup = [90.0, +30.0, -30.0, -30.0, +30.0, 90.0]
elif mat_label == 'EglassLY556':
    Layup = [90.0, +45.0, -45.0, +0.0, +0.0, -45.0, +45.0, 90.0]
elif mat_label == 'EglassMY750':
    Layup = [+55.0, -55.0, -55.0, +55.0]
#Layup = [+35.0, -35.0, -35.0, +35.0]
else:
    # INPUT YOUR LAYUP BELOW
    #Layup = [+45.0, +0.0, -45.0, 90.0, 90.0, -45.0, +0.0,+45.0]
    Layup = [90.0,90.0,+0.0,+0.0,90.0,90.0]
    #Layup = [+0.0,+0.0,90.0,90.0,+0.0,+0.0]

#
N0p = Layup.count(0)
N30p = Layup.count(30)
N_30p = Layup.count(-30)
N45p = Layup.count(45)
N_45p = Layup.count(-45)
N60p = Layup.count(60)
N_60p = Layup.count(-60)
N90p = Layup.count(90)
N55p = Layup.count(55)
N_55p = Layup.count(-55)
N35p = Layup.count(35)
N_35p = Layup.count(-35)
Ntot = len(Layup)

## ELASTIC PROPERTIES
E1 = material['E1']         #Longitudinal Tension Young Modulus
E2 = material['E2']         #Transverse Young Modulus
nu12 = material['v12']      #In-plane Poisson's coefficient
nu23 = material['v23']      #Transverse Poisson's coefficient
G12 = material['G12']       #In-plane Shear Modulus
v_f = material['v_f']       #Fibre volume fraction
GIC = material['GIC']
GIIC = material['GIIC']
KP = material['KP']
SLP = material['SLP']
eta_y = material['eta_y']
eta_s = material['eta_s']
Em = material['Em']
Em_star = material['Em_star']
t_ply = material['t']
#
if OmniLPF == 'YES' :
    E2,G12,nu12 = get_DegradedProp(
        E2,nu12,G12,v_f,Em_star,Em,eta_y,eta_s)
    G23 = E2/(2.0+2.0*nu23)
    nu23 = nu23 * Em_star
    Omni_label = 'Omni LPF'
else:
    G23 = E2/(2.0+2.0*nu23)
    Omni_label = 'Omni FPF'

nu21 = nu12*E2/E1
nu13 = nu12
nu31 = nu21
nu32 = nu23
#
XT = material['XT']
XC = material['XC']
YT = material['YT']
YBT = material['YBT']
YC = material['YC']
YBC = material['YBC']
SL = material['SL']
ST = material['ST']
beta = material['beta'] 

DEG2RAD = math.pi / 180.0

# Failure parameters (from 3D invariant-based failure theory)
a1 = 1.0/(ST*ST)
a2 = 1.0/(SL*SL)
a32T = (1.0 - YT/(2.0*YBT)- \
        a1*(YT*YT)/4.0)/(YT*YT - \
        2.0 * YBT * YT)
a3T  = 1.0/(2.0*YBT) - 2.0*a32T*YBT
a32C = (1.0 - YC/(2.0*YBC)- \
        a1*(YC*YC)/4.0)/(YC*YC - \
        2.0 * YBC * YC)
a3C  = 1.0/(2.0*YBC) - 2.0*a32C*YBC

# PHIC (assuming phi0 = 0)
sqr = np.sqrt(a1 - 4.0*a2 + a2*a2*XC*XC + a3C*a3C + \
       2.0*a2*a3C*XC + 4.0*a32C)
PHIC = 0.5*np.arccos(((a1+ 4.0*a32C)*XC + 4.0*a3C \
                     + 4.0*sqr) /((a1 \
                     - 4.0*a2 + 4.0*a32C) * XC))
Kcompl = np.array([[(1/E1), -(nu12/E1), -(nu12/E1),         0,      0,          0],
                   [-(nu12/E1), (1/E2), -(nu23/E2),         0,      0,          0],
                   [-(nu12/E1), -(nu23/E2), (1/E2),         0,      0,          0],
                   [    0,          0,          0, (1/(1*G23)),     0,          0],
                   [    0,          0,          0,          0, (1/(1*G12)),     0],
                   [    0,          0,          0,          0,      0, (1/(1*G12))]])

Kstif = np.linalg.inv(Kcompl)
#Delta_k = (1-2*nu12*nu21-nu23*nu23-2*nu21*nu12*nu23)/(E1*E2*E2)
#Kstif = np.array([
#    [(1-nu23*nu23)/(E2*E2*Delta_k), (nu21+nu23*nu21)/(E2*E2*Delta_k), (nu21+nu23*nu21)/(E2*E2*Delta_k),  0,  0, 0],
#    [(nu21+nu23*nu21)/(E2*E2*Delta_k), (1-nu12*nu21)/(E1*E2*Delta_k), (nu23+nu12*nu21)/(E1*E2*Delta_k),  0,  0, 0],
#    [(nu21+nu23*nu21)/(E2*E2*Delta_k), (nu23+nu12*nu21)/(E1*E2*Delta_k), (1-nu12*nu21)/(E1*E2*Delta_k),  0,  0, 0],
#    [                       0,                          0,                          0,                  G23, 0, 0],
#    [                       0,                          0,                          0,                   0, G12,0],
#    [                       0,                          0,                          0,                   0, 0,G12]])
#
C11 = Kstif[0,0]
C12 = Kstif[0,1]
C13 = Kstif[0,2]
C22 = Kstif[1,1]
C23 = Kstif[1,2]
C33 = Kstif[2,2]
C44 = Kstif[3,3]
C55 = Kstif[4,4]
C66 = Kstif[5,5]
#
U1 = (3*C11+3*C22+2*C12+4*C66)/8
U2 = (C11 - C22)/2
U3 = (C11+C22-2*C12-4*C66)/8
U4 = (C11+C22+6*C12-4*C66)/8
U5 = (C11+C22-2*C12+4*C66)/8
U6 = C33
U7 = (C44+C55)/2
U8 = (C13+C23)/2

Q11x = []
Q12x = []
Q13x = []
Q22x = []
Q23x = []
Q33x = []
Q44x = []
Q55x = []
Q66x = []
Q45x = []
Q36x = []
Q16x = []
Q26x = []

for csi in Layup:
    csi_dg = csi*DEG2RAD
    Tr2 = np.array([
        [cos(csi_dg)* cos(csi_dg), sin(csi_dg)*sin(csi_dg),0.0, 0.0, 0.0, 0.5 * sin(2*csi_dg)],
        [sin(csi_dg)* sin(csi_dg), cos(csi_dg) *cos(csi_dg), 0.0, 0.0, 0.0, -0.5 * sin(2*csi_dg)],
        [0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, cos(csi_dg), -sin(csi_dg), 1.0],
        [0.0, 0.0, 0.0, sin(csi_dg), cos(csi_dg), 1.0],
        [-sin(2.0 * csi_dg), sin(2*csi_dg), 0.0, 0.0, 0.0, cos(2*csi_dg)]])
    Kstif3D = np.matmul(Kstif, Tr2)
    #
    m = cos(csi_dg)
    m_2 = cos(csi_dg)*cos(csi_dg)
    m_3 = m_2*cos(csi_dg)
    m_4 = m_2*m_2
    n = sin(csi_dg)
    n_2 = sin(csi_dg)*sin(csi_dg)
    n_3 = n_2*sin(csi_dg)
    n_4 = n_2*n_2
    Tr3 = np.array([
        [m_4,       n_4,        2*m_2*n_2,      m_2*n_2],
        [n_4,       m_4,        2*m_2*n_2,      m_2*n_2],
        [m_2*n_2,   m_2*n_2,    m_4+n_4,        -m_2*n_2],
        [4*m_2*n_2, 4*m_2*n_2,  -8*m_2*n_2,     (m_2-n_2)**2],
        [2*m_3*n,   -2*m*n_3,   2*m*n_3-2*m_3*n, m*n_3-m_3*n],
        [2*m*n_3,   -2*m_3*n,   2*m_3*n-2*m*n_3, m_3*n-m*n_3]])
    Tr3_ = np.array([C11, C22, C12, 4*C66])
    Qmatr = np.matmul(Tr3,Tr3_)
    
    Q11x.append(Qmatr[0])
    Q22x.append(Qmatr[1])
    Q12x.append(Qmatr[2])
    Q66x.append(Qmatr[3]/4)
    Q16x.append(Qmatr[4]/2)
    Q26x.append(Qmatr[5]/2)

    Tr4 = np.array([
        [m_2,     n_2,     0,       0,    0],
        [n_2,     m_2,     0,       0,    0],
        [m*n,     -m*n,    0,       0,    0],
        [0,       0,       m_2,     n_2,  0],
        [0,       0,       n_2,     m_2,  0],
        [0,       0,       -m*n,    n*m,  0],
        [0,       0,       0,       0,    1]])
    Tr4_ = np.array([C13, C23, C44, C55, C33])
    Qmatr2 = np.matmul(Tr4,Tr4_)
    
    Q13x.append(Qmatr2[0])
    Q23x.append(Qmatr2[1])
    Q36x.append(Qmatr2[2])
    Q44x.append(Qmatr2[3])
    Q55x.append(Qmatr2[4])
    Q45x.append(Qmatr2[5])
    Q33x.append(Qmatr2[6])
    '''
    Q11 = Kstif3D[0,0]
    Q11x.append(Q11)
    Q12 = Kstif3D[0,1]
    Q12x.append(Q12)
    Q13 = Kstif3D[0,2]
    Q13x.append(Q13)
    Q22 = Kstif3D[1,1]
    Q22x.append(Q22)
    Q33 = Kstif3D[2,2]
    Q33x.append(Q33)
    Q44 = Kstif3D[3,3]
    Q44x.append(Q44)
    Q55 = Kstif3D[4,4]
    Q55x.append(Q55)
    Q66 = Kstif3D[5,5]
    Q66x.append(Q66)
    Q23 = Kstif3D[1,2]
    Q23x.append(Q23)
    '''

if Layup == [+35.0, -35.0, -35.0, +35.0]:
    A11_star = (N35p/Ntot)*Q11x[0]+(N_35p/Ntot)*Q11x[1]
    A12_star = (N35p/Ntot)*Q12x[0]+(N_35p/Ntot)*Q12x[1]
    A13_star = (N35p/Ntot)*Q13x[0]+(N_35p/Ntot)*Q13x[1]
    A22_star = (N35p/Ntot)*Q22x[0]+(N_35p/Ntot)*Q22x[1]
    A23_star = (N35p/Ntot)*Q23x[0]+(N_35p/Ntot)*Q23x[1]
    A33_star = (N35p/Ntot)*Q33x[0]+(N_35p/Ntot)*Q33x[1]
    A44_star = (N35p/Ntot)*Q44x[0]+(N_35p/Ntot)*Q44x[1]
    A55_star = (N35p/Ntot)*Q55x[0]+(N_35p/Ntot)*Q55x[1]
    A66_star = (N35p/Ntot)*Q66x[0]+(N_35p/Ntot)*Q66x[1]
    A16_star = (N35p/Ntot)*Q16x[0]+(N_35p/Ntot)*Q16x[1]
    A26_star = (N35p/Ntot)*Q26x[0]+(N_35p/Ntot)*Q26x[1]
    A36_star = (N35p/Ntot)*Q36x[0]+(N_35p/Ntot)*Q36x[1]
    A45_star = (N35p/Ntot)*Q45x[0]+(N_35p/Ntot)*Q45x[1]
elif Layup == [+55.0, -55.0, -55.0, +55.0]:
    A11_star = (N55p/Ntot)*Q11x[0]+(N_55p/Ntot)*Q11x[1]
    A12_star = (N55p/Ntot)*Q12x[0]+(N_55p/Ntot)*Q12x[1]
    A13_star = (N55p/Ntot)*Q13x[0]+(N_55p/Ntot)*Q13x[1]
    A22_star = (N55p/Ntot)*Q22x[0]+(N_55p/Ntot)*Q22x[1]
    A23_star = (N55p/Ntot)*Q23x[0]+(N_55p/Ntot)*Q23x[1]
    A33_star = (N55p/Ntot)*Q33x[0]+(N_55p/Ntot)*Q33x[1]
    A44_star = (N55p/Ntot)*Q44x[0]+(N_55p/Ntot)*Q44x[1]
    A55_star = (N55p/Ntot)*Q55x[0]+(N_55p/Ntot)*Q55x[1]
    A66_star = (N55p/Ntot)*Q66x[0]+(N_55p/Ntot)*Q66x[1]
    A16_star = (N55p/Ntot)*Q16x[0]+(N_55p/Ntot)*Q16x[1]
    A26_star = (N55p/Ntot)*Q26x[0]+(N_55p/Ntot)*Q26x[1]
    A36_star = (N55p/Ntot)*Q36x[0]+(N_55p/Ntot)*Q36x[1]
    A45_star = (N55p/Ntot)*Q45x[0]+(N_55p/Ntot)*Q45x[1]
elif Layup == [+45.0, +0.0, -45.0, 90.0, 90.0, -45.0, +0.0,+45.0]:
    A11_star = (N0p/Ntot)*Q11x[0]+(N45p/Ntot)*Q11x[1]+(N_45p/Ntot)*Q11x[2]+  \
                (N90p/Ntot)*Q11x[3]
    A12_star = (N0p/Ntot)*Q12x[0]+(N45p/Ntot)*Q12x[1]+(N_45p/Ntot)*Q12x[2]+  \
                (N90p/Ntot)*Q12x[3]
    A13_star = (N0p/Ntot)*Q13x[0]+(N45p/Ntot)*Q13x[1]+(N_45p/Ntot)*Q13x[2]+  \
                (N90p/Ntot)*Q13x[3]
    A22_star = (N0p/Ntot)*Q22x[0]+(N45p/Ntot)*Q22x[1]+(N_45p/Ntot)*Q22x[2]+  \
                (N90p/Ntot)*Q22x[3]
    A23_star = (N0p/Ntot)*Q23x[0]+(N45p/Ntot)*Q23x[1]+(N_45p/Ntot)*Q23x[2]+  \
                (N90p/Ntot)*Q23x[3]
    A33_star = (N0p/Ntot)*Q33x[0]+(N45p/Ntot)*Q33x[1]+(N_45p/Ntot)*Q33x[2]+  \
                (N90p/Ntot)*Q33x[3]
    A44_star = (N0p/Ntot)*Q44x[0]+(N45p/Ntot)*Q44x[1]+(N_45p/Ntot)*Q44x[2]+  \
                (N90p/Ntot)*Q44x[3]
    A55_star = (N0p/Ntot)*Q55x[0]+(N45p/Ntot)*Q55x[1]+(N_45p/Ntot)*Q55x[2]+  \
                (N90p/Ntot)*Q55x[3]
    A66_star = (N0p/Ntot)*Q66x[0]+(N45p/Ntot)*Q66x[1]+(N_45p/Ntot)*Q66x[2]+  \
                (N90p/Ntot)*Q66x[3]
    A16_star = (N0p/Ntot)*Q16x[0]+(N45p/Ntot)*Q16x[1]+(N_45p/Ntot)*Q16x[2]+  \
                (N90p/Ntot)*Q16x[3]
    A26_star = (N0p/Ntot)*Q26x[0]+(N45p/Ntot)*Q26x[1]+(N_45p/Ntot)*Q26x[2]+  \
                (N90p/Ntot)*Q26x[3]
    A36_star = (N0p/Ntot)*Q36x[0]+(N45p/Ntot)*Q36x[1]+(N_45p/Ntot)*Q36x[2]+  \
                (N90p/Ntot)*Q36x[3]
    A45_star = (N0p/Ntot)*Q45x[0]+(N45p/Ntot)*Q45x[1]+(N_45p/Ntot)*Q45x[2]+  \
                (N90p/Ntot)*Q45x[3]
elif Layup == [90.0, +30.0, -30.0, -30.0, +30.0, 90.0]:
    A11_star = (0.172)*Q11x[0]+(0.414)*Q11x[1]+(0.414)*Q11x[2]
    A22_star = (0.172)*Q22x[0]+(0.414)*Q22x[1]+(0.414)*Q22x[2]
    A12_star = (0.172)*Q12x[0]+(0.414)*Q12x[1]+(0.414)*Q12x[2]
    A66_star = (0.172)*Q66x[0]+(0.414)*Q66x[1]+(0.414)*Q66x[2]
    A16_star = (0.172)*Q16x[0]+(0.414)*Q16x[1]+(0.414)*Q16x[2]
    A26_star = (0.172)*Q26x[0]+(0.414)*Q26x[1]+(0.414)*Q26x[2]
elif Layup == [90.0,90.0,+0.0,+0.0,90.0,90.0]:
    A11_star = (N90p/Ntot)*Q11x[0]+(N0p/Ntot)*Q11x[2]
    A22_star = (N90p/Ntot)*Q22x[0]+(N0p/Ntot)*Q22x[2]
    A12_star = (N90p/Ntot)*Q12x[0]+(N0p/Ntot)*Q12x[2]
    A66_star = (N90p/Ntot)*Q66x[0]+(N0p/Ntot)*Q66x[2]
    A16_star = (N90p/Ntot)*Q16x[0]+(N0p/Ntot)*Q16x[2]
    A26_star = (N90p/Ntot)*Q26x[0]+(N0p/Ntot)*Q26x[2]
elif Layup == [+0.0,+0.0,90.0,90.0,+0.0,+0.0]:
    A11_star = (N0p/Ntot)*Q11x[0]+(N90p/Ntot)*Q11x[2]
    A22_star = (N0p/Ntot)*Q22x[0]+(N90p/Ntot)*Q22x[2]
    A12_star = (N0p/Ntot)*Q12x[0]+(N90p/Ntot)*Q12x[2]
    A66_star = (N0p/Ntot)*Q66x[0]+(N90p/Ntot)*Q66x[2]
    A16_star = (N0p/Ntot)*Q16x[0]+(N90p/Ntot)*Q16x[2]
    A26_star = (N0p/Ntot)*Q26x[0]+(N90p/Ntot)*Q26x[2]
#

TOLER = 1.0e-06
tol4 = 1.0e-05

#s11=s22=s33=s12=s13=s23=0.0

# PUT STRESS ON X AXIS SX and set the label SXlab
SX = 'e11'
SXlab = r'$\epsilon_{11}$'
# PUT STRESS ON Y AXIS SY and set the label SYlab
SY = 'e22'
SYlab = r'$\epsilon_{22}$'
# PUT STRESS ON Z AXIS SZ. IF A 2D ENVELOPE IS REQUIRED PUT None
SZ = None

# INPUT FOR THE OMNI STRAIN ENVELOPE

deltaTheta = 15 # choose a value from 1 to 90[deg] (default = 15) to set at which 
                # ply orientation you want a failure envelope in strain space. It highly
                # affects the accuracy of the resulting omni strain FPF envelope 
                # (e.g. with 15[deg] you get 7 failure envelopes, for the following 
                # ply orientations: 0-15-30-45-60-75-90[deg]).
                # If you just want a single failure envelope, set deltaTheta = 90
                # and comment: "# ThetaS.append(thetaRR)"

# --------------------------------------------------------------------------- #

# ---------------------------- local functions ------------------------------ # 

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


def get_SXX_SYY(R, GAMM, SX, SY):
    e11 = e22 = e33 = e12 = e23 = e13 = 0
    SXX = R * cos(GAMM)
    SYY = R * sin(GAMM)

    e11, e22, e33, e12, e23, e13 = get_sij(SX, SXX, e11, e22, e33, e12, e23, e13)
    e11, e22, e33, e12, e23, e13 = get_sij(SY, SYY, e11, e22, e33, e12, e23, e13)

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

def get_delta_F(R, GAMM, SX, SY,
                TOLER, DEG2RAD, a1, a2, a32T, a3T, a32C, a3C, PHIC, XT, XC,
                YC, SL, ST, beta, G12):
    SXX, SYY, s11, s22, s33, s12, s23, s13 = get_SXX_SYY(R, GAMM, SX, SY)

    F_MAX, FLAG, THETADEG, PHI_D, PSI_D = get_F_max(s11, s22, s33, s12, s13, s23, nu12,
                                TOLER, DEG2RAD, a1, a2, a32T, a3T, a32C, a3C, 
                                PHIC, XT, XC, YC, SL, ST, beta, G12)
    DELTA_F = 1.0 - F_MAX

    return DELTA_F, FLAG, THETADEG, PHI_D, PSI_D, SXX, SYY

# --------------------------------------------------------------------------- #

# ---------------------- Failure loci computation --------------------------- #
start = time.time()
thetaR = 0.0

if OmniLPF == 'YES':
    if INSITU == 'YES':
        RMIN = 0
        RMAX = 0.45
    else:
        RMIN = 0
        RMAX = 0.25
else:
    if INSITU == 'YES':
        RMIN = 0
        RMAX = 0.08
    else:
        RMIN = 0
        RMAX = 0.08

IMAX = 5000

# OUTFILE= open("outfile.txt","w+") #create, write, open and read a file
# OUTFILE.write("Envelope outfile")
# OUTFILE.close()
S_out = []
SXSY = []
SX_out = []
SY_out = []
S_out12 = []
IMAXfail = []

thetaRR = 90
ThetaS = list(range(0, thetaRR, deltaTheta))
ThetaS.append(thetaRR)
for thetaRdeg in ThetaS:
    S_out.append([])
    SXSY.append([])
    thetaR = thetaRdeg * DEG2RAD

    T1inv = np.array([[cos(thetaR) * cos(thetaR), sin(thetaR) * sin(thetaR),
                      0.0, 0.0, 0.0, -sin(2.0 * thetaR)],
                      [sin(thetaR) * sin(thetaR), cos(thetaR) * cos(thetaR),
                      0.0, 0.0, 0.0, sin(2.0 * thetaR)],
                      [0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
                      [0.0, 0.0, 0.0, cos(thetaR), -sin(thetaR), 1.0],
                      [0.0, 0.0, 0.0, sin(thetaR), cos(thetaR), 1.0],
                      [0.5 * sin(2.0 * thetaR), -0.5 * sin(2.0 * thetaR), 
                      0.0, 0.0, 0.0, cos(2.0 * thetaR)]])

    T2 = np.array([[cos(thetaR) * cos(thetaR), sin(thetaR) * sin(thetaR), 
                   0.0, 0.0, 0.0, 0.5 * sin(2.0 * thetaR)],
                   [sin(thetaR) * sin(thetaR), cos(thetaR) * cos(thetaR), 
                   0.0, 0.0, 0.0, -0.5 * sin(2.0 * thetaR)],
                   [0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
                   [0.0, 0.0, 0.0, cos(thetaR), -sin(thetaR), 1.0],
                   [0.0, 0.0, 0.0, sin(thetaR), cos(thetaR), 1.0],
                   [-sin(2.0 * thetaR), sin(2.0 * thetaR), 0.0, 0.0, 0.0, cos(2.0 * thetaR)]])

    ENNE = 3
    ENNE_deg = ENNE * 360
    GAMMAS = list(range(0, ENNE_deg, 1))
    GAMMAS.append(ENNE_deg)
    for gamma_i in GAMMAS:

        print('gamma_i=', gamma_i)
        GAMMA_DEG = gamma_i / ENNE
        GAMM = GAMMA_DEG * DEG2RAD

        DELTA_F = tol4 + 1.0
        iterI = 0

        a = RMIN
        DELTA_F_a, _, _, _, _, _, _ = get_delta_F(a, GAMM, SX, SY,
                                            TOLER, DEG2RAD, a1, a2, a32T, a3T,
                                            a32C, a3C, PHIC, XT, XC,
                                            YC, SL, ST, beta, G12)
        b = RMAX
        DELTA_F_b, _, _, _, _, _, _ = get_delta_F(b, GAMM, SX, SY,
                                            TOLER, DEG2RAD, a1, a2, a32T, a3T, 
                                            a32C, a3C, PHIC, XT, XC,
                                            YC, SL, ST, beta, G12)

        while np.abs(DELTA_F) > tol4 and iterI < IMAX:

            # compute the value for the mean interval value
            c = (a + b) / 2.0
            DELTA_F, FLAG, THETADEG, PHI_D, PSI_D, SXX, SYY = get_delta_F(c, GAMM, SX, SY,
                                                            TOLER, DEG2RAD, a1, a2, a32T, 
                                                            a3T, a32C, a3C, PHIC, XT, XC,
                                                            YC, SL, ST, beta, G12)

            # get new bounds
            if DELTA_F >= 0:
                a = c
            else:
                b = c

            # print('F_MAX=',F_MAX)
            #print('FLAG=', FLAG)
            iterI += 1

        if iterI < IMAX:

            if gamma_i>2:
                diff1 =abs(SYY)-np.abs(SY_out[-1])
                diff2 =abs(SXX)-np.abs(SX_out[-1])
                if diff1> 0.5 or diff2> 0.5:
                    S_out[-1].append([(SX_out[-1][0]),(SY_out[-1][0]),FLAG, PHI_D, PSI_D])
                    s11 = (SX_out[-1][0])*A11_star+(SY_out[-1][0])*A12_star
                    s22 = (SX_out[-1][0])*A12_star+(SY_out[-1][0])*A22_star
                    s12 = (SX_out[-1][0])*A16_star+(SY_out[-1][0])*A66_star
                    #
                    '''
                    Eps = np.array([(SX_out[-1][0]), (SY_out[-1][0]), 0.0, 0.0, 0.0, 0.0])
                    Sk = np.matmul(Kstif, Eps)
                    s11 = Sk[0]
                    s22 = Sk[1]
                    s33 = Sk[2]
                    s23 = Sk[3]
                    s13 = Sk[4]
                    s12 = Sk[5]
                    '''
                    SXSY[-1].append([s11, s22])
                else:
                    S_out[-1].append([SXX, SYY, FLAG, PHI_D, PSI_D])
                    SX_out.append([SXX])
                    SY_out.append([SYY])
                    s11 = SXX*A11_star+SYY*A12_star
                    s22 = SXX*A12_star+SYY*A22_star
                    s12 = SXX*A16_star+SYY*A66_star
                    '''
                    Eps = np.array([SXX, SYY, 0.0, 0.0, 0.0, 0.0])
                    Sk = np.matmul(Kstif, Eps)
                    s11 = Sk[0]
                    s22 = Sk[1]
                    s33 = Sk[2]
                    s23 = Sk[3]
                    s13 = Sk[4]
                    s12 = Sk[5]
                    '''
                    SXSY[-1].append([s11, s22, s12])
            else:
                S_out[-1].append([SXX, SYY, FLAG, PHI_D, PSI_D])
                SY_out.append([SYY])
                SX_out.append([SXX])
                s11 = SXX*A11_star+SYY*A12_star
                s22 = SXX*A12_star+SYY*A22_star
                s12 = SXX*A16_star+SYY*A66_star
                '''
                Eps = np.array([SXX, SYY, 0.0, 0.0, 0.0, 0.0])
                Sk = np.matmul(Kstif, Eps)
                s11 = Sk[0]
                s22 = Sk[1]
                s33 = Sk[2]
                s23 = Sk[3]
                s13 = Sk[4]
                s12 = Sk[5]
                '''
                SXSY[-1].append([s11, s22, s12])
        else:
            print('IMAX=', iterI)
            IMAXfail.append([SXX, SYY, FLAG, PHI_D, PSI_D])
            SY_out.append([SYY])
            SX_out.append([SXX])
            S_out[-1].append([(SX_out[-1][0]),(SY_out[-1][0]),FLAG, PHI_D, PSI_D])
            Eps = np.array([(SX_out[-1][0]), (SY_out[-1][0]), 0.0, 0.0, 0.0, 0.0])
            Sk =  1/(8*t_ply) * np.matmul(Kstif, Eps)
            s11 = Sk[0]
            s22 = Sk[1]
            s33 = Sk[2]
            s23 = Sk[3]
            s13 = Sk[4]
            s12 = Sk[5]
            SXSY[-1].append([s11, s22])
            #
            #Eps = np.array([SXX, SYY, 0.0, 0.0, 0.0, 0.0])
            #Sk = np.matmul(Kstif, Eps)
            #s11 = Sk[0]
            #s22 = Sk[1]
            #s33 = Sk[2]
            #s23 = Sk[3]
            #s13 = Sk[4]
            #s12 = Sk[5]
            #S_out12.append([s11, s22, s33, s23, s13, s12, SXX, SYY])

#----------------------------- PLOT -----------------------------------#
#print('S_out12=',S_out12)
#print('IMAXfail=', IMAXfail)

#S_out = [x for x in S_out if x != []]
#S_out = list(filter(None,S_out))

print(len(S_out),np.shape(S_out))
#print(len(SY_out),np.shape(SY_out))
SXX = []
SYY = []
FLAGS = []
PHIS = []
PSIS = []
#for SS_out in S_out:
#    if SS_out != [[None,None]]:
#        SXX.append([value[0] for value in SS_out if value != None])
#        SYY.append([value[1] for value in SS_out if value != None])
for SS_out in S_out:
    SXX.append([value[0] for value in SS_out])
    SYY.append([value[1] for value in SS_out])
    FLAGS.append([value[2] for value in SS_out])
    PHIS.append([round(value[3],1) for value in SS_out])
    PSIS.append([round(value[4],1) for value in SS_out])

#SXX = [x for x in SXX if (math.isnan(x[0]) == True)]
#SYY = [x for x in SYY if (math.isnan(x[1]) == True)]
SSXX = []
SSYY = []
SSXXmin = []
SSYYmin = []
for SXSY_out in SXSY:
    SSXX.append([value[0] for value in SXSY_out])
    SSYY.append([value[1] for value in SXSY_out])


SXXmin = []
SYYmin = []
FLAGcrit = []
PHIcrit = []
PSIcrit = []

for i in range(len(SXX[0])): #len(SXX[0])= failure points for each envelope
    vecXX = np.zeros(len(S_out)) #len(S_out) = number of envelopes
    vecYY = np.zeros(len(S_out))
    vecFLG = np.zeros(len(S_out))
    vecPHI = np.zeros(len(S_out))
    vecPSI = np.zeros(len(S_out))
    vecSXX = np.zeros(len(S_out))
    vecSYY = np.zeros(len(S_out))
    for j in range(len(S_out)):
        vecXX[j] = SXX[j][i]
        vecYY[j] = SYY[j][i]
        vecFLG[j] = FLAGS[j][i]
        vecPHI[j] = PHIS[j][i]
        vecPSI[j] = PSIS[j][i]
        vecSXX[j] = SSXX[j][i]
        vecSYY[j] = SSYY[j][i]
    iminXX = np.argmin(abs(vecXX))# argmin returns the indices of the min value along an axis
    iminYY = np.argmin(abs(vecYY))
    SXXmin.append(vecXX[iminXX])
    SYYmin.append(vecYY[iminYY])
    SSXXmin.append(vecSXX[iminXX])
    SSYYmin.append(vecSYY[iminYY])
    FLAGcrit.append(vecFLG[np.minimum(iminXX,iminYY)])
    PHIcrit.append(vecPHI[np.minimum(iminXX,iminYY)])
    PSIcrit.append(vecPSI[np.minimum(iminXX,iminYY)])

colours = ['r', 'g', 'b', 'yellowgreen', 'm', 'c', 'navy', 'peru',
           'darkslateblue', 'goldenrod']

FLAGmarker = []
FLAGcolour = []
'''
for x in FLAGcrit:
    if x == -1 or x ==  1:
        # MATRIX-DOMINATED FAILURE
        FLAGmarker.append('circle')# Type of markers:automatic, none, square, diamond
                                   # triangle, x, star, short_dash, long_dash, circle, plus
        FLAGcolour.append('green')# List of colours: black, blue, brown, cyan, gray, green, lime,
                                   # magenta, navy, orange, pink, purple, red, silver, white, yellow
    else:
        # FIBER-DOMINATED FAILURE
        FLAGmarker.append('diamond')
        FLAGcolour.append('red')
'''
for x in FLAGcrit:
    if x ==  1:
        # MATRIX-DOMINATED FAILURE (Transverse cracking) when I3>0
        FLAGmarker.append('circle')# Type of markers:automatic, none, square, diamond
                                   # triangle, x, star, short_dash, long_dash, circle, plus
        FLAGcolour.append('green')# List of colours: black, blue, brown, cyan, gray, green, lime,
                                   # magenta, navy, orange, pink, purple, red, silver, white, yellow
    elif x ==  -1:
        # MATRIX-DOMINATED FAILURE (Transverse cracking) when I3<0
        FLAGmarker.append('circle')
        FLAGcolour.append('red')
    elif x ==  3:
        # FIBRE-DOMINATED FAILURE when I3>0
        FLAGmarker.append('x')
        FLAGcolour.append('green')
    elif x ==  -3:
        # FIBRE-DOMINATED FAILURE when I3<0
        FLAGmarker.append('x')
        FLAGcolour.append('red')
    elif x ==  2:
        # FIBRE-TENSILE FAILURE
        FLAGmarker.append('triangle')
        FLAGcolour.append('cyan')

if PlotFigure == 'YES':

    fig = plt.figure(figsize=(5.5,4.0))
    ax = fig.add_subplot(1, 1, 1)
    ax.grid(True, which='both')
    #ax.set_aspect('equal')
    plt.xlabel(SXlab, fontsize=13)
    plt.ylabel(SYlab, fontsize=13)
    ax.axhline(y=0, color='k', linewidth=1.0)
    ax.axvline(x=0, color='k', linewidth=1.0)
    for i in range(len(SXX)):
        plt.plot(SXX[i], SYY[i], colours[i], linewidth=2.0,
                label=' ' + str(i * deltaTheta) + 'º ply')
    if INSITU == 'YES':
        plt.plot(SXXmin, SYYmin, color='k', linewidth=3.0, label=' '+Omni_label + '_insitu_t0'+str(t_ply)[2:])
    else:
        plt.plot(SXXmin, SYYmin, color='k', linewidth=3.0, label=' '+Omni_label)

    ax.legend(loc='upper center', bbox_to_anchor=(1.25, 0.75))

    if SAVEFIG == 'YES':
        fig.savefig('Figure_'+Omni_label+'_Inv3DFC_'+SX+SY+'_'+mat_label+'.png', dpi=300, bbox_inches="tight")
    else:
        pass
        
    plt.show()

#------------------------- SAVE DATA -------------------------#
#out2 = open("outfile2.txt", "w+")
#out2.write(str(S_out12))
#out2.close()

#for i in range(len(SXX)):
#    saveSXY = np.savetxt('Env'+str(i)+'Data.txt', (SXX[i],SYY[i]), delimiter =', ')
#saveOmni = np.savetxt('EnvOmniData.txt', (SXXmin,SYYmin), delimiter =', ')

#------------------- EXCEL Postprocessing --------------------#

if ExcelOutput == 'YES':
    workbook   = xlsxwriter.Workbook(Omni_label+'_Inv3DFC_'+mat_label+'_[90,0,90].xlsx')
    worksheet1 = workbook.add_worksheet()
    worksheet2 = workbook.add_worksheet('OmniEnv')
    chart =  workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
    if INSITU == 'YES':
        plylabel = 'º ply_insitu_t0'+str(t_ply)[2:]
    else:
        plylabel = 'º ply'
    for i in range(len(SXX)):
        worksheet2.write_string(2, int(4*i), SX+'_Env'+str(i))
        worksheet2.write_column(3, int(4*i), SXX[i])
        worksheet2.write_string(2, int(1+4*i), SY+'_Env'+str(i))
        worksheet2.write_column(3, int(1+4*i), SYY[i])
        worksheet2.write_string(2, int(2+4*i), 'PHI-kink.angle')
        worksheet2.write_column(3, int(2+4*i), PHIS[i])
        worksheet2.write_string(2, int(3+4*i), 'PSI-angle of kink.plane')
        worksheet2.write_column(3, int(3+4*i), PSIS[i])
        worksheet2.write_string(1, int(4*i), str(i*deltaTheta) + plylabel)
        chart.add_series({
            'name':       ['OmniEnv', 1, int(4*i)]  ,
            'categories': ['OmniEnv', 3, int(4*i), int(ENNE_deg+3), int(4*i)],
            'values':     ['OmniEnv', 3, int(1+4*i), int(ENNE_deg+3), int(1+4*i)],
            })

    worksheet2.write_string(2, int(4*len(SXX)+1), SX+'_OmniEnv')
    worksheet2.write_column(3, int(4*len(SXX)+1), SXXmin)
    worksheet2.write_string(2, int(4*len(SXX)+2), SY+'_OmniEnv')
    worksheet2.write_column(3, int(4*len(SXX)+2), SYYmin)
    worksheet2.write_string(1, int(4*len(SXX)+1), Omni_label+'_Inv3DFC')
    worksheet2.write_string(2, int(4*len(SXX)+3), 'Flags_OmniEnv')
    worksheet2.write_column(3, int(4*len(SXX)+3), FLAGcrit)
    worksheet2.write_column(3, int(4*len(SXX)+4), FLAGmarker)
    worksheet2.write_column(3, int(4*len(SXX)+5), FLAGcolour)
    worksheet2.write_string(2, int(4*len(SXX)+7), 'S'+SX[1:]+'_OmniEnv')
    worksheet2.write_column(3, int(4*len(SXX)+7), SSXXmin)
    worksheet2.write_string(2, int(4*len(SXX)+8), 'S'+SY[1:]+'_OmniEnv')
    worksheet2.write_column(3, int(4*len(SXX)+8), SSYYmin)
    worksheet2.write_string(1, int(4*len(SXX)+7), Omni_label+'_Inv3DFC')
    worksheet2.write_column(3, int(4*len(SXX)+9), PHIcrit)
    worksheet2.write_string(2, int(4*len(SXX)+9), 'PHI-kink.angle')
    worksheet2.write_column(3, int(4*len(SXX)+10), PSIcrit)
    worksheet2.write_string(2, int(4*len(SXX)+10), 'PSI-angle of kink.plane')

    chart.add_series({
        'name':       ['OmniEnv',1, int(4*len(SXX)+1)] ,  
        'categories': ['OmniEnv', 3, int(4*len(SXX)+1), int(ENNE_deg+3), int(4*len(SXX)+1)],
        'values':     ['OmniEnv', 3, int(4*len(SXX)+2), int(ENNE_deg+3), int(4*len(SXX)+2)],
        'line':       {'color': 'black'},
        })
    
    chart.set_x_axis({'name': 'ε₁₁',
                      'name_font': {
                            'name': 'Calibri'},
                      'label_position': 'low',
                      'line': {'color': 'black','width': 1.00},
                      'major_gridlines': {
                            'visible': True,
                            'line': {'width': 0.75, 'dash_type': 'solid'}},
                      'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                      'minor_tick_mark': 'none'
                        })
    chart.set_y_axis({'name': 'ε₂₂',
                      'label_position': 'low',
                      'line': {'color': 'black','width': 1.00},
                      'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                      'minor_tick_mark': 'none'})
    worksheet2.insert_chart(3, int(4*len(SXX)+12), chart,{'x_scale': 1.5, 'y_scale': 1.5} )
    ######### CHART 1 ##########
    chart1 =  workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
    chart1.add_series({
        'name':       ['OmniEnv',1, int(4*len(SXX)+7)] ,  
        'categories': ['OmniEnv', 4, int(4*len(SXX)+7), int(ENNE_deg+3), int(4*len(SXX)+7)],
        'values':     ['OmniEnv', 4, int(4*len(SXX)+8), int(ENNE_deg+3), int(4*len(SXX)+8)],
        'line':       {'color': 'black'},
        })
    
    chart1.set_x_axis({'name': 'σ₁₁',
                      'name_font': {
                            'name': 'Calibri'},
                      'label_position': 'low',
                      'line': {'color': 'black','width': 1.00},
                      'major_gridlines': {
                            'visible': True,
                            'line': {'width': 0.75, 'dash_type': 'solid'}},
                      'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                      'minor_tick_mark': 'none'
                        })
    chart1.set_y_axis({'name': 'σ₂₂',
                      'label_position': 'low',
                      'line': {'color': 'black','width': 1.00},
                      'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                      'minor_tick_mark': 'none'})
    worksheet2.insert_chart(25, int(4*len(SXX)+12), chart1,{'x_scale': 1.5, 'y_scale': 1.5} )

    ######### CHART 2 ##########
    worksheet3 = workbook.add_worksheet('OmniEnv_Chart')
    worksheet3.write_string(0, 15, 'LEGEND:')
    cell_format = workbook.add_format({'bold': True, 'font_color': 'cyan'})
    worksheet3.write_string(1, 15, 'triangle(cyan): Fibre tensile failure',cell_format)
    cell_format = workbook.add_format({'bold': True, 'font_color': 'green'})
    worksheet3.write_string(2, 15, 'x(green): Fibre kinking when I3>0',cell_format)
    worksheet3.write_string(4, 15, 'o(green): Matrix cracking when I3>0',cell_format)
    cell_format = workbook.add_format({'bold': True, 'font_color': 'red'})
    worksheet3.write_string(3, 15, 'x(red): Fibre kinking when I3<0',cell_format)
    worksheet3.write_string(5, 15, 'o(red): Matrix cracking when I3<0',cell_format)
    chart2 =  workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
    chart2.add_series({
        'name':       ['OmniEnv',1, int(4*len(SXX)+1)] ,  
        'categories': ['OmniEnv', 3, int(4*len(SXX)+1), int(ENNE_deg+3), int(4*len(SXX)+1)],
        'values':     ['OmniEnv', 3, int(4*len(SXX)+2), int(ENNE_deg+3), int(4*len(SXX)+2)],
        'line':       {'color': 'black','width': 1.00},
        })
    delta_value = 4*ENNE
    for i in range(0,len(FLAGmarker),delta_value):
        if FLAGmarker[i] == 'x':
            chart2.add_series({
                'categories': ['OmniEnv', 3+i, int(4*len(SXX)+1), int(i+3), int(4*len(SXX)+1)],
                'values':     ['OmniEnv', 3+i, int(4*len(SXX)+2), int(i+3), int(4*len(SXX)+2)],
                'marker':     {
                     'type': FLAGmarker[i],
                     'size': 5,
                     'border': {'color': FLAGcolour[i]},
                     'fill':   {'none': True},},})
        elif FLAGmarker[i] == 'triangle':
            chart2.add_series({
                'categories': ['OmniEnv', 3+i, int(4*len(SXX)+1), int(i+3), int(4*len(SXX)+1)],
                'values':     ['OmniEnv', 3+i, int(4*len(SXX)+2), int(i+3), int(4*len(SXX)+2)],
                'marker':     {
                     'type': FLAGmarker[i],
                     'size': 5,
                     'border': {'none': True},
                     'fill':   {'color': FLAGcolour[i]},},})
        else:
            chart2.add_series({
                'categories': ['OmniEnv', 3+i, int(4*len(SXX)+1), int(i+3), int(4*len(SXX)+1)],
                'values':     ['OmniEnv', 3+i, int(4*len(SXX)+2), int(i+3), int(4*len(SXX)+2)],
                'marker':     {
                     'type': FLAGmarker[i],
                     'size': 5,
                     'border': {'color': FLAGcolour[i]},
                     'fill':   {'none': True},},})
                     #'border': {'none': True},
                     #'fill':   {'color': FLAGcolour[i]},},})
    chart2.set_x_axis({'name': 'ε₁₁',
                      'name_font': {
                            'name': 'Calibri'},
                      'label_position': 'low',
                      'line': {'color': 'black','width': 1.00},
                      'major_gridlines': {
                            'visible': True,
                            'line': {'width': 0.75, 'dash_type': 'solid'}},
                      'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                      'minor_tick_mark': 'none',
                        })
    chart2.set_y_axis({'name': 'ε₂₂',
                      'label_position': 'low',
                      'line': {'color': 'black','width': 1.00},
                      'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                      'minor_tick_mark': 'none'})
    chart2.set_legend({'position': 'none'})
    worksheet3.insert_chart(1, 3, chart2,{'x_scale': 1.5, 'y_scale': 1.5} )

    ######### CHART 3 ##########
    chart3 =  workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
    delta_value=3*ENNE
    for i in range(0,len(FLAGmarker),delta_value):
        chart3.add_series({
            'categories': ['OmniEnv', 3+i, int(4*len(SXX)+1), int(i+3+delta_value), int(4*len(SXX)+1)],
            'values':     ['OmniEnv', 3+i, int(4*len(SXX)+2), int(i+3+delta_value), int(4*len(SXX)+2)],
            'line':       {'color': FLAGcolour[i],'width': 1.00},
            })
    chart3.set_x_axis({'name': 'ε₁₁',
                      'name_font': {
                            'name': 'Calibri'},
                      'label_position': 'low',
                      'line': {'color': 'black','width': 1.00},
                      'major_gridlines': {
                            'visible': True,
                            'line': {'width': 0.75, 'dash_type': 'solid'}},
                      'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                      'minor_tick_mark': 'none',
                        })
    chart3.set_y_axis({'name': 'ε₂₂',
                      'label_position': 'low',
                      'line': {'color': 'black','width': 1.00},
                      'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                      'minor_tick_mark': 'none'})
    chart3.set_legend({'position': 'none'})
    worksheet3.insert_chart(30, 3, chart3,{'x_scale': 1.5, 'y_scale': 1.5} )
    workbook.close()
else:
    pass

print('IMAXfail=', IMAXfail)
elapsed = (time.time() - start)
print('elapsed time=',round(elapsed,2),'sec')
print('END')
