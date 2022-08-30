# --------------------------- Imports -------------------------------- #
from numpy import *
import numpy as np
import math
import time
from matplotlib import pyplot as plt
import xlsxwriter

print('loading your inputs...')
from your_inputs import *
#
from _includes.Inv3DFCFunctions import FailureMatrix
from _includes.Inv3DFCFunctions import FailureFiber
from _includes.DegradedPropEvaluation import get_DegradedProp
from _includes.omni_env_script import *
#
# -------------------------------------------------------------------- #
#
print('computing the omni-strain failure envelopes...')
#
if Symmetric == 'YES':
    N_sym = len(Layup)
    for xnum in range(N_sym):
        xnum = xnum + 1
        Layup.append(Layup[N_sym-xnum])
else:
    pass
Ntot = len(Layup)
#
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
#
if PlaneStress_2D == 'YES':

    # 2D Stiffness matrix
    Qxx = E1/(1-(nu12*nu21))
    Qyy = E2/(1-(nu12*nu21))
    Qxy = nu12*Qyy
    Qss = G12

    Kstif2D = np.array([
        [   Qxx,    Qxy,       0],
        [   Qxy,    Qyy,       0],
        [     0,      0,     Qss]])
    Kstif2D_array = np.array([Qxx, Qyy, Qxy, Qss])
    Kstif = np.array([
        [Qxx,   Qxy,   0,   0,   0,   0],
        [Qxy,   Qyy,   0,   0,   0,   0],
        [  0,     0,   0,   0,   0,   0],
        [  0,     0,   0,   0,   0,   0],
        [  0,     0,   0,   0,   0,   0],
        [  0,     0,   0,   0,   0, Qss]])
    #print('K_2D = ',Kstif2D)
    U1 = (3*Qxx+3*Qyy+2*Qxy+4*Qss)/8
    U2 = (Qxx - Qyy)/2
    U3 = (Qxx+Qyy-2*Qxy-4*Qss)/8
    U4 = (Qxx+Qyy+6*Qxy-4*Qss)/8
    U5 = (Qxx+Qyy-2*Qxy+4*Qss)/8
else:
    # 3D Stiffness matrix for transversely isotropic materials
    Delta_k = (1-2*nu12*nu21-nu23*nu23-2*nu21*nu12*nu23)/(E1*E2*E2)
    C11 = (1-nu23*nu23)/(E2*E2*Delta_k)
    C12 = (nu21+nu23*nu21)/(E2*E2*Delta_k)
    C13 = C12
    C22 = (1-nu12*nu21)/(E1*E2*Delta_k)
    C33 = C22
    C23 = (nu23+nu12*nu21)/(E1*E2*Delta_k)
    C44 = G23
    C55 = G12
    C66 = G12
    Kstif = np.array([
        [C11, C12, C13,   0,   0,   0],
        [C12, C22, C23,   0,   0,   0],
        [C13, C23, C33,   0,   0,   0],
        [  0,   0,   0, C44,   0,   0],
        [  0,   0,   0,   0, C55,   0],
        [  0,   0,   0,   0,   0, C66]])
    #
    U1 = (3*C11+3*C22+2*C12+4*C66)/8
    U2 = (C11 - C22)/2
    U3 = (C11+C22-2*C12-4*C66)/8
    U4 = (C11+C22+6*C12-4*C66)/8
    U5 = (C11+C22-2*C12+4*C66)/8
    U6 = C33
    U7 = (C44+C55)/2
    U8 = (C13+C23)/2
    #
    Tr_C = C11+C22+C33 + 2*(C44+C55+C66)

Q11x,Q12x,Q13x,Q22x,Q23x,Q33x,Q44x,Q55x,Q66x,Q45x,Q36x,Q16x,Q26x = \
     ([] for i in range(13))

Orientations = []
Layup_new = []
Z_delta = []

for itheta in Layup:
    if (-90 < itheta <= 90):
        pass
    elif itheta == -90:
        itheta = 90
    elif 90 < itheta <= 180:
        itheta = itheta - 90
    elif -180 <= itheta < -90:
        itheta = itheta + 90
    else:
        print('Error in Layup definition! Please make sure all angles',
        ' are within a range of -90 < theta < 180')
    if itheta not in Orientations:
        Orientations.append(itheta)
    else:
        pass
    Layup_new.append(itheta)

for itheta in Orientations:
    Ncount = Layup_new.count(itheta)
    Z_delta.append(Ncount/Ntot)

A11,A12,A13,A22,A23,A33,A44,A55,A66,A16,A26,A36,A45 = (0 for i in 
                                                            range(13))
x_ori = -1

for csi in Orientations:
    x_ori = x_ori + 1
    csi_dg = csi*DEG2RAD
    #
    m = cos(csi_dg)
    m_2 = cos(csi_dg)*cos(csi_dg)
    m_3 = m_2*cos(csi_dg)
    m_4 = m_2*m_2
    n = sin(csi_dg)
    n_2 = sin(csi_dg)*sin(csi_dg)
    n_3 = n_2*sin(csi_dg)
    n_4 = n_2*n_2
    #
    if  PlaneStress_2D == 'YES':

        Tr2 = np.array([
            [m_4,       n_4,        2*m_2*n_2,        4*m_2*n_2],
            [n_4,       m_4,        2*m_2*n_2,        4*m_2*n_2],
            [m_2*n_2,   m_2*n_2,     m_4+n_4,        -4*m_2*n_2],
            [m_2*n_2,   m_2*n_2,  -8*m_2*n_2,       (m_2-n_2)**2],
            [m_3*n,     -m*n_3,   2*m*n_3-2*m_3*n, 2*(m*n_3-m_3*n)],
            [m*n_3,     -m_3*n,   2*m_3*n-2*m*n_3, 2*(m_3*n-m*n_3)]])
        Qmatr2D = np.matmul(Tr2, Kstif2D_array)
        Q11x.append(Qmatr2D[0])
        Q22x.append(Qmatr2D[1])
        Q12x.append(Qmatr2D[2])
        Q66x.append(Qmatr2D[3])
        Q16x.append(Qmatr2D[4])
        Q26x.append(Qmatr2D[5])
        Q13x.append(0)
        Q23x.append(0)
        Q36x.append(0)
        Q44x.append(0)
        Q55x.append(0)
        Q45x.append(0)
        Q33x.append(0)
    else:
        Tr3 = np.array([
            [m_4,       n_4,        2*m_2*n_2,      m_2*n_2],
            [n_4,       m_4,        2*m_2*n_2,      m_2*n_2],
            [m_2*n_2,   m_2*n_2,    m_4+n_4,        -m_2*n_2],
            [4*m_2*n_2, 4*m_2*n_2,  -8*m_2*n_2,     (m_2-n_2)**2],
            [2*m_3*n,   -2*m*n_3,   2*m*n_3-2*m_3*n, m*n_3-m_3*n],
            [2*m*n_3,   -2*m_3*n,   2*m_3*n-2*m*n_3, m_3*n-m*n_3]])
        Tr3_ = np.array([C11, C22, C12, 4*C66])
        Qmatr = np.matmul(Tr3,Tr3_)
        #
        Q11x.append(Qmatr[0])
        Q22x.append(Qmatr[1])
        Q12x.append(Qmatr[2])
        Q66x.append(Qmatr[3]/4)
        Q16x.append(Qmatr[4]/2)
        Q26x.append(Qmatr[5]/2)
        #
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
        #
        Q13x.append(Qmatr2[0])
        Q23x.append(Qmatr2[1])
        Q36x.append(Qmatr2[2])
        Q44x.append(Qmatr2[3])
        Q55x.append(Qmatr2[4])
        Q45x.append(Qmatr2[5])
        Q33x.append(Qmatr2[6])
    #
    A11 = A11 + Z_delta[x_ori]*Q11x[-1]
    A12 = A12 + Z_delta[x_ori]*Q12x[-1]
    A13 = A13 + Z_delta[x_ori]*Q13x[-1]
    A22 = A22 + Z_delta[x_ori]*Q22x[-1]
    A23 = A23 + Z_delta[x_ori]*Q23x[-1]
    A33 = A33 + Z_delta[x_ori]*Q33x[-1]
    A44 = A44 + Z_delta[x_ori]*Q44x[-1]
    A55 = A55 + Z_delta[x_ori]*Q55x[-1]
    A66 = A66 + Z_delta[x_ori]*Q66x[-1]
    A16 = A16 + Z_delta[x_ori]*Q16x[-1]
    A26 = A26 + Z_delta[x_ori]*Q26x[-1]
    A36 = A36 + Z_delta[x_ori]*Q36x[-1]
    A45 = A45 + Z_delta[x_ori]*Q45x[-1]

if PlaneStress_2D == 'YES':
    A = np.array([
        [A11,   A12,    A16],
        [A12,   A22,    A26],
        [A16,   A26,    A66]])
else:
    A = np.array([
        [A11,   A12,    A13,      0,      0,    A16],
        [A12,   A22,    A23,      0,      0,    A26],
        [A13,   A23,    A33,      0,      0,    A36],
        [  0,     0,      0,    A44,    A45,      0],
        [  0,     0,      0,    A45,    A55,      0],
        [A16,   A26,    A36,      0,      0,    A66]])

#
TOLER = 1.0e-06
tol4 = 1.0e-05
#
if SX == 'e11':
    SXlab = 'ε₁₁'
    SSXlab = 'σ₁₁'
    SX_index = 0
elif SX == 'e22':
    SXlab = 'ε₂₂'
    SSXlab = 'σ₂₂'
    SX_index = 1
elif SX == 'e33':
    SXlab = 'ε₃₃'
    SSXlab = 'σ₃₃'
    SX_index = 2
elif SX == 'e13':
    SXlab = 'ε₁₃'
    SSXlab = 'σ₁₃'
    SX_index = 3
elif SX == 'e23':
    SXlab = 'ε₂₃'
    SSXlab = 'σ₂₃'
    SX_index = 4
elif SX == 'e12':
    SXlab = 'ε₁₂'
    SSXlab = 'σ₁₂'
    if PlaneStress_2D == 'YES':
        SX_index = 2
    else:
        SX_index = 5

if SY == 'e11':
    SYlab = 'ε₁₁'
    SSYlab = 'σ₁₁'
    SY_index = 0
elif SY == 'e22':
    SYlab = 'ε₂₂'
    SSYlab = 'σ₂₂'
    SY_index = 1
elif SY == 'e33':
    SYlab = 'ε₃₃'
    SSYlab = 'σ₃₃'
    SY_index = 2
elif SY == 'e13':
    SYlab = 'ε₁₃'
    SSYlab = 'σ₁₃'
    SY_index = 3
elif SY == 'e23':
    SYlab = 'ε₂₃'
    SSYlab = 'σ₂₃'
    SY_index = 4
elif SY == 'e12':
    SYlab = 'ε₁₂'
    SSYlab = 'σ₁₂'
    if PlaneStress_2D == 'YES':
        SY_index = 2
    else:
        SY_index = 5

# ------------------ Failure loci computation ------------------------ #
start = time.time()
thetaR = 0.0

if OmniLPF == 'YES':
    RMIN = 0
    RMAX = 0.35
else:
    RMIN = 0
    RMAX = 0.30

IMAX = 5000

S_out = []
SX_out = []
SY_out = []
IMAXfail = []

thetaRR = 90
ThetaS = list(range(0, thetaRR, deltaTheta))
ThetaS.append(thetaRR)
for thetaRdeg in ThetaS:
    S_out.append([])
    print('failure envelope for ply orientation = ',thetaRdeg, 'deg ...')
    thetaR = thetaRdeg * DEG2RAD
    #
    c_theta = cos(thetaR)
    s_theta = sin(thetaR)
    c_2theta = cos(2*thetaR)
    s_2theta = sin(2*thetaR)
    c_theta_2 = c_theta*c_theta
    s_theta_2 = s_theta*s_theta
    T1inv = np.array([[c_theta_2, s_theta_2, 0, 0, 0, -s_2theta],
                      [s_theta_2, c_theta_2, 0, 0, 0, s_2theta],
                      [0, 0, 1.0, 0, 0, 1.0],
                      [0, 0, 0, c_theta, -s_theta, 1.0],
                      [0, 0, 0, s_theta, c_theta, 1.0],
                      [0.5*s_2theta, -0.5*s_2theta, 0, 0, 0, c_2theta]])
    #
    T2 = np.array([[c_theta_2, s_theta_2, 0, 0, 0, 0.5*s_2theta],
                   [s_theta_2, c_theta_2, 0, 0, 0, -0.5*s_2theta],
                   [0, 0, 1, 0, 0, 1],
                   [0, 0, 0, c_theta, -s_theta, 1],
                   [0, 0, 0, s_theta, c_theta, 1],
                   [-s_2theta, s_2theta, 0, 0, 0, c_2theta]])
    #
    ENNE = 3
    ENNE_deg = ENNE * 360
    GAMMAS = list(range(0, ENNE_deg, 1))
    GAMMAS.append(ENNE_deg)
    for gamma_i in GAMMAS:
        #print('gamma_i=', gamma_i)
        GAMMA_DEG = gamma_i / ENNE
        GAMM = GAMMA_DEG * DEG2RAD

        DELTA_F = tol4 + 1.0
        iterI = 0

        a = RMIN
        DELTA_F_a, _, _, _, _, _, _ = get_delta_F(a, GAMM, SX, SY,
                                TOLER, DEG2RAD, a1, a2, a32T, a3T,
                                a32C, a3C, PHIC, XT, XC,
                                YC, SL, ST, beta, G12, nu12, T2, Kstif)
        b = RMAX
        DELTA_F_b, _, _, _, _, _, _ = get_delta_F(b, GAMM, SX, SY,
                                TOLER, DEG2RAD, a1, a2, a32T, a3T, 
                                a32C, a3C, PHIC, XT, XC,
                                YC, SL, ST, beta, G12, nu12, T2, Kstif)

        while np.abs(DELTA_F) > tol4 and iterI < IMAX:

            # compute the value for the mean interval value
            c = (a + b) / 2.0
            DELTA_F, FLAG, THETADEG, PHI_D, PSI_D, SXX, SYY = get_delta_F(c, 
                                GAMM, SX, SY, TOLER, DEG2RAD, a1, a2, 
                                a32T, a3T, a32C, a3C, PHIC, XT, XC,
                                YC, SL, ST, beta, G12, nu12, T2, Kstif)

            # get new bounds
            if DELTA_F >= 0:
                a = c
            else:
                b = c
            iterI += 1
        #
        if iterI < IMAX:

            if gamma_i>2:
                diff1 =abs(SYY)-np.abs(SY_out[-1])
                diff2 =abs(SXX)-np.abs(SX_out[-1])
                if diff1> 0.5 or diff2> 0.5:
                    S_out[-1].append([(SX_out[-1][0]),(SY_out[-1][0]),
                                                    FLAG, PHI_D, PSI_D])
                    #
                else:
                    S_out[-1].append([SXX, SYY, FLAG, PHI_D, PSI_D])
                    SX_out.append([SXX])
                    SY_out.append([SYY])
            else:
                S_out[-1].append([SXX, SYY, FLAG, PHI_D, PSI_D])
                SY_out.append([SYY])
                SX_out.append([SXX])
        else:
            print('Maximum number of iterations reached!')
            print('IMAX=', iterI)
            IMAXfail.append([SXX, SYY, FLAG, PHI_D, PSI_D])
            SY_out.append([SYY])
            SX_out.append([SXX])
            S_out[-1].append([(SX_out[-1][0]),(SY_out[-1][0]),
                                                    FLAG, PHI_D, PSI_D])
#
#----------------------------- PLOT -----------------------------------#

if PlaneStress_2D == 'YES':
    load_label = 'InPlaneLoading'
    load_label_long = 'in-plane loading conditions'
else:
    load_label = '3DLoading'
    load_label_long = 'general 3D loading conditions'

SXX = []
SYY = []
FLAGS = []
PHIS = []
PSIS = []

for SS_out in S_out:
    SXX.append([value[0] for value in SS_out])
    SYY.append([value[1] for value in SS_out])
    FLAGS.append([value[2] for value in SS_out])
    PHIS.append([round(value[3],1) for value in SS_out])
    PSIS.append([round(value[4],1) for value in SS_out])


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
    for j in range(len(S_out)):
        vecXX[j] = SXX[j][i]
        vecYY[j] = SYY[j][i]
        vecFLG[j] = FLAGS[j][i]
        vecPHI[j] = PHIS[j][i]
        vecPSI[j] = PSIS[j][i]
    # argmin returns the indices of the min value along an axis
    iminXX = np.argmin(abs(vecXX))
    iminYY = np.argmin(abs(vecYY))
    SXXmin.append(vecXX[iminXX])
    SYYmin.append(vecYY[iminYY])
    FLAGcrit.append(vecFLG[np.minimum(iminXX,iminYY)])
    PHIcrit.append(vecPHI[np.minimum(iminXX,iminYY)])
    PSIcrit.append(vecPSI[np.minimum(iminXX,iminYY)])

SSXXmin = []
SSYYmin = []

for i in range(len(SXXmin)):
    if PlaneStress_2D == 'YES':
        S_vec = np.zeros(3)
    else:
        S_vec = np.zeros(6)
    #
    S_vec[SX_index]=SXXmin[i]
    S_vec[SY_index]=SYYmin[i]
    SS_vec = np.matmul(A, S_vec)
    if PlaneStress_2D == 'YES':
        s11 = SS_vec[0]
        s22 = SS_vec[1]
        s12 = SS_vec[2]
        SSXXmin.append(SS_vec[SX_index])
        SSYYmin.append(SS_vec[SY_index])
    else:   
        s11 = SS_vec[0]
        s22 = SS_vec[1]
        s33 = SS_vec[2]
        s23 = SS_vec[3]
        s13 = SS_vec[4]
        s12 = SS_vec[5]
        SSXXmin.append(SS_vec[SX_index])
        SSYYmin.append(SS_vec[SY_index])


colours = ['r', 'g', 'b', 'yellowgreen', 'm', 'c', 'navy', 'peru',
           'darkslateblue', 'goldenrod']

FLAGmarker = []
FLAGcolour = []
#
for x in FLAGcrit:
    if x ==  1:
        # MATRIX-DOMINATED FAILURE (Transverse cracking) when I3>0
        FLAGmarker.append('circle')# Type of markers: 
        #                          automatic, none, square, diamond
        #                          triangle, x, star, short_dash, 
        #                          long_dash, circle, plus
        FLAGcolour.append('green')# List of colours: 
        #                          black, blue, brown, cyan, gray, green,
        #                          lime, magenta, navy, orange, pink,
        #                          purple, red,  silver, white, yellow
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
                label='[' + str(i * deltaTheta) + ']')
    plt.plot(SXXmin, SYYmin, color='k', linewidth=3.0, label=' '+Omni_label)

    ax.legend(loc='upper center', bbox_to_anchor=(1.25, 0.75))

    if SAVEFIG == 'YES':
        fig.savefig('_results/' +
                'Figure_'+Omni_label+'_Inv3DFC_'+SX+SY+'_'+mat_label+
                '.png', dpi=300, bbox_inches="tight")
    else:
        pass
        
    plt.show()

#------------------- EXCEL Postprocessing --------------------#

if ExcelOutput == 'YES':
    workbook   = xlsxwriter.Workbook('_results/' +
        Omni_label + '_' + SX + SY + '_Inv3DFC_' + load_label +
        '_' + mat_label + '_test.xlsx')
    #
    worksheet1 = workbook.add_worksheet('Overview')
    worksheet1.write_string(1, 1, 'BRIEF OVERVIEW OF THE WORKBOOK:',)
    worksheet1.write_string(2, 1, '-In this workbook you will find the '+
                            Omni_label + ' envelope for '+
                            mat_label + '.',)
    worksheet1.write_string(3, 1, '-Under the assumption of "'+ load_label_long + '"' 
                            ' , load was applied in the following axis:  '+
                            SXlab + ' - ' + SYlab,)
    worksheet1.write_string(4, 1, '-The laminate layup considered for the '+
                            ' the generation of ' + Omni_label +
                            ' in stress space is: ' + str(Layup),)
    #
    worksheetName = 'Omni-strain failure data'
    worksheet2 = workbook.add_worksheet(worksheetName)
    chart =  workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
    #
    cell_format = workbook.add_format()
    cell_format.set_bottom()
    cell_format.set_bold()
    cell_format.set_font_color('navy')
    header_format = workbook.add_format({'bold': True,
                'font_name': 'Calibri',
                'font_size':12,
                'align': 'center',
                'valign': 'vcenter',
                'fg_color': '#D7E4BC',
                'border': 1,
                'bottom': 1,
                'border_color': 'black',
                'bottom_color': 'black', 
                                        })
    header_format2 = workbook.add_format({'bold': True,
                'font_name': 'Calibri',
                'font_size':11, 
                                })
    #
    cell_format2 = workbook.add_format()
    cell_format2.set_right()

    for i in range(len(SXX)):
        worksheet2.write_string(2, int(4*i), SXlab, header_format)
        worksheet2.write_column(3, int(4*i), SXX[i])
        worksheet2.write_string(2, int(1+4*i), SYlab, header_format)
        worksheet2.write_column(3, int(1+4*i), SYY[i])
        worksheet2.write_string(2, int(2+4*i), 'ϕ (deg)', header_format)
        worksheet2.data_validation(2, int(2+4*i),2, int(2+4*i),{
            'validate': 'any','input_message': 'Kinking angle (deg)',})
        worksheet2.write_column(3, int(2+4*i), PHIS[i])
        worksheet2.write_string(2, int(3+4*i), 'ψ (deg)', header_format)
        worksheet2.data_validation(2, int(3+4*i),2, int(3+4*i),{
            'validate': 'any','input_message': 
            'Angle of the kinking plane (deg)',})
        worksheet2.write_column(3, int(3+4*i), PSIS[i], cell_format2)
        worksheet2.write_string(1, int(4*i), '['+str(i*deltaTheta) + ']', 
                                header_format)
        chart.add_series({
            'name':       [worksheetName, 1, int(4*i)]  ,
            'categories': [worksheetName, 3, int(4*i), int(ENNE_deg+3), int(4*i)],
            'values':     [worksheetName, 3, int(1+4*i), int(ENNE_deg+3), int(1+4*i)],
            'line':       {'width': 1.20,},
            })
    #
    worksheet2.write_string(2, int(4*len(SXX)+1), SXlab, header_format)
    worksheet2.write_column(3, int(4*len(SXX)+1), SXXmin)
    worksheet2.write_string(2, int(4*len(SXX)+2), SYlab, header_format)
    worksheet2.write_column(3, int(4*len(SXX)+2), SYYmin)
    worksheet2.write_string(1, int(4*len(SXX)+1), Omni_label + ' in strain space', header_format2)
    worksheet2.write_string(2, int(4*len(SXX)+3), 'Flags', header_format)
    worksheet2.write_column(3, int(4*len(SXX)+3), FLAGcrit)
    worksheet2.write_column(3, int(4*len(SXX)+4), FLAGmarker)
    worksheet2.write_column(3, int(4*len(SXX)+5), FLAGcolour)
    worksheet2.write_string(2, int(4*len(SXX)+7), SSXlab, header_format)
    worksheet2.write_column(3, int(4*len(SXX)+7), SSXXmin)
    worksheet2.write_string(2, int(4*len(SXX)+8), SSYlab, header_format)
    worksheet2.write_column(3, int(4*len(SXX)+8), SSYYmin)
    worksheet2.write_string(1, int(4*len(SXX)+7), Omni_label+ ' in stress space ', header_format2)
    worksheet2.write_string(0, int(4*len(SXX)+7), 'Layup : '+ str(Layup), header_format2)
    worksheet2.write_column(3, int(4*len(SXX)+9), PHIcrit)
    worksheet2.write_string(2, int(4*len(SXX)+9), 'ϕ (deg)', header_format)
    worksheet2.data_validation(2, int(4*len(SXX)+9),2, int(4*len(SXX)+9),{
            'validate': 'any','input_message': 'Kinking angle (deg)',})
    worksheet2.write_column(3, int(4*len(SXX)+10), PSIcrit)
    worksheet2.write_string(2, int(4*len(SXX)+10), 'ψ (deg)', header_format)
    worksheet2.data_validation(2, int(4*len(SXX)+10),2, int(4*len(SXX)+10),{
            'validate': 'any','input_message': 
            'Angle of the kinking plane (deg)',})
    #
    ##########################################################################
    #
    worksheet3 = workbook.add_worksheet('Charts')        

    chart.add_series({
        'name':       [worksheetName, 1, int(4*len(SXX)+1)] ,  
        'categories': [worksheetName, 3, int(4*len(SXX)+1), int(ENNE_deg+3), int(4*len(SXX)+1)],
        'values':     [worksheetName, 3, int(4*len(SXX)+2), int(ENNE_deg+3), int(4*len(SXX)+2)],
        'line':       {'color': 'black', 'dash_type':'dash'},
        #dash_type can be: 'solid','round_dot','square_dot','dash','dash_dot'
        #                  'long_dash','long_dash_dot','long_dash_dot_dot'
        })
    
    chart.set_x_axis({'name': SXlab,
                'name_font': {
                    'name': 'Segoe UI',
                    'size': 14,
                    'bold': True},
    #                      'label_position': 'low',
                'num_font':  {'name': 'Segoe UI','size':12 },
                'line': {'color': 'black','width': 1.00},
                'major_gridlines': {
                    'visible': True,
                    'line': {'color': '#D9D9D9','width': 0.75, 
                            'dash_type': 'solid'}},
                'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                'minor_tick_mark': 'none'})
    chart.set_y_axis({'name': SYlab,
                'name_font': {
                    'name': 'Segoe UI',
                    'size': 14,
                    'bold': True},
                'name_layout':{
                            'x':      0.0,
                            'y':      0.45},
                'num_font':  {'name': 'Segoe UI','size':12 },
                'line': {'color': 'black','width': 1.00},
                'major_gridlines': {
                        'visible': True,
                        'line': {'color': '#D9D9D9','width': 0.75, 
                                'dash_type': 'solid'}},
                'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                'minor_tick_mark': 'none'})
    chart.set_title({'name': mat_label,
                    'name_font': {
                        'name': 'Segoe UI',
                        'size': 14,
                        'bold': True}},)
    chart.set_plotarea({
                        'layout': {
                            'x':      0.08,
                            'y':      0.13,
                            'width':  0.65,
                            'height': 0.75,}})
    chart.set_legend({'font': {'name': 'Segoe UI','size': 10.5},
                    'layout': {'x':      0.76,
                               'y':      0.12,
                               'width':  0.24,
                               'height': 0.80,}})
    worksheet3.insert_chart(3+18*0, 3, chart,{'x_scale': 1.194, 'y_scale': 1.227} )
    ######### CHART 1 ##########
    chart1 =  workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
    #scatter has the following subtypes: straight_with_markers, straight
    #                                    smooth_with_markers, smooth
    chart1.add_series({
        'name':       [worksheetName, 1, int(4*len(SXX)+7)] ,  
        'categories': [worksheetName, 3, int(4*len(SXX)+7), int(ENNE_deg+3), int(4*len(SXX)+7)],
        'values':     [worksheetName, 3, int(4*len(SXX)+8), int(ENNE_deg+3), int(4*len(SXX)+8)],
        'line':       {'color': 'black'},
        })
    chart1.set_x_axis({'name': SSXlab + ' [MPa]',
              'name_font': {
                    'name': 'Segoe UI',
                    'size': 14,
                    'bold': True},
              'num_font':  {'name': 'Segoe UI','size':12 },
              'line': {'color': 'black','width': 1.00},
              'major_gridlines': {
                    'visible': True,
                    'line': {'color': '#D9D9D9','width': 0.75, 
                            'dash_type': 'solid'}},
              'major_tick_mark': 'none', # choose between: none, inside, outside, cross
              'minor_tick_mark': 'none'})
    chart1.set_y_axis({'name': SSYlab  + ' [MPa]',
                'name_font': {
                    'name': 'Segoe UI',
                    'size': 14,
                    'bold': True},
                'name_layout':{
                            'x':      0.0,
                            'y':      0.40},
                'num_font':  {'name': 'Segoe UI','size':12 },
                'line': {'color': 'black','width': 1.00},
                'major_gridlines': {
                        'visible': True,
                        'line': {'color': '#D9D9D9','width': 0.75, 
                                'dash_type': 'solid'}},
                'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                'minor_tick_mark': 'none'})
    chart1.set_title({'name': str(Layup)+' '+ mat_label,
                    'name_font': {
                        'name': 'Segoe UI',
                        'size': 14,
                        'bold': True}},)
    chart1.set_legend({'font': {'name': 'Segoe UI','size': 10.5}, 
                    'position': 'right'})
    worksheet3.insert_chart(3+18*0, 13, chart1,{'x_scale': 1.194, 'y_scale': 1.227} )

    ######### CHART 2 ##########
    worksheet3.write_string(6+18*1, 11, 'LEGEND:')
    cell_format = workbook.add_format({'bold': True, 'font_color': 'cyan'})
    worksheet3.write_string(7+18*1, 11, 'triangle(cyan): Fibre tensile failure',cell_format)
    cell_format = workbook.add_format({'bold': True, 'font_color': 'green'})
    worksheet3.write_string(8+18*1, 11, 'x(green): Fibre kinking when I3>0',cell_format)
    worksheet3.write_string(10+18*1, 11, 'o(green): Matrix cracking when I3>0',cell_format)
    cell_format = workbook.add_format({'bold': True, 'font_color': 'red'})
    worksheet3.write_string(9+18*1, 11, 'x(red): Fibre kinking when I3<0',cell_format)
    worksheet3.write_string(11+18*1, 11, 'o(red): Matrix cracking when I3<0',cell_format)
    chart2 =  workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
    chart2.add_series({
        'name':       [worksheetName, 1, int(4*len(SXX)+1)] ,  
        'categories': [worksheetName, 3, int(4*len(SXX)+1), int(ENNE_deg+3), int(4*len(SXX)+1)],
        'values':     [worksheetName, 3, int(4*len(SXX)+2), int(ENNE_deg+3), int(4*len(SXX)+2)],
        'line':       {'color': 'black','width': 1.00},
        })
    delta_value = 4*ENNE
    for i in range(0,len(FLAGmarker),delta_value):
        if FLAGmarker[i] == 'x':
            chart2.add_series({
                'categories': [worksheetName, 3+i, int(4*len(SXX)+1), int(i+3), int(4*len(SXX)+1)],
                'values':     [worksheetName, 3+i, int(4*len(SXX)+2), int(i+3), int(4*len(SXX)+2)],
                'marker':     {
                     'type': FLAGmarker[i],
                     'size': 5,
                     'border': {'color': FLAGcolour[i]},
                     'fill':   {'none': True},},})
        elif FLAGmarker[i] == 'triangle':
            chart2.add_series({
                'categories': [worksheetName, 3+i, int(4*len(SXX)+1), int(i+3), int(4*len(SXX)+1)],
                'values':     [worksheetName, 3+i, int(4*len(SXX)+2), int(i+3), int(4*len(SXX)+2)],
                'marker':     {
                     'type': FLAGmarker[i],
                     'size': 5,
                     'border': {'none': True},
                     'fill':   {'color': FLAGcolour[i]},},})
        else:
            chart2.add_series({
                'categories': [worksheetName, 3+i, int(4*len(SXX)+1), int(i+3), int(4*len(SXX)+1)],
                'values':     [worksheetName, 3+i, int(4*len(SXX)+2), int(i+3), int(4*len(SXX)+2)],
                'marker':     {
                     'type': FLAGmarker[i],
                     'size': 5,
                     'border': {'color': FLAGcolour[i]},
                     'fill':   {'none': True},},})
                     #'border': {'none': True},
                     #'fill':   {'color': FLAGcolour[i]},},})
    chart2.set_x_axis({'name': SXlab,
                'name_font': {
                    'name': 'Segoe UI',
                    'size': 14,
                    'bold': True},
                'num_font':  {'name': 'Segoe UI','size':12 },
                'line': {'color': 'black','width': 1.00},
                'major_gridlines': {
                    'visible': True,
                    'line': {'color': '#D9D9D9','width': 0.75, 
                            'dash_type': 'solid'}},
                'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                'minor_tick_mark': 'none'})
    chart2.set_y_axis({'name': SYlab,
                'name_font': {
                        'name': 'Segoe UI',
                        'size': 14,
                        'bold': True},
                'name_layout':{
                            'x':      0.0,
                            'y':      0.45},
                'num_font':  {'name': 'Segoe UI','size':12 },
                'line': {'color': 'black','width': 1.00},
                'major_gridlines': {
                        'visible': True,
                        'line': {'color': '#D9D9D9','width': 0.75, 
                                'dash_type': 'solid'}},
                'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                'minor_tick_mark': 'none'})
    chart2.set_title({'none': True})
    chart2.set_legend({'none': True})
    worksheet3.insert_chart(5+18*1, 3, chart2,{'x_scale': 1.058, 'y_scale': 1.227} )

    ######### CHART 3 ##########
    chart3 =  workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
    delta_value=3*ENNE
    for i in range(0,len(FLAGmarker),delta_value):
        chart3.add_series({
            'categories': [worksheetName, 3+i, int(4*len(SXX)+1), int(i+3+delta_value), int(4*len(SXX)+1)],
            'values':     [worksheetName, 3+i, int(4*len(SXX)+2), int(i+3+delta_value), int(4*len(SXX)+2)],
            'line':       {'color': FLAGcolour[i],'width': 1.00},
            })
    chart3.set_x_axis({'name': SXlab,
                'name_font': {
                    'name': 'Segoe UI',
                    'size': 14,
                    'bold': True},
                'num_font':  {'name': 'Segoe UI','size':12 },
                'line': {'color': 'black','width': 1.00},
                'major_gridlines': {
                    'visible': True,
                    'line': {'color': '#D9D9D9','width': 0.75, 
                            'dash_type': 'solid'}},
                'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                'minor_tick_mark': 'none'})
    chart3.set_y_axis({'name': SYlab,
                'name_font': {
                    'name': 'Segoe UI',
                    'size': 14,
                    'bold': True},
                'name_layout':{
                            'x':      0.0,
                            'y':      0.45},
                'num_font':  {'name': 'Segoe UI','size':12 },
                'line': {'color': 'black','width': 1.00},
                'major_gridlines': {
                        'visible': True,
                        'line': {'color': '#D9D9D9','width': 0.75, 
                                'dash_type': 'solid'}},
                'major_tick_mark': 'none', # choose between: none, inside, outside, cross
                'minor_tick_mark': 'none'})
    chart3.set_title({'none': True})
    chart3.set_legend({'none': True})
    worksheet3.insert_chart(5+18*2, 3, chart3,{'x_scale': 1.058, 'y_scale': 1.227} )
    workbook.close()
else:
    pass

elapsed = (time.time() - start)
print('elapsed time=',round(elapsed,2),'sec')
print('END')