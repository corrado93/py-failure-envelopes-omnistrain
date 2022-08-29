########################################################################
#                                                                      #
#     3D INVARIANT-BASED FAILURE CRITERIA FOR LAMINATED COMPOSITES     #
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
#
# -------------------------------------------------------------------- #
from numpy import *
import math
import numpy as np

def FailureMatrix(s22,s33,s12,s13,s23,TOLER,DEG2RAD,YC,SL,ST,
                  a1,a2,a32T,a3T,a32C,a3C,PHIC):
    
# Initialize variables
    FM_MAX  = 0.0
    FLAG_M  = 0.0
    
# Computing the stress invariants 

    I1 = 0.25*s22*s22 - 0.5*s22*s33 + 0.25*s33*s33 + s23*s23
    I2 = s12*s12 + s13*s13
    I3 = s22 + s33
    
# TRANSVERSELY ISOTROPIC FAILURE CRITERIA (F):
    if I3 <= 0.0 :
       F = a1*I1 + a2*I2 + a3C*I3 + a32C*I3*I3 - 1.0
       FM_MAX = F + 1.0
       FLAG_M = -1
    else:
       F = a1*I1 + a2*I2 + a3T*I3 + a32T*I3*I3 - 1.0
       FM_MAX = F + 1.0
       FLAG_M = 1
    
    if FM_MAX >= (1.0-TOLER):
       THETADEG=Fractangle(s22,s33,s12,s13,s23,TOLER,DEG2RAD,I1,I2,I3,
                            YC,ST,a1,a2,a32T,a3T,a32C,a3C,PHIC)
    else :
       THETADEG=None

    return FM_MAX, FLAG_M, THETADEG

def Fx(x,chi,s11psi,s22psi,s12psi):
    sin2x=sin(2*x)
    cos2x=cos(2*x)
    Fxx = x*chi + 0.5*sin2x*(s11psi - s22psi) - abs(s12psi)*cos2x
    return Fxx
    

def FailureFiber(s11,s22,s33,s12,s13,s23,nu12,TOLER,DEG2RAD,
                 a1,a2,a32T,a3T,a32C,a3C,PHIC,XT,XC,YC,SL,ST,beta,G12):
    # Initialize variables
    FF_MAX = 0.0
    PHI_D = 0.0
    PSI_D = 0.0
    FLAG_F = 0.0
    THETADEG = 0.0

    if s11 >= 0.0:   # fiber tension (maximum strain criterion)

       FF_MAX = (s11 - nu12*(s22+s33))/XT   
       FLAG_F = 2
    else:            # fiber compression
       if abs(s12) < TOLER and abs(s13) < TOLER:
          psi = 0.5*math.atan2(2.0*s23,(s22 - s33))
       else:
          psi = math.atan2(s13,s12)
       PSI_D  = psi / DEG2RAD
       # Find the stresses on the kinking plane: psi rotation
       s11psi = s11
       s22psi = s22*cos(psi)*cos(psi) + 2*s23*sin(psi)*cos(psi) \
                 + s33*sin(psi)*sin(psi)
       s33psi = s22*sin(psi)*sin(psi) - 2*s23*sin(psi)*cos(psi) \
                 + s33*cos(psi)*cos(psi)
       s12psi = s12*cos(psi) + s13*sin(psi)
       s13psi = -s12*sin(psi) + s13*cos(psi)
       s23psi = (s33 - s22)*sin(psi)*cos(psi) + s23*(cos(psi)*cos(psi) \
                 - sin(psi)*sin(psi))
       
       # Calculation of gamma_m
       #  Bisection method with "chi" parameter

       maxiter = 70
       tol = 1e-07
       gammamToler = tol * PHIC
       error = 1
       # Calculation of gamma_m
       chi = sin(PHIC+PHIC)*(abs(XC))/(PHIC+PHIC)
       a = 10e-05
       F_a = Fx(a,chi,s11psi,s22psi,s12psi)
       b = pi
       F_b = Fx(b,chi,s11psi,s22psi,s12psi)
       if (F_a*F_b) > 0:
          gamma_m = 0
       else:    
          iter = 0
          while abs(error) > gammamToler and iter <= maxiter :   
             c = ((a + b) / 2)
             F_c = Fx(c,chi,s11psi,s22psi,s12psi)
             # get new bounds
             if F_a*F_c > 0:
                error = c-a
                a = c
             elif F_a*F_c <= 0:
                error = c-b
                b = c
             gamma_m = max(0,c)
             iter = iter + 1
             if iter > maxiter:
                print('Bisection method! gammam =',gamma_m)
                gamma_m = 0
       # Calculation of phi
       if s12psi >= 0.0:
          PHI = abs(gamma_m)
       else:
          PHI = -abs(gamma_m)
      
       # Preferred direction (A):
       A = array([
         [cos(PHI)], [cos(psi)*sin(PHI)], [sin(psi)*sin(PHI)]])
       # Structural tensor (AA)
       N = 2
       AT = np.transpose(A)
       
       AA = np.matmul(A,AT)
       
       #Stress tensor S
       
       S = np.array([
                [s11,s12,s13], [s12,s22,s23], [s13,s23,s33]
                ])
       #Reaction stress tensor SR
       TRS = np.trace(S)
       ATS = np.matmul(AT,S)
       
       ATSA = np.matmul(ATS,A)
       SR1 = np.identity(3)*(TRS - ATSA[(0,0)])  
       SR2=(TRS -3.0*ATSA[(0,0)])*AA
       SR = 0.5*SR1 - 0.5*SR2       

       # Plasticity inducing stresses (SPL)
       SPL = S - SR       
       # Transversely isotropic invariants (I's)
       
       SPL2 = np.matmul(SPL,SPL)
       TRSPL2 = np.trace(SPL2)
       ATSPL2 = np.matmul(AT,SPL2)
       ATSPL2A = np.matmul(ATSPL2,A)
       #
       I1 = 0.5 * TRSPL2 - ATSPL2A[(0,0)]
       I2 = ATSPL2A[(0,0)]
       I3 = TRS - ATSA[(0,0)]
       # TRANSVERSELY ISOTROPIC FAILURE CRITERIA (F):
       if I3 <= 0.0 :
          F = a1*I1 + a2*I2 + a3C*I3 + \
              a32C*I3*I3 - 1.0
          FF_MAX = F + 1.0
          FLAG_F = -3
       else:
          F = a1*I1 + a2*I2 + a3T*I3 + \
              a32T*I3*I3 - 1.0
          FF_MAX = F + 1.0
          FLAG_F = 3
       
       PHI_D  = PHI / DEG2RAD

       
       
    return FF_MAX, FLAG_F, THETADEG, PHI_D, PSI_D

def Fractangle(s22,s33,s12,s13,s23,TOLER,DEG2RAD,I1,I2,I3,YC,ST,
                a1,a2,a32T,a3T,a32C,a3C,PHIC):    
    # Calculate the angle of fracture plane
    # Note: ST = R23 is assumed for the transverse compression-
    #       dominated failure modes
    theta = 0.0
       
    if I3 < 0.0 :
       if s22 < 0.0 and abs(s22)>= ST :
            if  abs(s33)>abs(s22):
                theta = np.arcsin(np.sqrt(ST/min(abs(s33),abs(YC))))
            else:
                theta = np.arccos(np.sqrt(ST/min(abs(s22),abs(YC))))
                
       elif s33 < 0.0 and abs(s33)>= ST :
            theta = np.arcsin(np.sqrt(ST/min(abs(s33),abs(YC)))) 
       elif abs(s23)<TOLER and abs(s22-s33)<TOLER :
            if abs(s12)<TOLER and abs(s13)<TOLER :   
               theta = np.arccos(np.sqrt(ST/min(abs(s22),abs(YC))))
            else:
               theta = np.arctan2(s13,s12)
       else:
            if abs(a1*I1 + a3C*I3 + a32C*I3*I3)>= abs(a2*I2):
               theta = 0.5*np.arctan2((s23+s23),(s22-s33))
            else:
               theta = np.arctan2(s13,s12)
    else:
       if abs(s23)<TOLER and abs(s22-s33)<TOLER :
            if abs(s12)<TOLER and abs(s13)<TOLER :
               theta = math.pi/4.0
            else:
               theta = np.arctan2(s13,s12)
       else:
            if abs(a1*I1 + a3T*I3 + a32T*I3*I3)>= abs(a2*I2):
               theta = 0.5*np.arctan2((s23+s23),(s22-s33))
            else:
               theta = np.arctan2(s13,s12)
    THETADEG = theta/DEG2RAD
   
    return THETADEG