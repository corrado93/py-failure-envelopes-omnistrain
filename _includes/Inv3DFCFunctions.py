################################################################################
#                                                                              #
#      3D INVARIANT-BASED FAILURE CRITERIA FOR LAMINATED COMPOSITES            #
#                                                                              #
################################################################################
#                                                                              #
#   Implemented by                                                             #
#   Giuseppe Corrado, PhD candidate FEUP (Porto) - gcorrado@fe.up.pt           #
#                     and ITN-EID Fellow (https://optimacs.net/)               #
#                                                                              #
################################################################################
#
#   Main references:
#
#   - P.P. Camanho et al. /Int. Journal of Solids and Structures 55(2015) 92–107
#
# ---------------------------------------------------------------------------- #
from numpy import *
import math
import pickle
import numpy as np
from matplotlib import pyplot as plt

def FailureMatrix(s22,s33,s12,s13,s23,TOLER,DEG2RAD,YC,SL,ST,
                  a1,a2,a32T,a3T,a32C,a3C,PHIC):
    
# Initialize variables
    FM_MAX  = 0.0
    FLAG_M  = 0.0
    
############### Matrix failure ################

    I1 = 0.25*s22*s22 - 0.5*s22*s33 + \
         0.25*s33*s33 + s23*s23
    I2 = s12*s12 + s13*s13
    I3 = s22 + s33
    
# TRANSVERSELY ISOTROPIC FAILURE CRITERIA (F):
    if I3 <= 0.0 :
       F = a1*I1 + a2*I2 + a3C*I3 + a32C*I3*I3 - 1.0
       FM_MAX = F + 1.0
       FLAG_M = -1
#       print('FM_MAX_c=',FM_MAX)
#       print('s22=',s22)
#       print('s12=',s12)
    else:
       F = a1*I1 + a2*I2 + a3T*I3 + a32T*I3*I3 - 1.0
       FM_MAX = F + 1.0
       FLAG_M = 1
#       print('FM_MAX_t=',FM_MAX)
#       print('s22=',s22)
#       print('s12=',s12)
    
    
    if FM_MAX >= (1.0-TOLER):
       THETADEG=Fractangle(s22,s33,s12,s13,s23,TOLER,DEG2RAD,I1,I2,I3,YC,ST,
                            a1,a2,a32T,a3T,a32C,a3C,PHIC)
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
       #FF_MAX = s11/XT # Max stress failure criterion     
       FLAG_F = 2
    else:            # fiber compression
       if abs(s12) < TOLER and abs(s13) < TOLER:
#          if abs(s23) < TOLER:
#            psi = 0.5*pi
#          else:
            #psi = 0.5*np.arctan2((2.0*s23),(s22 - s33))
          psi = 0.5*math.atan2(2.0*s23,(s22 - s33))
       else:
          #psi = np.arctan2(s13,s12)
          psi = math.atan2(s13,s12)
#       print('psi=',psi)
       PSI_D  = psi / DEG2RAD
#       print('psi_deg=',PSI_D)
       # Find the stresses on the kinking plane: psi rotation
       s11psi = s11
       s22psi = s22*cos(psi)*cos(psi) + 2.0*s23*sin(psi)*cos(psi) \
                 + s33*sin(psi)*sin(psi)
       s33psi = s22*sin(psi)*sin(psi) - 2.0*s23*sin(psi)*cos(psi) \
                 + s33*cos(psi)*cos(psi)
       s12psi = s12*cos(psi) + s13*sin(psi)
       s13psi = -s12*sin(psi) + s13*cos(psi)
       s23psi = (s33 - s22)*sin(psi)*cos(psi) + s23*(cos(psi)*cos(psi) - sin(psi)*sin(psi))
#       Spsi= np.array([
#                [s11psi,s12psi,s13psi], [s12psi,s22psi,s23psi], 
#                [s13psi,s23psi,s33psi]])
       #print('Spsi=',Spsi)
       # Calculation of gamma_m
       #  materialLaw = 0 (Newton-Raphson with chi)
       #  materialLaw = 1 (Simple linear behavior)
       #  materialLaw = 2 (Newton-Raphson with phi0)
       #  materialLaw = 3 (Bisection with chi)
       materialLaw = 3
       
       if materialLaw == 0:
       	  maxiter = 100
          tol = 1e-07
          gamma_m = 15*PHIC
          gammamToler = tol * PHIC
          error = 1
          # Calculation of gamma_m
          chi = sin(PHIC+PHIC)*(abs(XC))/(PHIC+PHIC)
          iter = 0  
          while abs(error) > gammamToler and iter <= maxiter :       
               F = gamma_m*chi + 0.5*sin(2.0*gamma_m)*(s11psi - s22psi) \
                  - abs(s12psi)*cos(2.0*gamma_m)
               dF = chi + cos(2.0*gamma_m)*(s11psi - s22psi) + \
                  2.0*abs(s12psi)*sin(2.0*gamma_m)
               delta = -F/dF
               error = delta
               gamma_m = max(gamma_m+delta,0)
               iter = iter + 1
               if iter >= maxiter:
               		print('kInv3DFC_FK Newton method! gammam =',gamma_m)
               		#gamma_m = abs(s12psi)/G12
               		#print('NEW gammam =',gamma_m)  
          # Calculation of phi
          if s12psi >= 0.0:
               PHI = abs(gamma_m)
          else:
               PHI = -abs(gamma_m)
          #if s22 >(1.5*ST) :
               #PHI = 0
               #print('s22>ST')
       elif materialLaw == 1:
          #gamma_m = phiC*s11/XC  
          # Calculate gamma_mc
          gamma_mc = sin(2.0*PHIC)*(abs(XC))/2.0/G12
          # Define phi_0
          phi_0_theta = PHIC - gamma_mc
          # Compute gamma_m (linear shear behavior)
          gamma_m = (phi_0_theta*G12 + abs(s12psi))/(G12 + s11psi - \
                    s22psi) - phi_0_theta
          # Calculate phi
          if s12psi >= 0.0:
             PHI = phi_0_theta + gamma_m
          else:
             PHI = -(phi_0_theta + gamma_m)
       elif materialLaw == 2:
          tol = PHIC*TOLER
          maxiter = 1000
          iter = 0
          error = PHIC
          PHI0 = PHIC*5.0
          while np.abs(error) > tol and iter < maxiter  : 
             Fphi0 = PHIC - PHI0 - np.abs(XC*sin(2*PHI0)/(2*G12)+ \
                      beta*XC*XC*XC*sin(2*PHI0)*sin(2*PHI0)*sin(2*PHI0)/8)
             dFphi0 = -1-np.abs(XC*cos(2*PHI0)/(G12)+ \
                      3*beta*XC*XC*XC*sin(2*PHI0)*sin(2*PHI0)*cos(2*PHI0)/4)
             PHI0_plus = PHI0 - Fphi0/dFphi0
             error = np.abs(Fphi0/dFphi0)
             PHI0 = PHI0_plus
             iter = iter + 1 
             if iter >= maxiter:
                  print('kInv3DFC_FK Newton method! gammam =',gamma_m)
                  gamma_m = (abs(XC))*abs(s12psi)/G12
                  print('NEW gammam =',gamma_m) 
          s12R = 0.5*sin(2*PHI0)*(s22*cos(psi)*cos(psi)-s11+s33*sin(psi)*sin(psi)+ \
                s23*sin(2*psi))+cos(2*PHI0)*(s12*cos(psi)+s13*sin(psi))
          PHI = np.sign(s12R)*(PHI0+np.abs(s12R/G12+beta*s12R*s12R*s12R))
       
       elif materialLaw == 3:
          maxiter = 70
          tol = 1e-07
          gammamToler = tol * PHIC
          error = 1
          # Calculation of gamma_m
          chi = sin(PHIC+PHIC)*(abs(XC))/(PHIC+PHIC)
          a = 10E-05
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
                  #gamma_m = abs(s12psi)/G12
                  #print('NEW gammam =',gamma_m)
          # Calculation of phi
          if s12psi >= 0.0:
               PHI = abs(gamma_m)
          else:
               PHI = -abs(gamma_m)

      
       # Preferred direction (A):
       A = array([
         [cos(PHI)], [cos(psi)*sin(PHI)], [sin(psi)*sin(PHI)]])
#       print('A=',A)
       # Structural tensor (AA)
       N = 2
       AT = np.transpose(A)
       
       AA = np.matmul(A,AT)
       
       #Stress tensor S
       
       S = np.array([
                [s11,s12,s13], [s12,s22,s23], [s13,s23,s33]
                ])
#       print('S=',S)
       #Reaction stress tensor SR
       TRS = np.trace(S)
       ATS = np.matmul(AT,S)
       
       ATSA = np.matmul(ATS,A)
       SR1 = np.identity(3)*(TRS - ATSA[(0,0)])  
       #print('ATSA[(0,0)]=',ATSA[(0,0)])
       SR2=(TRS -3.0*ATSA[(0,0)])*AA
       SR = 0.5*SR1 - 0.5*SR2       

       # Plasticity inducing stresses (SPL)
       SPL = S - SR       
       # Transversely isotropic invariants (I's)
       
       SPL2 = np.matmul(SPL,SPL)
       TRSPL2 = np.trace(SPL2)
       #TRSPL2 = SPL2[(0,0)] + SPL2[(1,1)] + SPL2[(2,2)]
       ATSPL2 = np.matmul(AT,SPL2)
       ATSPL2A = np.matmul(ATSPL2,A)
       #
       I1 = 0.5 * TRSPL2 - ATSPL2A[(0,0)]
       I2 = ATSPL2A[(0,0)]
       I3 = TRS - ATSA[(0,0)]
#       print('I1=',I1)
#       print('I2=',I2)
#       print('I3=',I3)
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