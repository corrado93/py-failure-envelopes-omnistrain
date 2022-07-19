################################################################################
#                                                                              #
#      Function to compute the DEGRADED MATERIAL PROPERTIES                    #
#      required for the OMNI LA            #
#                                                                              #
################################################################################
#      
#      Implemented by Giuseppe Corrado, PhD candidate
#      
#
################################################################################

from numpy import *
import math
import numpy as np

def get_DegradedProp(E2,nu12,G12,v_f,Em_star,Em,eta_y,eta_s):
    Gm = Em/2.7
    #Computate the degraded material properties
    nu12 = Em_star * nu12
    k_y = eta_y*((1.0-v_f)/v_f)
    E2 = (1.0+k_y)/ ((1.0/E2)+(k_y*(Em_star*Em+E2*(1.0-Em_star)))/(E2*Em_star*Em))
    k_s = eta_s*((1-v_f)/v_f)
    G12 = (1.0+k_s)/((1.0/G12)+(k_s*(Em_star*Gm+G12*(1.0- Em_star)))/(G12*Em_star*Gm))

    return E2, G12, nu12