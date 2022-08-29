# ------------------------- User INPUT ------------------------------- #

OmniLPF = 'NO'         #'YES' or 'NO': 
#                       'YES' for Omni Last Ply Failure (LPF) envelopes,
#                             generated using degraded material proper-
#                             ties.
#                       'NO' for Omni First Ply Failure (FPF) envelopes, 
#                            generated using intact material properties.
#
PlaneStress_2D = 'NO'  #'YES' or 'NO':
#                       'YES' for failure envelopes under plane stress 
#                             conditions (in-plane loading)
#                       'NO' for failure envleopes under general loading 
#                             conditions (3D)   
# 
ExcelOutput = 'YES'    #'YES' or 'NO':
#                       'YES' to create a ".xlsx" file with results and  
#                             charts
#                       'NO' to not create ".xlsx" files
PlotFigure = 'NO'     #'YES' or 'NO':
#                       'YES' to plot the failure envelopes with pyplot 
#                       'NO' to not output any figure
SAVEFIG = 'NO'         #'YES' or 'NO' 
#                       'YES' to save the figure created with pyplot 
#                       'NO' to not save any figure
#     Note: This is a suboption of PlotFigure! If you want to save the 
#           figure, you need to set PlotFigure = 'YES'.
#
## Material properties
#
import material_data
# Choose the material from the database - "material_data.py"
# or add a new one with the same number of variable (or keys) contained
# in the dictionary.
#   materials = IM78552, EglassMY750, custom ...
material = material_data.EglassMY750

# Laminate design 
# (Needed to represent failure envelopes in STRESS SPACE)

Symmetric = 'YES'       #'YES' or 'NO'
# INPUT YOUR LAYUP BELOW [use orientation from -90 to 90]
Layup = [45, 0, -45, 90]
#Layup = [55, -55]
#Layup = [0, 0, 90, 90, 0, 0]

# Omni-strain envelope
# Choose the strain loading on X axis 
# (e11,e22,e33,e12,e13,e23,e12)
SX = 'e11'
# Choose the strain loading on Y axis
#(e11,e22,e33,e12,e13,e23,e12)
SY = 'e22'
#
# INPUT FOR THE OMNI-STRAIN FAILURE ENVELOPE
#
deltaTheta = 15 # Choose a value from 1 to 90[deg] (default = 15) to set
# at which ply orientation you want a failure envelope in strain space. 
# It highly affects the accuracy of the resulting omni strain FPF  
# envelope (e.g. with 15[deg] you get 7 failure envelopes, for the  
# following ply orientations: 0-15-30-45-60-75-90[deg]).
#----------------------------------------------------------------------#
#----------------------------------------------------------------------#
# Material properties will be read from the database
# Please apply modifications to a material directly from the database or
# add your custom material and load it with:
#                              "material = material_data.custom"
#----------------------------------------------------------------------#
#----------------------------------------------------------------------#
#
E1 = material['E1']           # Longitudinal Tension Young Modulus [MPa]
E2 = material['E2']           # Transverse Young Modulus [MPa]
nu12 = material['v12']        # In-plane Poisson's coefficient [-]
nu23 = material['v23']        # Transverse Poisson's coefficient [-]
G12 = material['G12']         # In-plane Shear Modulus [MPa]
v_f = material['v_f']         # Fibre volume fraction [-]
GIC = material['GIC']         # Mode I Fracture Toughness [kJ/m2]
GIIC = material['GIIC']       # Mode II Fracture Toughness [kJ/m2]
KP = material['KP']           # Shear incremental stiffness under 
#                               plastic flow [-]
SLP = material['SLP']         # Shear Stress that activates plastic 
#                               flow [MPa]
eta_y = material['eta_y']     # Stress partitioning parameter (transver-
#                               se). For more info, check Chapter 2 of
#                               "S.W. Tsai et al. (2017). Composite La-
#                               minates: Theory and practice of analy- 
#                               sis, design and automated layup. Stan-
#                               ford Aeronautics & Astronautics."
eta_s = material['eta_s']     # Stress partitioning parameter (shear)
Em = material['Em']           # Matrix modulus [MPa]
Em_star = material['Em_star'] # Matrix degradation factor [-]
t_ply = material['t']         # Ply thickness [mm]
mat_label = material['label'] # Material name
#
from _includes.DegradedPropEvaluation import get_DegradedProp
if OmniLPF == 'YES' :
    E2,G12,nu12 = get_DegradedProp(E2,nu12,G12,v_f,Em_star,
                                    Em,eta_y,eta_s)
    nu23 = nu23 * Em_star
    G23 = E2/(2.0+2.0*nu23) # Transverse Shear Modulus [MPa]
    Omni_label = 'Omni LPF'
else:
    G23 = E2/(2.0+2.0*nu23) # Transverse Shear Modulus [MPa]
    Omni_label = 'Omni FPF'
#
nu21 = nu12*E2/E1
nu13 = nu12
nu31 = nu21
nu32 = nu23
#
XT = material['XT']     # Longitudinal tensile strength [MPa]
XC = material['XC']     # Longitudinal compressive strength [MPa]
YT = material['YT']     # Transverse tensile strength [MPa]
YBT = material['YBT']   # Biaxial transverse tensile strength [MPa]
YC = material['YC']     # Longitudinal compressive strength [MPa]
YBC = material['YBC']   # Biaxial transverse compressive strength [MPa]
SL = material['SL']     # In-plane shear strength [MPa]
ST = material['ST']     # Transverse shear strength [MPa]
beta = material['beta'] # Parameter defining the nonlinearity of the
#                         shear stress-shear strain relation [-]
#
# ---------------------- End of User INPUT --------------------------- #