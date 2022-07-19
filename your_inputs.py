# ------------------------------ User INPUT ----------------------------------- #

OmniLPF = 'NO'          # 'YES' or 'NO': 
#                         Omni last-ply failure (LPF) envelopes are generated 
#                         using degraded properties. If you set 'NO', omni  
#                         first-ply failure (FPF) envelopes will be generated.
ExcelOutput = 'YES'     # 'YES' or 'NO'
PlotFigure = 'NO'       # 'YES' or 'NO'
SAVEFIG = 'NO'          # 'YES' or 'NO' 
#                         Note: This is a suboption of PlotFigure! If you want to save the figure,
#                               you need to set PlotFigure = 'YES'
from numpy import *
import _includes.material_data
#
## Material properties
#
# Choose the material from the database (_include/material_data.py) or 
# set your custom material.
#          Materials = IM78552, custom ...
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