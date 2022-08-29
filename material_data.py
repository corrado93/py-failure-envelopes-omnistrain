## MATERIAL DICTIONARY
##
IM78552 = {
    'label': 'IM78552',
    'E1': 171420.0,
    'E2': 9080.0,
    'v12': 0.32,
    'v23': 0.487,
    'G12': 5290.0,
    'G13': 5290.0,
    'XC': -1200.1,
    'XT': 2806.0,
    'YC': -253.7,
    'YT': 62.3,
    'SL': 89.6,
    'YBC': -600.0,
    'YBT': 38.7,
    'ST': 62.3,
    'SLP': 66.9,
    'KP': 0.05,
    'beta':4.525e-08,
    'GIC': 0.277,
    'GIIC': 0.788,
    'density': 1.58e-09,
    't': 0.125,
    'densCoh': 0.00268,
    'Knn' : 1000000.0,
    'Kss' : 1000000.0,
    'Ktt' : 1000000.0,
    'BK_exp' : 2.0,
    'Coh_t3' : 46.0,
    'Coh_t1' : 92.3,
    'Coh_GIC' : 0.35,
    'Coh_GIIC' : 1.5,
    'v_f' : 0.57,
    'eta_y' : 0.52,
    'eta_s' : 0.32,
    'Em' : 3400.0,
    'Em_star' : 0.15,
}
custom = {
    'label': 'custom_name',
    'E1': 171420.0,
    'E2': 9080.0,
    'v12': 0.32,
    'v23': 0.487,
    'G12': 5290.0,
    'G13': 5290.0,
    'XC': -1200.1,
    'XT': 2806.0,
    'YC': -437.6,
    'YT': 160.2,
    'SL': 179.2,
    'YBC': -600.0,
    'YBT': 160.0,
    'ST': 124.6,
    'SLP': 66.9,
    'KP': 0.05,
    'beta':4.525e-08,
    'GIC': 0.277,
    'GIIC': 0.788,
    'density': 1.58e-09,
    't': 0.125,
    'densCoh': 0.00268,
    'Knn' : 1000000.0,
    'Kss' : 1000000.0,
    'Ktt' : 1000000.0,
    'BK_exp' : 2.0,
    'Coh_t3' : 46.0,
    'Coh_t1' : 92.3,
    'Coh_GIC' : 0.35,
    'Coh_GIIC' : 1.5,
    'v_f' : 0.57,
    'eta_y' : 0.52,
    'eta_s' : 0.32,
    'Em' : 3400.0,
    'Em_star' : 0.15,
}
MASTERPLY = {
    'label': 'MASTERPLY',
    'E1': 160800.0,
    'E2': 9150.0,
    'v12': 0.29,
    'v23': 0.47,
    'G12': 6370.0,
    'G13': 6370.0,
    'XC': -1669.0,
    'XT': 2530.0,
    'YC': -220.0,
    'YT': 66.0,
    'SL': 93.0,
    'YBC': -600.0,      #Assumed the same as IM7/8552  
    'YBT': 38.7,        #Assumed the same as IM7/8552  
    'ST': 66.0,         #Assumed equal to YT
    'beta':2.12e-8,     #Assumed the same as IM7/8552  
    'density': 1.6e-09,
    'densCoh': 0.00268, #Assumed the same as IM7/8552  
    'Knn' : 1000000.0,
    'Kss' : 1000000.0,
    'Ktt' : 1000000.0,
    'BK_exp' : 2.0,     #Assumed the same as IM7/8552  
    'Coh_t3' : 46.0,    #Assumed the same as IM7/8552  
    'Coh_t1' : 92.3,    #Assumed the same as IM7/8552  
    'Coh_GIC' : 0.35,   #Assumed the same as IM7/8552  
    'Coh_GIIC' : 1.5,   #Assumed the same as IM7/8552  
    'v_f' : 0.55,
    'eta_y' : 0.52,    
    'eta_s' : 0.32,
    'Em' : 3400.0,
    'Em_star' : 0.15,
}
#Fibre: Silenka E-Glass 1200tex, Matrix: MY750/HY917/DY063 epoxy (Ciba Geigy)
EglassMY750 = {
    'label': 'EglassMY750',
    'E1': 45600.0,
    'E1c': 45600.0,
    'E2': 16200.0,
    'E3': 16200.0,
    'v12': 0.28,
    'v23': 0.4,
    'G12': 5830.0,
    'G13': 5830.0,
    'XC': -800.0,
    'XT': 1280.0,
    'YC': -145.0,
    'YT': 40.0,
    'SL': 73.0,
    'YBC': -145.0*(600.0/253.7),     #Scaled with the values of IM7/8552  
    'YBT': 40.0*(38.7/62.3),         #Scaled with the values of IM7/8552  
    'ST': 40.0,             #Assumed equal to YT
    'SLP': 66.9,            #Assumed the same as IM7/8552  
    'KP': 0.05,             #Assumed the same as IM7/8552  
    'beta':4.525e-08,       #Assumed the same as IM7/8552  
    'GIC': 0.165,   
    'GIIC': 0.788,          #Assumed the same as IM7/8552  
    'density': 1.58e-09,    #Assumed the same as IM7/8552  
    't': 0.250,             
    'densCoh': 0.00268,     #Assumed the same as IM7/8552  
    'Knn' : 1000000.0,      #Assumed the same as IM7/8552  
    'Kss' : 1000000.0,      #Assumed the same as IM7/8552  
    'Ktt' : 1000000.0,      #Assumed the same as IM7/8552  
    'BK_exp' : 2.0,         #Assumed the same as IM7/8552  
    'Coh_t3' : 46.0,        #Assumed the same as IM7/8552  
    'Coh_t1' : 92.3,        #Assumed the same as IM7/8552  
    'Coh_GIC' : 0.35,       #Assumed the same as IM7/8552  
    'Coh_GIIC' : 1.5,       #Assumed the same as IM7/8552  
    'v_f' : 0.6,
    'eta_y' : 0.52,         #Assumed the same as IM7/8552  
    'eta_s' : 0.32,         #Assumed the same as IM7/8552  
    'Em' : 3350.0,
    'Em_star' : 0.15,
}