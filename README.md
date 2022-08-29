## Description of the script

This script was developed to generate *extended omni-strain failure envelopes* (G. Corrado et al., An extended invariant approach to laminate failure of fibre-reinforced polymer structures. The Aeronautical Journal, https://doi.org/10.1017/aer.2021.121, available in U.Porto's repository: https://hdl.handle.net/10216/139464).
In this script, the failure loci are computed using the 3D invariant-based failure theory for all possible ply orientations, where the number of ply orientations (from 0º to 90º) can be modified by the user. The omni-strain failure envelope is the minimum envelope among the ones generated.
Results can be visualized and saved in different ways: using pyplot or excel. User inputs allow to opt for the generation of these outputs.

## Requirements

1. Python 3
2. Office package installed (to save results in an excel file)

## Instructions to use this file

1. Check the requirements above.
2. Modify the script **your_inputs.py**, based on the your curiosity/needs. Additional info can be found in the script and in the papers in reference.
3. Run the script **filetorun.py**. This can be done on the Terminal, for instance using PowerShell if you have Windows.

## Detailed description of the workflow

Herein, for those interested in more details on the script, there is a sequence of steps to describe *'how the script is able to do what it does'*.

1. Firstly, it loads all the user inputs, including the material properties. If omni-strain LPF envelopes are selected, it calculates the degraded material properties. Otherwise, we skip this point and the intact material properties are use instead.
2. Next, it computes the failure parameters used in the 3D invariant-based failure theory, namely: a1, a2, a32T, a32C, a3T, a3C.
3. Then, it calculates the stiffness matrix of the lamina of the selected material. Since the 3D invariant-based failure theory is formulated in stress space, failure envelopes are firstly generated in stress space and sequently expressed in strain space to find omni-strain FPF/LPF envelopes.
4. Next, the laminate stiffness matrix [A] is computed by using the laminate layup. Please, note that omni-strain FPF/LPF envelopes are invariants with respect the laminate stacking sequence composition, since they are expressed in strain space. However, to transform the obtained omni-strain envelopes in stress space, the laminate stiffness is required. 
5. Then, the script computes the failure loci in strain space based on the 3D invariant-based failure theory, using the bisection method as root-finding method.
In particular, in order to consider all the possible combinations of applied loading in stress or strain space, the script is set to search each failure point using a polar coordinate system and, thus, to apply the bisection method for each angular position. The number of failure loci will coincide with the number of angular increments, by which the stress or strain space is divided. The bisection method is applied through the function **get_delta_F()**
6. Finally, the results are locally saved in lists and, based on the user inputs, figures can be saved in a .png file and/or the obtained failure loci can be saved and plotted in a .xlsx file. 
