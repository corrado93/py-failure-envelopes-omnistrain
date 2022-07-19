## Description of the script

This script was developed to generate extended omni-strain failure envelopes (G. Corrado et al., An extended invariant approach to laminate failure of fibre-reinforced polymer structures. The Aeronautical Journal, https://doi.org/10.1017/aer.2021.121, available in U.Porto's repository: https://hdl.handle.net/10216/139464).
In this script, the failure loci are computed using the 3D invariant-based failure theory for all possible ply orientations, where the number of ply orientations (from 0º to 90º) can be modified by the user.
Results can be visualized in 3 ways: using pyplot, using matplotlib and using excel.

## Requirements

1. Python, version ...
2. Python libraries installed (?)

## Instructions to use this file

1. Check which version of Python you have (*you need a version more recent than py ---*). If you don't have python, you need to install it and verify that it works properly before going to the next step.
2. Modify the script **omni_env_script.py**, based on your needs. Additional info can be found in the script.
3. Run the script. This can be done on the Terminal, for instance using PowerShell if you have Windows.

## Detailed description of the script

Herein, for those interested, there is a sequence of steps to describe *'how the script is able to do what it does'*.

1. Firstly, we compute the failure parameters used in the 3D invariant-based failure theory, namely: a1, a2, a32T, a32C, a3T, a3C.  
2. Then, we compute the stiffness matrix for each orientation of the studied laminate.
3. Then, we compute the following invariants: U1, U2, U4, tr[Q]
4. 

