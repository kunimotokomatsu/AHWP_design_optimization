# AHWP design optimization
The detail of the AHWP and the functions used in this code can be found in [AHWP related functions](https://github.com/kunimotokomatsu/AHWP_related_functions).  
This code optimize the AHWP design with rondomizing the optic axis angles with a unity distribution under two conditions.  
About this conditions and optimization methods, please see [K. Komatsu et al. (2021)](http://dx.doi.org/10.1117/1.JATIS.7.3.034005) [[arXiv](https://arxiv.org/abs/2105.05561)].
**Please edit line 24-44 to suit your case.**   

## Parameters

num: number of randomizing (integer)  
patterns: conditions for the optimization called anti-symmetric (frequency independent phase) and symmetric (conventional): ["asym","sym"]  
layers: integer array of number of layers you use for the optimization: [3,5,7,9]  
no: efractive index for ordinary ray of the birefringent material used for AHWP (sapphire at 1.5K: 3.047)  
ne: efractive index for extraordinary ray of the birefringent material used for AHWP (sapphire at 1.5K: 3.361)  
a_in: angle between a incident polarization and a detector sensitive angle [rad] (this value not affect to optimized angle set and only change a offset in a phase plot)  
p_in: ratio of polarization in a incidnet radiation  
f_range: array of initial and end frequency sets [Hz] you use in the optimization: [[4.e9,191.e9],[34.e9,161.e9,],[84.e9,111.e9]]
f_step: freqency step [Hz] in each frequency range: 1.e9
thickness: thickness of each plate (in this code, we use same thickness for all plates and fix them): 4.9e-3  
output_dir: output directory name (if there is no directory, it is automatically made): 'output' 

## Contents included in npz file  

This code output figures and npz file for each case you set.  
They are named as '''opt_'+pattern+'_'+str(layer)+'_run'+str(num)+'_'+str(int(f_i*1.e-9))+'_'+str(int(f_f*1.e-9))+'_odd'''.  
layer, num, no, ne, a_in, p_in, freq_opt, thickness: parameters you input  
poleff, phi4: polarization efficiency and phase calculated using a optimized angle set  
angle_best: a optimized angle set  
ave_best: band averaged polarization efficiency calculated using a optimized angle set  
ave_arr: band averaged polarization efficiency calculated using each optimized angle set  
random_seed: seed of randomizing  
random_angle_arr: angle set used for each randomizing  