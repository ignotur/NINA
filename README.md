NINA
====

Nova Investigii Neutronicorum Astrorum

A computer code to follow the evolution of large ensamble of neutron stars forming 
in isolated supernova explosions.

The code uses the electron density model YMW16 http://www.xao.ac.cn/ymw16/ http://www.atnf.csiro.au/research/pulsar/ymw16/
Yao, Manchester and Wang, 2017, Astrophys. J., 835, 29; arXiv:1610.09448

The population synthesis code is mostly based on Faucher-Giguere & Kaspi (2006) procedure with small modifications


//--------------------------------------------------------------------------------//
//                              Version agreement                                 //
//--------------------------------------------------------------------------------//

The following numeration system of version was suggested for application:

1) All version which was appeared before 12.03.2011 get versions like 0.0n, where n - 
is letter number with compilation of program:

	0.01 alpha file puls_popul_18092010.rar   - classes description
	0.02 alpha file puls_popul_02102010.rar   - elements of kinematic
	0.03 alpha file puls_popul_02102010_2.rar - life time on MS was added
	0.04 alpha file puls_popul_11102010.rar   - kinematics of galaxy was corrected
	0.05 alpha file puls_popul_14102010.rar   - life time on giant branch was added
	0.06 alpha file puls_popul_19102010.rar   - period, period derivity (standart model) 
	0.07 alpha file puls_popul_20102010.rar   - beaming
	0.08 alpha file puls_popul_15112010.rar   - new luminosity model (old)
	0.09 alpha file puls_popul_27112010.rar   - new braking model
 	0.10 alpha file puls_popul_13122010.rar   - frequenceses and their derivity
	0.11 alpha file puls_popul_12032011.rar   - new luminosity model (new)
	
Other, if it was founded, would get number like 0.0nm, where n is numbers between declared, and m order
The classes have the same numbers as packet

	0.2  alpha file puls_popul_0.2_12032011.rar  - the first in the system of version, input and output in files
        0.21 alpha file puls_popul_0.21_01052011.rar - the brightness models doesn't depend on physical properties of pulsars
	0.22 alpha file puls_popul_0.22_20072011.rar - different way to calculate the quantity of steps in the integrator, little bit quicker, random birth's time in the millenium
	0.50 alpha file puls_popul_0.50_29072011.rar - NE2001, T_sky, S_min were added i.e. this version has more complicated and more accurate radiotelescope model 
	0.70 test  file puls_popul_0.70_testing_05082011.rar - the way putting function in the files was changed, user's model parametrs classes were added, the integrator was tested, two suns, the functions out of report about testing of the pulsar population were tested, however final test haven't been done.
        0.75 alpha file puls_popul_0.75_23112011.rar - parrallel, accurate reproduction of Faucher-Giguere
	0.80 alpha file puls_popul_0.80_31102012.rar - new integrator Runge-Kutta-Fehlberg, fixed too long calculations of DM, two models of magnetic fields decay are added, trying to reproduce Parker and Swinburne surweys.
	0.83 alpha file puls_popul_0.83_12032017.rar - new electron density model ymw16

2) Version 1.0 - after final correction
3) Version 2.0 - with X-ray luminosity model

