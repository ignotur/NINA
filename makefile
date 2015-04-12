#####################################################################
##            Makefile for population synthesis code               ##
#####################################################################
# The name of model contains a letter and roman digit 
# e.g.  make AIII
# Буквы - модель светимости
# Цифры - модель убывание поля
# A   - модель со сгустками плазмы
# B   - модель с нормальным распределением энергии в конусе излучения
# C   - модель с плоским распределением излучения в конусе радио пульсара, как было у Faucher
# I   - модель без убывания магнитного поля на поверхности
# II  - модель со ступенчатым убыванием поля
# III - модель с убыванием поля по статье Geppert
# IV  - модель с первоначальным экспоненциальным спадом и дальнейшим неубыванием поля
# V   - модель двухкомпонентного убывания напряжённости магнитного поля за счёт эффекта Холла и омического распада, с учётом изменения угла меджу осями в соотвествии с вакуумной магнитосферой

list = density.NE2001.o  dmdsm.NE2001.o  neclumpN.NE2001.o  neLISM.NE2001.o     nevoidN.NE2001.o  scattering98.o
#list_mono  = class.cpp new_main.cpp math_func.cpp astro_func.cpp Runge_Kutta.cpp class_sun.cpp class_neutron.cpp class_T_map.cpp
list_mono  = class.cpp comm_main.cpp math_func.cpp astro_func.cpp Runge_Kutta_Fehlberg.cpp class_sun.cpp class_neutron.cpp class_T_map.cpp bp.cpp messages.cpp syntax.cpp without_field_decay.cpp  flat_distr_lum_model_mod.cpp ns_methods.cpp step_field_decay.cpp Geppert_field_decay.cpp expon_field_decay.cpp norm_distr_lum_model.cpp  old_pons_field_decay.cpp 
list_para  = class.cpp omp_main.cpp math_func.cpp astro_func.cpp Runge_Kutta_Fehlberg.cpp class_sun.cpp class_neutron.cpp class_T_map.cpp bp.cpp
list_otl   = class.cpp otl.cpp math_func.cpp astro_func.cpp Runge_Kutta.cpp class_sun.cpp class_neutron.cpp class_T_map.cpp bp.cpp
list_param = param_main.cpp class_T_map.cpp class.cpp math_func.cpp astro_func.cpp Runge_Kutta.cpp class_sun.cpp class_neutron.cpp bp.cpp
list_A     = sparks_fun.cpp
list_B     = norm_distr_lum_model.cpp 
list_C     = 
#list_C     = flat_distr_lum_model.cpp
list_I     =  
list_II    = step_field_decay.cpp
list_III   = Geppert_field_decay.cpp 
list_IV    = expon_field_decay.cpp
#list_V     = two_comp_field_decay.cpp
list_V     = two_comp_field_decay_mod.cpp
keys = 
main:         ./ne2001_f/libNE2001.a 	$(list_mono) $(list_C) $(list_I)
		g++  -O3 $(list_mono) -o population.out  -L./ne2001_f/ -lNE2001 -L/opt/local/lib/ -lf95 -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas -lm $(keys)
./ne2001_f/libNE2001.a: 
		cd ne2001_f/; make
