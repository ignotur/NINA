comp=gcc
files=dora.c frb_d.c gum.c lmc.c nps.c spiral.c thin.c ymw16par.c fermibubble.c galcen.c lib_ymw16.c localbubble.c smc.c thick.c 
flags= -dynamiclib -flat_namespace -O3
libname=-o libymw16.so

main: 
		${comp} ${flags} ${files} ${libname}
		g++  test.cpp -L./ -lymw16
#		${comp} ${flags} lib_simpl.c ${libname}
clean:
		rm *.so 
