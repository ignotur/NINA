fortran=g95

list = density.NE2001.o  dmdsm.NE2001.o  neclumpN.NE2001.o  neLISM.NE2001.o	nevoidN.NE2001.o  scattering98.o 
lib: $(list)
#		g++ $(list) -lf2c -lm
		ar r libNE2001.a $(list)
		ranlib libNE2001.a
main.o: main.cpp
		g++ -O3 -c main.cpp
%.o:%.f
#		g++ -c $^
		${fortran} -O3 -c $^ -o $@
%.o:%.c
		g++ -c $^
clean:
		rm $(list)
prg: libNE2001.a main.o
		g++ main.cpp -L. -lNE2001 -lf2c -lm
