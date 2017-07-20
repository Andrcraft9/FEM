CC=gfortran

all: compile run

clean:
	rm src/*.mod && rm src/*.o && rm program

compile:
	cd src; $(CC) -c solver.f90 -o solver.o -O3
	cd src; $(CC) -c fem.f90 -o fem.o -O3
	cd src; $(CC) -c main.f90 -o main.o -O3
	cd src; $(CC) main.o solver.o fem.o -o ../program -L /home/andr/inm_ras/ani2D-3.1/lib -lrcb2D-3.1 -lview2D-3.1 -laft2D-3.1 -lfem2D-3.1 -lmba2D-3.1 -lilu-3.1 -lblas 
run:
	./program

