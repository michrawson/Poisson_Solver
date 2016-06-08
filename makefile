
FFLAGS=-freal-4-real-8 -fPIC -O3
#FFLAGS=-freal-4-real-8 -fPIC -fbounds-check -fbacktrace -Wall -Wextra  -pedantic  -fcheck=all -Warray-temporaries -fimplicit-none -ffree-line-length-0 -ffpe-trap=zero,overflow -finit-real=nan -Wconversion -pg -g -O0    

SRC08=$(wildcard *.f08)

SRC03=$(wildcard *.f03)

SRC95=$(wildcard *.f95)

SRC90=$(wildcard *.f90)

SRC77=$(wildcard *.f)

.PHONY: fmodule.so main.tsk

all : main.tsk fmodule.so

main.tsk :
	gfortran $(FFLAGS) $(SRC77) $(SRC90) $(SRC95) $(SRC03) $(SRC08) -o main.tsk

fmodule.so :
	f2py --fcompiler=gfortran -c $(SRC77) $(SRC90) -m fmodule
# --debug-capi

%.o : %.f95
	gfortran -c $(FFLAGS) $< -o $@

%.o : %.f
	gfortran -c $(FFLAGS) $< -o $@

clean:
	rm -f *.mod *.o *.so *.tsk *.tsk.dSYM *.so.dSYM
	
	