
FC=gfortran #Fortran compiler

all:
	@$(FC) primitives.f90 decomp.f90 minila.f90 -o minila.x

primitives:
	@$(FC) -c primitives.f90

decomp:
	@$(FC) -c decomp.f90
