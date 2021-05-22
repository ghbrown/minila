
FC=gfortran #Fortran compiler

PRIMITIVE_MODULES=primitives.f90
DECOMP_MODULES=decomp.f90


all:
	@$(FC) primitives.f90 decomp.f90 minila.f90 -o minila.x

primitives:
	@$(FC) -c $(PRIMITIVE_MODULES)

decomp:
	@$(FC) -c $(PRIMITIVE_MODULES) $(DECOMP_MODULES)
test:
	@$(FC) $(PRIMITIVE_MODULES) test_primitives.f90 -o testprim; ./testprim; rm testprim
clean:
	@rm ./*.o ./*.mod

