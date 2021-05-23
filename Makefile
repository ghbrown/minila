
FC=gfortran #Fortran compiler

#Source programs
PRIM_PROGS= src/prim.f90
DECOMP_PROGS=src/decomp.f90
ALL_PROGS=src/minila.f90

#Module (.mod file location)
MOD_DIR=src/mod
OBJ_DIR=src/obj
EXE_DIR=src/exe

#Test programs
PRIM_TEST=tests/test_prim.f90
DECOMP_TEST=tests/test_decomp.f90


all:
	@make prim
	@make decomp
	@$(FC) $(ALL_PROGS) 
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -o $(EXE_DIR)/minila.x

clean:
	@rm $(MOD_DIR)/* $(OBJ_DIR)/*


#Compile a subset of the modules
prim:
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(PRIM_PROGS) -o $(OBJ_DIR)/prim.o

decomp:
	@make prim
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(DECOMP_PROGS) -o $(OBJ_DIR)/decomp.o


#Tests
test:
	@make test-prim

test-prim:
	@make prim
	@$(FC) -I $(MOD_DIR) $(PRIM_TEST) -o test_prim.x
	@./test_prim.x
	@rm testprim.x



