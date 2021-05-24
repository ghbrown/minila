
FC=gfortran #Fortran compiler

#Source programs
PRIM_PROGS= src/prim.f90
DECOMP_PROGS=src/decomp.f90
ALL_PROGS=src/minila.f90

#Module (.mod file location)
MOD_DIR=src/mod
OBJ_DIR=src/obj
EXE_DIR=exe

#Test programs
PRIM_TEST=tests/test_prim.f90
DECOMP_TEST=tests/test_decomp.f90


all:
	@$(MAKE) -s prim
	@$(MAKE) -s decomp
	@$(FC) $(ALL_PROGS) 
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -o $(EXE_DIR)/minila.x

clean:
	@rm -f ./*.mod ./*.o ./*.x
	@rm -f $(MOD_DIR)/*.mod $(OBJ_DIR)/*.o $(EXE_DIR)/*.x


#Compile a subset of the modules
prim:
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(PRIM_PROGS) -o $(OBJ_DIR)/prim.o

decomp:
	@$(MAKE) -s prim
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(DECOMP_PROGS) -o $(OBJ_DIR)/decomp.o


#Tests
test:
	@$(MAKE) -s test-prim

test-prim:
	@$(MAKE) -s prim
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(PRIM_TEST) -o $(OBJ_DIR)/test_prim.o
	@$(FC) -I $(MOD_DIR) -o $(EXE_DIR)/test_prim.x $(OBJ_DIR)/* 
	@./$(EXE_DIR)/test_prim.x



