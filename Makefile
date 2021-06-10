
FC=gfortran #Fortran compiler

#Source programs
PRIM_PROGS= src/prim.f90
DECOMP_PROGS=src/decomp.f90
LINSYS_PROGS=src/linsys.f90 src/linsys_linsol.f90 src/linsys_inv.f90 src/linsys_sor.f90
ALL_PROGS=src/minila.f90

#Module (.mod file location)
MOD_DIR=src/mod
OBJ_DIR=src/obj
EXE_DIR=exe

#Test programs and folders
PRIM_TEST=tests/test_prim.f90
DECOMP_TEST=tests/test_decomp.f90
LINSYS_TEST=tests/test_linsys.f90


all:
	@$(MAKE) -s prim
	@$(MAKE) -s decomp
	@$(MAKE) -s linsys
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(ALL_PROGS) -o $(OBJ_DIR)/minila.o


clean:
	@rm -f ./*.mod ./*.o ./*.x
	@rm -f $(MOD_DIR)/*.mod $(OBJ_DIR)/*.o $(EXE_DIR)/*.x


#Compile a subset of the modules
prim:
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(PRIM_PROGS) -o $(OBJ_DIR)/prim.o

decomp:
	@$(MAKE) -s prim
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(DECOMP_PROGS) -o $(OBJ_DIR)/decomp.o

linsys:
	@$(MAKE) -s decomp
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(LINSYS_PROGS)


#Tests
test:
	@$(MAKE) -s test-prim
	@$(MAKE) -s test-decomp
	@$(MAKE) -s test-linsys

#all of these tests-* blocks are essentially the same up to a single string (prim, decomp, linsys) so perhaps they could be intelligently condensed?
test-prim:
	@$(MAKE) -s prim
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(PRIM_TEST) -o $(OBJ_DIR)/test_prim.o
	@$(FC) -I $(MOD_DIR) -o $(EXE_DIR)/test_prim.x $(OBJ_DIR)/* 
	@./$(EXE_DIR)/test_prim.x
	@$(MAKE) -s clean

test-decomp:
	@$(MAKE) -s decomp
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(DECOMP_TEST) -o $(OBJ_DIR)/test_decomp.o
	@$(FC) -I $(MOD_DIR) -o $(EXE_DIR)/test_decomp.x $(OBJ_DIR)/*
	@./$(EXE_DIR)/test_decomp.x
	@$(MAKE) -s clean

test-linsys:
	@$(MAKE) -s linsys
	@$(FC) -I $(MOD_DIR) -J $(MOD_DIR) -c $(LINSYS_TEST) -o $(OBJ_DIR)/test_linsys.o
	@$(FC) -I $(MOD_DIR) -o $(EXE_DIR)/test_linsys.x $(OBJ_DIR)/*
	@./$(EXE_DIR)/test_linsys.x
	@$(MAKE) -s clean

