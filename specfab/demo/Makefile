# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020-2021

#COMPILER=ifort -free -xHost -shared-intel -align all  -debug all -qopt-report-phase=vec -diag-disable 8291 -diag-disable 8290  
COMPILER=gfortran -ffree-line-length-none -m64 -Wall -fPIC
OPTS=-O1 -mcmodel=medium -lm -llapack -lblas -lfftw3 -I../src -L../src
OPTSNETCDF=-lnetcdff -I/usr/include -L/usr/lib

SPECFAB=specfab
MOMENTS=moments
GAUNT=gaunt
TENPROD=tensorproducts

DEMOSRCS=$(wildcard demo*.f90)
DEMOS=$(DEMOSRCS:.f90=)

SRCDIR=../src/
ALLOBJS=$(SRCDIR)$(SPECFAB).o $(SRCDIR)$(TENPROD).o $(SRCDIR)$(MOMENTS).o $(SRCDIR)$(GAUNT).o
ALLMODS=$(ALLOBJS:.o=.mod)


.PHONY: all clean clear demo demoso elmer lib
########################
# Phony targets
########################
lib:
	@echo "No lib in demo"

all: $(DEMO)

clean:
	rm -f $(DEMOS)

clear:
	rm -f $(DEMOS)

demo: $(DEMOS) $(ALLOBJS)
	mkdir -p solutions
	@echo "\n*** To get going, try the demos:"
	@echo "\n--- Lattice rotation ---"
	@echo "cd demo; ./demo_LATROT uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd demo; ./demo_LATROT ss_xz :: simple shear (ss) along the x--z plane"
	@echo "\n--- Dynamic recrystallization ---"
	@echo "cd demo; ./demo_DRX uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd demo; ./demo_DRX ss_xz :: simple shear (ss) along the x--z plane"
	@echo "\n--- Both Dynamic recrystallization and lattice rotation ---"
	@echo "cd demo; ./demo_FULL uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd demo; ./demo_FULL ss_xz :: simple shear (ss) along the x--z plane"
########################
# Real targets
########################

demo%: demo%.f90 | $(ALLOBJS) $(ALLMODS)
	$(COMPILER) $< $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $@
