# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

#COMPILER=ifort -free -xHost -shared-intel -align all  -debug all -qopt-report-phase=vec -diag-disable 8291 -diag-disable 8290  
COMPILER=gfortran -ffree-line-length-none -m64 -Wall -fPIC
OPTS=-O1 -mcmodel=medium -lm -llapack -lblas -lfftw3
OPTSNETCDF=-lnetcdff -I/usr/include -L/usr/lib 

SPECFAB=specfab
MOMENTS=moments
GAUNT=gaunt
TENPROD=tensorproducts

DEMO=demo
DEMOSRCS=$(wildcard $(DEMO)/demo*.f90)
DEMOS=$(DEMOSRCS:.f90=)

ELMERSRCS=$(wildcard $(DEMO)/elmer*.f90)
ELMERS=$(ELMERSRCS:.f90=)

ALLOBJS=$(SPECFAB).o $(TENPROD).o $(MOMENTS).o $(GAUNT).o
ALLMODS=$(ALLOBJS:.o=.mod)
ALLSRCS=$(ALLOBJS:.o=.f90)

.PHONY: all clean clear demo demoso elmer lib
########################
all: lib$(SPECFAB).so $(ALLOBJS) $(ALLMODS) $(DEMO)

lib: lib$(SPECFAB).so

elmer: $(ELMERS)

$(DEMO)/elmer%: $(DEMO)/elmer%.f90 | lib$(SPECFAB).so $(ALLMODS)
	$(COMPILER) $< $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $@

demo: $(DEMOS) $(SPECFAB).o
	mkdir -p demo/solutions
	@echo "\n*** To get going, try the demos:"
	@echo "\n--- Lattice rotation ---"
	@echo "cd demo; ./demo_ROT uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd demo; ./demo_ROT ss_xz :: simple shear (ss) along the x--z plane"
	@echo "\n--- Dynamic recrystallization ---"
	@echo "cd demo; ./demo_DRX uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd demo; ./demo_DRX ss_xz :: simple shear (ss) along the x--z plane"
	@echo "\n--- Both Dynamic recrystallization and lattice rotation ---"
	@echo "cd demo; ./demo_FULL uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd demo; ./demo_FULL ss_xz :: simple shear (ss) along the x--z plane"

$(DEMO)/demo%: $(DEMO)/demo%.f90 | $(ALLOBJS) $(ALLMODS)
	$(COMPILER) $< $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $@

# DL: I do not see a need for this one?
demoso: lib$(SPECFAB).so
	$(COMPILER) demo/demo_ROT.f90 -L./ -lspecfab $(OPTS) $(OPTSNETCDF) -o demo/demo_ROT
	mkdir -p demo/solutions

$(SPECFAB)py: $(SPECFAB).o
	f2py -lm -llapack -lblas -lfftw3 -I. $(ALLOBJS) -c -m specfabpy specfabpy.f90 --f90flags="-ffree-line-length-none -mcmodel=medium" --quiet 
	mkdir -p demo/solutions
	cp specfabpy.cpython* demo/
	@echo "\n*** To get going, try the demo:"
	@echo "cd demo; python3 demo_ROT.py uc_zz ::: for uniaxial compression (uc) in the vertical (z)"
	@echo "cd demo; python3 demo_ROT.py ss_xz ::: for simple shear (ss) along the x--z plane"

########################

$(SPECFAB).o: $(MOMENTS).o $(GAUNT).o $(SPECFAB).f90
	$(COMPILER) $(OPTS) -c $(TENPROD).f90 
	$(COMPILER) $(OPTS) -c $(SPECFAB).f90

lib$(SPECFAB).so: $(ALLSRCS) $(ALLOBJS)
	$(COMPILER) $(OPTS) -shared $(ALLSRCS) -o $@


########################

$(MOMENTS).o: $(MOMENTS).f90
	@echo "***************************************************************************************************"
	@echo "*** Compiling structure tensor expressions... this may take some time but is required only once ***"
	@echo "***************************************************************************************************"
	$(COMPILER) $(OPTS) -c $(MOMENTS).f90

$(GAUNT).o: $(GAUNT).f90
	@echo "*****************************************************************************************"
	@echo "*** Compiling gaunt coefficients... this may take some time but is required only once ***"
	@echo "*****************************************************************************************"
	$(COMPILER) $(OPTS) -c $(GAUNT).f90

clean:
	rm -f demo/demo_ROT demo/demo_DRX demo/demo_FULL *.o *.mod *.so
	
clear:
	rm -f demo/demo_ROT demo/demo_DRX demo/demo_FULL $(SPECFAB).o $(SPECFAB).mod *.so

