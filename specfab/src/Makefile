# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020-2021

#COMPILER=ifort -free -xHost -shared-intel -align all  -debug all -qopt-report-phase=vec -diag-disable 8291 -diag-disable 8290  
COMPILER=gfortran -ffree-line-length-none -m64 -Wall -fPIC
OPTS=-O1 -mcmodel=medium -lm -llapack -lblas -lfftw3
OPTSNETCDF=-lnetcdff -I/usr/include -L/usr/lib 

SPECFAB=specfab
MOMENTS=moments
GAUNT=gaunt
TENPROD=tensorproducts

ALLOBJS=$(SPECFAB).o $(TENPROD).o $(MOMENTS).o $(GAUNT).o
ALLMODS=$(ALLOBJS:.o=.mod)
ALLSRCS=$(ALLOBJS:.o=.f90)

.PHONY: all clean clear demo lib
########################
# Phony targets
########################
all: lib$(SPECFAB).so $(ALLOBJS) $(ALLMODS)

lib: lib$(SPECFAB).so

clean:
	rm -f *.o *.mod *.so
	
clear:
	rm -f $(SPECFAB).o $(SPECFAB).mod *.so
########################
# Top-level targets
########################

$(SPECFAB).o: $(MOMENTS).o $(GAUNT).o $(SPECFAB).f90
	$(COMPILER) $(OPTS) -c $(TENPROD).f90 
	$(COMPILER) $(OPTS) -c $(SPECFAB).f90

lib$(SPECFAB).so: $(ALLSRCS) $(ALLOBJS)
	$(COMPILER) $(OPTS) -shared $(ALLSRCS) -o $@


########################
# Dependencies
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
