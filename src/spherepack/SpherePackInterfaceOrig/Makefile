# on 64 bit machines (like alphas, itaniums, crays) use FPPFLAGS=-DINT64
# flags for Sun Forte compiler
FFLAGS= -e -fast
#FFLAGS= -fast -xarch=v9a -xtypemap=real:64
LINKS= -L/home/jsw/lib -lspherepack
#LINKSDP= -L/data/jsw/lib -lspherepack_dp
FC = f95
PROMOTEDP = -xtypemap=real:64
# flags for alpha linux
#FFLAGS= -fast -O4
#LINKS= -L/home/jsw/lib -lspherepack
#LINKSDP= -L/home/jsw/lib -lspherepack_dp
#FC = fort
#PROMOTEDP = -r8
# flags for intel linux with intel compiler.
#FFLAGS= -O3 -tpp7 -xKW -pad -ip -Vaxlib
#LINKS= -L/home/whitaker/lib -lspherepack
#LINKSDP= -L/home/whitaker/lib -lspherepack_dp
#FC = ifc
#PROMOTEDP = -r8

all: 
	$(FC) $(FFLAGS) -c sphere.f90

dp: 
	$(FC) $(FFLAGS) $(PROMOTEDP) -c sphere.f90

clean:
	-rm sphere.o test_sphere

test:
	$(FC) $(FFLAGS) -o test_sphere test_sphere.f sphere.o $(LINKS)
	./test_sphere      

testdp:
	$(FC) $(FFLAGS) $(PROMOTEDP) -o test_sphere test_sphere.f sphere.o $(LINKSDP)
	./test_sphere      
