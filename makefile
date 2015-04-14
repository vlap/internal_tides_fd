# Path for source files
VPATH = src/baro_utils:src/defs:src/generate:src/misc:src/shtools:src/sparse:src/spherepack
# Defaults
compiler  = gnu#		intel or gnu
debug     = work#		work or optim
architecture= intel64#	intel64 or ia32
interface = lp64#		lp64 or ilp64 (long integers)
threading = iomp# 	iomp, gomp or sequential
sd        = dynamic#		static or dynamic


# -----------------------------------------------------------------------------
# Significant suffixes [assuming Fortran 90 (.f90) and possibly Fortran 77
# (.f) source code]:

.SUFFIXES : .o .mod .f90 .f

# -----------------------------------------------------------------------------
# Check on which machine the file is being compiled

hostname | { read machine; echo machine=$machine; }

# -----------------------------------------------------------------------------
# (a) Set Fortran version and compiler options:               (COMPILER MACROS)
ifeq ($(compiler),gnu) #******************* GFORTRAN ************************
  FC = gfortran
  OPTS = -fimplicit-none -ffree-form -ffixed-line-length-none -mcmodel=medium #-frecord-marker=8-fcray-pointer  -fno-second-underscore #-x f77-cpp-input
  #FPPSTOP=-x none
  
#Include directory with mod files
  MOD_INC = -Jmod/
  OPTS += $(MOD_INC)
  
ifeq ($(debug),work)
  OPTS += -g -w -fbounds-check -fbacktrace -finit-real=nan # -Wall (includes -Wunused -Wuninitialized) -Wextra
else
  OPTS += -w -DNO_DEBUG -malign-double -ffast-math -fno-cx-limited-range -O3 -funroll-loops --param max-unroll-times=4 # ! -w  :: Inhibit all warning messages.
endif 

ifeq ($(architecture),ia32) # This option tells compiler to generate code for IA-32 architecture.
  ifeq ($(interface),ilp64)
    $(warning  *** ILP64 interface is not available for MKL-IA32)
    $(error Try >make help)
  else
    OPTS += -m32  
  endif
  IFACE_SUFF=
else # This option tells compiler to generate code for Intel64 architecture.
  ifeq ($(interface),ilp64)
    OPTS += -fdefault-integer-8 -m64 
  else
    OPTS += -m64 
  endif
  IFACE_SUFF=_$(interface)
endif

ifeq ($(threading),gomp)
  OPTS += -fopenmp
  IFACE_THREADING_PART = gnu
else
ifeq ($(threading),iomp)
  IFACE_THREADING_PART = intel
endif
endif 
  
  IFACE_COMP_PART = gf
  
else #******************* IFORT ************************
  FC = ifort
#Include directory with mod files
  MOD_INC = -module mod/
  
  OPTS = -w -132 -shared-intel -mcmodel=medium -i-dynamic#-fpp
  OPTS += -implicitnone -free $(MOD_INC) -sox

ifeq ($(debug),work)
  OPTS += -g -check all -fpe0 -warn -mp -traceback -debug extended
else
  OPTS += -O2 -fast -ipo  -IPF-fp-relaxed
endif

ifeq ($(architecture),ia32)
  ifeq ($(interface),ilp64)
    $(warning  *** ILP64 interface is not available for MKL-IA32)
    $(error Try >make help)
  else
    # This option tells the compiler to generate optimized code for Pentium or later processor.
    # If you don't need it, you can remove this option.
    OPTS += -mia32
  endif
  IFACE_SUFF=
else
  IFACE_SUFF=_$(interface)
endif

ifeq ($(threading),gomp)
    $(warning  *** ILP64 interface is not available for MKL-IA32)
    $(error Try >make help)
else
ifeq ($(threading),iomp)
  OPTS += -openmp -parallel -par-report2
  IFACE_THREADING_PART = intel
endif
endif
  
  IFACE_COMP_PART = intel
endif

#******** Select dynamic or static linking: ************
ifeq ($(sd),static)
  EXT=a
  RES_EXT=lib
  S_GRP=-Wl,--start-group
  E_GRP=-Wl,--end-group
else
  EXT=so
  RES_EXT=so
  S_GRP=
  E_GRP=
endif

# -----------------------------------------------------------------------------
# (b) Set library paths & names, and include file paths:  (PATH AND LIB MACROS)
#	IMPORTANT:
#	These variables are already set up as environmental variables in bash_rc

#MY_LIB_PATH = /home/amt6/amtvl/Desktop/scratch/lib
#MY_INC_PATH = /home/amt6/amtvl/Desktop/scratch/lib



#Place of MKL libraries
ifndef MKLROOT
#---------------------------
ifeq ($(HOSTNAME),theon)
MKLROOT=/opt/intel/composerxe-2011.5.220/mkl
else
#***************************
ifeq ($(HOSTNAME),xigor)
MKLROOT=/scratch/a/users/amtvl/lib/intel/composerxe-2011.5.220/mkl

#MKLROOT=/opt/intel/Compiler/11.1/064/mkl
# somehow mkl directory is a bit wield on xigor
#ifeq ($(architecture),ia32)
#architecture=32
#else
#ifeq ($(architecture),intel64)
#architecture=em64t
#endif
#endif
endif
#***************************
endif
#---------------------------
endif
MKL_INC = -I$(MKLROOT)/include
#MKL_INC_EXT = -I$(MKLROOT)/include/$(architecture)/$(interface)
MKL_LIB = -L$(MKLROOT)/lib/$(architecture)
MKL_PATH = $(MKLROOT)/lib/$(architecture)
CMPLR_PATH = $(MKLROOT)/../compiler/lib/$(architecture)

IFACE_LIB=$(MKL_PATH)/libmkl_$(IFACE_COMP_PART)$(IFACE_SUFF).$(EXT) # IFACE_LIB=$(MKL_PATH)/libmkl_intel_lp64.a
CORE_LIB=$(MKL_PATH)/libmkl_core.$(EXT) # CORE_LIB=$(MKL_PATH)/libmkl_core.a
ifeq ($(threading),sequential)
  THREADING_LIB=$(MKL_PATH)/libmkl_sequential.$(EXT)
  OMP_LIB =
else
  THREADING_LIB=$(MKL_PATH)/libmkl_$(IFACE_THREADING_PART)_thread.$(EXT) # THREADING_LIB=$(MKL_PATH)/libmkl_intel_thread.a
  OMP_LIB = -L$(CMPLR_PATH) -liomp5
  #OMP_LIB = $(CMPLR_PATH)/libiomp5.$(EXT)
#***************************
#ifeq ($(HOSTNAME),xigor)
#ifeq ($(architecture),32)
#  OMP_LIB = -L$(MKLROOT)/../lib/ia32/ -liomp5
#else
#ifeq ($(architecture),em64t)
#OMP_LIB = -L$(MKLROOT)/../lib/intel64/ -liomp5
#endif
#endif
#endif
#***************************
endif

#Place of Sparse BLAS objects and modules
#SB_LIB = -L${MY_LIB_PATH}/spblas/SOFTWARE -lSparseBLAS
#SB_INCL = -I${MY_INC_PATH}/spblas/SOFTWARE

#Place of numeric libraries
MATH_LIBS = -llapack -lblas -lgfortran -lm
 
#Place of UMFPACK libraries
L_UMF = -L${MY_LIB_PATH}/SuiteSparse/UMFPACK/Lib -lumfpack
L_AMD = -L${MY_LIB_PATH}/SuiteSparse/AMD/Lib -lamd
L_CHOL = -L${MY_LIB_PATH}/SuiteSparse/CHOLMOD/Lib -lcholmod
L_COL = -L${MY_LIB_PATH}/SuiteSparse/COLAMD/Lib -lcolamd
UMFPACK_LIB = ${MY_LIB_PATH}/SuiteSparse/UMFPACK/Demo/umf4_f77zwrapper.o \
          -lblas -lgfortran -lm $(L_UMF) $(L_AMD) $(L_CHOL) $(L_COL)

I_UMF_F = -I${MY_INC_PATH}/SuiteSparse/UMFPACK/Lib -lumfpack

#Place of NetCDF libraries and modules
ifeq ($(compiler),gnu)
NETCDF_LIB = -L${MY_LIB_PATH}/netcdf/gfortran_v/lib/ -lnetcdf
NETCDF_INC = -I${MY_INC_PATH}/netcdf/gfortran_v/include/
else
NETCDF_LIB = -L${MY_LIB_PATH}/netcdf/ifort_v/lib/ -lnetcdf
NETCDF_INC = -I${MY_INC_PATH}/netcdf/ifort_v/include/
endif

#Place of SPARSKIT library
SKIT_LIB = -L${MY_LIB_PATH}/SPARSKIT2 -lskit
#Place of Ngmath library (part of the NCAR Graphics (NCL))
#NGMATH_LIB = -L${MY_LIB_PATH}/NCL/lib -lngmath

#L_NLC = ${MY_LIB_PATH}/NCL/lib -lncarg -lncarg_gks -lncarg_c
#L_X11 = -L/usr/lib64/ -lX11

#------------------------------------------------------------------------------
# (c) Set program (executable file) name:                      (PROGRAM TARGET)

        PROGRAM = igw_v1

# -----------------------------------------------------------------------------
# (d) Create expected object file list, with .o and .mod extension: [TARGETS OF PRE-
#     REQUISITE (.f90,.f) SOURCE FILES]

SOURCES=igw_v1.f90 
OBJECTS=$(SOURCES:%.f90=obj/%.o)

# high level modules
MODULES_H_S=projection.f90 baro_domain_integrals.f90 itm_solver.f90 baro_solver.f90
MODULES_H_O=$(MODULES_H_S:%.f90=obj/%.o)
MODULES_H_M=$(MODULES_H_S:%.f90=mod/%.mod)
# medium level modules
MODULES_M_S=read_etopo.f90 generate_global_grid.f90 itd.f90 generate_matrices.f90 sal.f90 unsym_solvers.f90
MODULES_M_O=$(MODULES_M_S:%.f90=obj/%.o)
MODULES_M_M=$(MODULES_M_S:%.f90=mod/%.mod)
# low level modules
MODULES_L_S=save_load.f90 interpolate.f90 control_file.f90
MODULES_L_O=$(MODULES_L_S:%.f90=obj/%.o)
MODULES_L_M=$(MODULES_L_S:%.f90=mod/%.mod)
# my sparse subroutines
MYSPARSE_L_S= my_blas.f90 my_lapack.f90 my_sparse.f90 my_sparse_aggregate.f90
MYSPARSE_L_O=$(MYSPARSE_L_S:%.f90=obj/%.o)
MYSPARSE_L_M=$(MYSPARSE_L_S:%.f90=mod/%.mod)
# my sparse subroutines
MISC_L_S= my_trigs.f90 err_manager.f90
MISC_L_O=$(MISC_L_S:%.f90=obj/%.o)
MISC_L_M=$(MISC_L_S:%.f90=mod/%.mod)
# parameter level modules
ifeq ($(HOSTNAME),theon)
MODULES_P_S=precisions.f90 iso_varying_string.f90 my_paths_theon.f90
else
ifeq ($(HOSTNAME),xigor)
MODULES_P_S=precisions.f90 iso_varying_string.f90 my_paths_xigor.f90
endif
endif 
MODULES_P_O=$(MODULES_P_S:%.f90=obj/%.o)
MODULES_P_M=$(MODULES_P_S:%.f90=mod/%.mod)
#-----------------------------------------------------------------------------#
# tmp mkl
FGMRES_MKL_S=iter_solvers.f90
FGMRES_MKL_O=$(FGMRES_MKL_S:%.f90=obj/%.o)
FGMRES_MKL_M=$(FGMRES_MKL_S:%.f90=mod/%.mod)
# sparskit subroutines
SPARSKIT_S=blassm.f formats.f unary.f ilut.f
SPARSKIT_O=$(SPARSKIT_S:%.f=obj/%.o)
# spherepack subroutines
SPHEREPACK_IFACE_S=SHRotate.f90 spherepack_iface.f90 testshrotate.f90
SPHEREPACK_IFACE_O=$(SPHEREPACK_IFACE_S:%.f90=obj/%.o)
SPHEREPACK_IFACE_M=$(SPHEREPACK_IFACE_S:%.f90=mod/%.mod)
SPHEREPACK_S=sphcom.f hrfft.f shaes.f shaec.f shses.f shsec.f shags.f shsgs.f gaqd.f
SPHEREPACK_O=$(SPHEREPACK_S:%.f=obj/%.o)
# dispmodule subroutines
DISPMODULE_S=dispmodule.f90
DISPMODULE_O=$(DISPMODULE_S:%.f90=obj/%.o)
# spline subroutines
SPLINE_S=spline.f90 fd.f90 arrays.f90
SPLINE_O=$(SPLINE_S:%.f90=obj/%.o)

#SOURCES2=bye.f90 module1.f90
#OBJECTS2=$(SOURCES2:.f90=.o)
EXECUTABLE=global_topo

# -----------------------------------------------------------------------------
# (g) Assign dependents to target on dependency line & link options on command
#     line. Command line is initiated with a tab. ($@ is an internal macro for
#     the current target.)

#all: $(SOURCES) $(EXECUTABLE)

#${PROGRAM}: $(OBJECTS)
#
#	$(FC) -o $@ $(OBJECTS) $(MODULES_H_O) $(MODULES_M_O) $(MODULES_L_O) $(SPARSKIT_O) \
#	            ${COMPILER_OPTS} \
#	            $(NETCDF_LIB) \
#	            $(SKIT_LIB) \
#	            $(UMFPACK_LIB)

${PROGRAM}: $(OBJECTS)

	$(FC) $(M32_64) $(OPTS) -I$(MKLROOT)/include \
	$(OBJECTS) \
	$(MODULES_H_O) $(MODULES_M_O) $(MODULES_L_O) $(MYSPARSE_L_O) \
	$(FGMRES_MKL_O) $(SPHEREPACK_O) $(SPHEREPACK_IFACE_O) $(SPARSKIT_O) \
	$(DISPMODULE_O) $(SPLINE_O) $(MISC_L_O) $(MODULES_P_O) \
	$(SKIT_LIB) \
	$(NETCDF_LIB) \
	$(UMFPACK_LIB) \
	$(FPPSTOP) $(S_GRP) \
	$(IFACE_LIB) \
	$(THREADING_LIB) \
	$(CORE_LIB) \
	$(E_GRP) $(OMP_LIB) -lpthread -lm -openmp -o $@

# -----------------------------------------------------------------------------
# Individual compilation using significant suffix. (`make' uses the suffix rule
# .f90.o: f90 -c ... $< to compile .f90 prerequisite files into their .o target
# files if any of the .f90 files were modified since the last make. Same for
# .f.o: f90 -c ... $<.)
#
#.f90.o :
#	f${FORTRAN_VER} -c ${OPTS} $<
#
#.f.o :
#	f${FORTRAN_VER} -c ${OPTS} ${INCLUDE_PATHS} $<
	
obj/igw_v1.o: igw_v1.f90 $(MODULES_H_O) $(MODULES_M_O) $(MODULES_L_O) $(MODULES_P_O) $(FGMRES_MKL_O)
	$(FC) -c ${OPTS} $< -o $@
	
$(MODULES_H_O): obj/%.o: %.f90 $(MODULES_M_O) $(SPARSKIT_O) $(FGMRES_MKL_O)
	$(FC) -c ${OPTS} $(MKL_LIB) $< -o $@
$(MODULES_M_O): obj/%.o: %.f90 $(MODULES_L_O) $(SPHEREPACK_IFACE_O)
	$(FC) -c ${OPTS} $< -o $@
obj/unsym_solvers.o: unsym_solvers.f90 $(MODULES_L_O)
	$(FC) -c ${OPTS} $(MKL_INC) $< -o $@
$(FGMRES_MKL_O): obj/%.o: %.f90 $(MODULES_L_O)
	$(FC) -c ${OPTS} $(MKL_INC) $< -o $@
	
obj/read_etopo.o: read_etopo.f90  $(MODULES_L_O)
	$(FC) -c ${OPTS} $(NETCDF_INC) $< -o $@

$(MODULES_L_O): obj/%.o: %.f90 $(MYSPARSE_L_O) $(SPLINE_O)
	$(FC) -c ${OPTS} $< -o $@
	
$(MYSPARSE_L_O): obj/%.o: %.f90 $(MODULES_P_O) $(DISPMODULE_O) $(MISC_L_O)
	$(FC) -c ${OPTS} $(MKL_INC) $< -o $@

$(MISC_L_O): obj/%.o: %.f90 $(MODULES_P_O) $(DISPMODULE_O) 
	$(FC) -c ${OPTS} $< -o $@
	
$(MODULES_P_O): obj/%.o: %.f90
	$(FC) -c ${OPTS} $< -o $@
	
#-----------------------------------------------------------------------------#

$(SPARSKIT_O): obj/%.o: %.f $(MODULES_P_O)
	$(FC) -c -O2 $(MOD_INC) $< -o $@

# calculations with spherepack routines are always in double precision regardless of wp/cwp
# all that are declared real become real*8
$(SPHEREPACK_IFACE_O): obj/%.o: %.f90 $(SPHEREPACK_O) $(MODULES_P_O) $(MODULES_L_O)
	$(FC) -c -fdefault-real-8 ${OPTS} $< -o $@ 
$(SPHEREPACK_O): obj/%.o: %.f
	$(FC) -c -O2 -fdefault-real-8 $(MOD_INC) $< -o $@ 

$(DISPMODULE_O): obj/%.o: %.f90
	$(FC) -c ${OPTS} $< -o $@
	
$(SPLINE_O): obj/%.o: %.f90
	$(FC) -c ${OPTS} $< -o $@

#umf_test: /home/amt6/amtvl/Desktop/scratch/lib/SuiteSparse/UMFPACK/Demo/umf4zhb.f
#	$(FC) -O -o umf $< $(UMFPACK_LIB)

# -----------------------------------------------------------------------------
# (e, f) SPECIAL TARGET for compiling files containing modules so that
#        module_name.kmo files and associated object files are created or
#        updated before all other source code object files. Also, a target
#        is provided for removing object files and the executable. Invoke
#        this target as make CLEAN. (A tab must proceed the /bin/rm.)
 
#modules:
#	$(FC) -c ${OPTS} MOD_*.f90

test: test.f90
	$(FC) $(OPTS) -I$(MKLROOT)/include \
	$< $(FPPSTOP) $(S_GRP) \
	$(IFACE_LIB) \
	$(THREADING_LIB) \
	$(CORE_LIB) \
	$(E_GRP) $(OMP_LIB) -lpthread -lm -o test
	
clean:
#	rm -f *.o *~ *.mod ${PROGRAM}
	rm -f obj/*.o
	rm -f mod/*.mod		

show:
	echo $(MODULES_L_S) $(MODULES_M_S) $(MODULES_H_S)
	
	
