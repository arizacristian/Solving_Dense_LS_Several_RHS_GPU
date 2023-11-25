# This is a simple standalone example. See README.txt
# Initially it is setup to use OpenBLAS.
# See magma/make.inc for alternate BLAS and LAPACK libraries,
# or use pkg-config as described below.
##
#  watch -n0.5 nvidia-smi    look use of nvidia
##  set nthead   
#		export OMP_NUM_THREADS=6
#                      OPENBLAS_NUM_THREADS
##


#-include make.check-mkl
#-include make.check-cuda

SYSTEM = marreca_gpu
# Paths where MAGMA, CUDA, and OpenBLAS are installed.
# MAGMADIR can be .. to test without installing.
#MAGMADIR     ?= ..
#MAGMADIR     ?= /usr/local/magma
#CUDADIR      ?= /usr/local/cuda
#OPENBLASDIR  ?= /usr/local/openblas

##   export

MACHINE=marreca_mkl
ifeq ($(MACHINE),marreca_openblas)
##   marreca
OPENBLASDIR ?= /home/cdaa/user2/OpenBLAS
MAGMADIR     ?= /home/cdaa/user2/magma
CUDADIR ?= /usr/local/cuda-10.0

else ifeq ($(MACHINE),ogun_openblas)

#OPENBLASDIR ?= /home/milton.porsani/user/OpenBLAS
#MAGMADIR     ?= /home/milton.porsani/user/magma
#CUDADIR ?= /usr/local/cuda-9.2/
OPENBLASDIR ?= /scratch/milton.porsani/user/OpenBLAS
MAGMADIR     ?= /scratch/milton.porsani/user/magma
CUDADIR ?= /usr/local/cuda-9.2/

else ifeq ($(MACHINE),ogun_mkl)


#MKLDIR ?= /opt/soft/intel/compilers_and_libraries_2018.1.163/linux/mkl
MKLDIR ?= $(MKLROOT)
#MKLDIR ?=/opt/soft/intel/2019/compilers_and_libraries_2019.0.117/linux/mkl
#MKLDIR = /exports/soft/intel/compilers_and_libraries_2018.1.163/linux/mkl
MAGMADIR     ?= /scratch/milton.porsani/user/magma_mkl
CUDADIR ?= /usr/local/cuda-9.2/


else ifeq ($(MACHINE),marreca_mkl)

MKLDIR ?= /home/cdaa/intel2019/compilers_and_libraries_2019.5.281/linux/mkl
#MAGMADIR     ?= /home/cdaa/user2/magma_mkl_icpc
MAGMADIR     ?=/home/cdaa/user2/magma-2-5-3
CUDADIR ?= /usr/local/cuda-10.0

endif


##  check mkl and cuda
-include make.check-mkl
-include make.check-cuda


ifeq ($(MACHINE),marreca_mkl)

CC        = icc
CXX       = icpc
NVCC      = nvcc
FORT      = ifort
#FORT      = pgf90

else ifeq ($(MACHINE),ogun_mkl)

CC        = icc
CXX       = icpc
NVCC      = nvcc
FORT      = ifort

else

CC            = gcc
CXX	      = g++
FORT          = gfortran
#LD            = gcc

endif

libName ?=mkl

# compiler= gnu, intel
compiler ?=intel
ifeq ($(compiler),gnu)
	F90 = gfortran
	MAGMADIR     =/home/cdaa/user2/magma-2-5-3_gnu_mkl
else ifeq ($(compiler),intel)
	F90 = ifort
	MAGMADIR     =/home/cdaa/user2/magma-2-5-3
endif





linking ?=dynamic
ifeq ($(linking),dynamic)
    lstr = lmkl
else ifeq ($(linking),static)
    lstr = libmkl
endif

interface ?=ilp64
ifeq ($(compiler),gnu)
    lib_mkl = $(lstr)_gf_$(interface)
else ifeq ($(compiler),intel)
    lib_mkl = $(lstr)_intel_$(interface)
endif


    ifeq ($(compiler),gnu)
        lib_omp = $(lstr)_gnu_thread
        Compiler_opt +=  -fopenmp
    else ifeq ($(compiler),intel)
        lib_omp = $(lstr)_intel_thread
        Compiler_opt += -qopenmp
    endif
  
#$(info lib omp= $(lib_omp))  

    ifeq ($(compiler),gnu)
	        Compiler_opt += -fdefault-integer-8 -m64
    else ifeq ($(compiler),intel)
	        Compiler_opt += -integer-size 64
    endif



# Use -fPIC to make shared (.so) and static (.a) library;
# can be commented out if making only static library.
FPIC      = -fPIC
#LD	 ?= -O3 $(FPIC)

#################################################### compiler flags
F90OPT = -c -O3 $(Compiler_opt)

LFOPT = -O3 $(Compiler_opt)



# ----- include paths
MAGMA_INC  = -I$(MAGMADIR)/include
MAGMA_INC += -I./control
MAGMA_INC += -I./testing


CFLAGS        = -Wall
INC       = -I$(CUDADIR)/include \
            -I$(MKLROOT)/include
CPPFLAGS   = $(INC) $(MAGMA_INC)

# may want -std=c99 for CFLAGS, -std=c++11 for CXXFLAGS
CFLAGS     ?= -O3 $(FPIC) -DADD_ -Wall -fopenmp -std=c++11 -DHAVE_CUBLAS -DMIN_CUDA_ARCH=300
CXXFLAGS   ?= $(CFLAGS) -std=c++11
#F90FLAGS  = -O3 $(FPIC) -Wall -Wno-unused-dummy-argument -x f95-cpp-input



############  DEBUG
#ifor  intel
#FFLAGS = -g -CB -ftrapuv -fp-model=strict -fpe-all=0 -traceback -warn
#GNU  GFORTRAN
#FFLAGS = -g –fbounds-check –ffpe-trap= -finit-real=snan -Wall
## pgi  pgfortran
#FFLAGS = -g –C –Ktrap=fp –Minform=inform

#################################################### debug flags
# -CB checks array bounds: very slow!  -fcheck=bounds  -fcheck=all   -Wno-tabs
ifeq ($(compiler),gnu)
	DEBUG_opt= -O0 -g3 -fcheck=bounds -fcheck=all -ffpe-trap= -finit-real=snan -Wall -Wno-tabs
    #DEBUG_opt += -fno-stack-arrays #-Ofast -O3
	F90DEBUG =   $(DEBUG_opt) $(Compiler_opt)
    LFDEBUG  =   $(DEBUG_opt) $(Compiler_opt)
else ifeq ($(compiler),intel)
	## mpiifort debug options
	# -O0       disable optimizations ; -g[level] 3  - Emit extra information which may be useful for some tools.
	DEBUG_opt = -O0 -g3
	DEBUG_opt += -traceback -check all -fp-stack-check
	#  -CB: -check bounds; -ftrapuv  trap uninitialized variables;  Check floating point exception ; -warn level of warning messages
	DEBUG_opt += -CB -ftrapuv -fp-model strict -fpe-all=0 -no-ftz  -warn all  
	##        += -g -CB -ftrapuv -fp-model=strict -fpe-all=0 -traceback -warn       
	#DEBUG_opt += -debug extended
	#DEBUG_opt += -fno-merge-debug-strings -fno-merge-constants
	 #-print-multi-lib :print information about libraries being used
	#DEBUG_opt += -print-multi-lib
	# -heap-arrays [n] temporary arrays of minimum size n (in kilobytes) are allocated in heap memory rather than on the stack.  
	#                  If n is not specified, all temporary arrays are allocated in heap memory.
    #DEBUG_opt += -heap-arrays
          
        #DEBUG_opt = -O0 -g3 -debug all -check all -check bounds -traceback -warn

	F90DEBUG =   $(DEBUG_opt) $(Compiler_opt)
	LFDEBUG  =   $(DEBUG_opt) $(Compiler_opt)

endif


#ifeq ($(MACHINE),marreca_mkl)
##  ###      gfortran
##F90FLAGS  = -O3 $(FPIC)          -DNDEBUG -DADD_ -Wall -Wno-unused-dummy-argument -x f95-cpp-input
##F90FLAGS  += -fdefault-integer-8
##LDFLAGS   =     $(FPIC) -fopenmp

##      ifort
##F90FLAGS  = -O3 $(FPIC)          -DNDEBUG -DADD_ -warn all -warn nounused
##F90FLAGS  += -integer-size 64
##LDFLAGS   =     $(FPIC) -qopenmp
#endif



#use the debug flags if $(DEBUG)='yes':
ifeq ($(DEBUG),yes)
   F90FLAGS = $(F90DEBUG)
   LDFLAGS = $(LFDEBUG)
else
   F90FLAGS += $(F90OPT)
   LDFLAGS += $(LFOPT)
endif





########################################################## PGO flags
ifeq ($(MACHINE),marreca_mkl)
#use the pgo flags if $(prof_gen)='yes':
ifeq ($(prof_gen),yes)
	F90FLAGS      +=-prof_gen
endif
ifeq ($(prof_use),yes)
	F90FLAGS      +=-prof_use
endif
endif



# needs -fopenmp if MAGMA was compiled with OpenMP
#ifeq ($(MACHINE),marreca_mkl)
#LDFLAGS       = $(FPIC)  -qopenmp -W1
#else ifeq ($(MACHINE),ogun_mkl)
#LDFLAGS       = $(FPIC)  -qopenmp -W1
#else
#LDFLAGS       = $(FPIC)  -fopenmp -Wall
#endif



# Extension for object files: o for unix, obj for Windows?
o_ext      ?= o

# ------------------------------------------------------------------------------
# Define the pointer size for fortran compilation
# If there's an issue compiling sizeptr, assume 8 byte (64 bit) pointers
PTRFILE = control/sizeptr.c
PTROBJ  = control/sizeptr.o
PTREXEC = control/sizeptr
PTRSIZE = $(shell if [ -x $(PTREXEC) ]; then $(PTREXEC); else echo 8; fi)
PTROPT  = -Dmagma_devptr_t="integer(kind=$(PTRSIZE))"



$(PTREXEC): $(PTROBJ)
	-$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<
	touch $@

$(PTROBJ): $(PTRFILE)
	-$(CC) $(CFLAGS) -c -o $@ $<
	touch $@

###$(F90) $(F90FLAGS) $(CPPFLAGS) $(MAGMA_F90FLAGS) -c -o $@ $<

# ----------------------------------------
# Flags and paths to MAGMA, CUDA, and LAPACK/BLAS
 MAGMA_CFLAGS     := -DADD_ \
                     -I$(MAGMADIR)/include \
                     -I$(CUDADIR)/include \
		     -I./testing -ltest
#                      -I$(MAGMADIR)/sparse/include \
 MAGMA_F90FLAGS   := -Dmagma_devptr_t="integer(kind=8)" \
                     -I$(MAGMADIR)/include

# MAGMA_F90FLAGS   := -I$(CUDADIR)/include \
#		     -I$(MAGMADIR)/include \
#                     -I./testing 

# 
# 

# libraries MKL
#LIB       = -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lstdc++ -lm  ###  -liomp5 -lpthread -lm -ldl
###LIB       = -lmkl_gf_ilp64 -lmkl_gnu_thread -lmkl_core -lpthread -lstdc++ -lm -lgfortran  ###ilp64 and gcc
ifeq ($(compiler),gnu)
    LIB       = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -$(lib_mkl) -$(lib_omp) -lmkl_core -lgomp -lpthread -lm -ldl
else ifeq ($(compiler),intel)
    #LIB       = -$(lib_mkl) -$(lib_omp) -lmkl_core -lpthread -lstdc++ -lm -ldl
    LIB       = -$(lib_mkl) -$(lib_omp) -lmkl_core -lpthread  -lm -ldl
endif



# libraries cublas
LIB      += -lcublas -lcusparse -lcudart -lcudadevrt




# # may be lib instead of lib64 on some systems
ifeq ($(MACHINE),marreca_mkl)
 MAGMA_LIBS       := -L./testing -ltest \
		     -L./testing/lin -llapacktest \
		     -L$(MAGMADIR)/lib -lmagma \
                     -L$(CUDADIR)/lib64  \
                     -L$(MKLDIR)/lib/intel64 $(LIB)

else ifeq ($(MACHINE),ogun_mkl)
 MAGMA_LIBS       := -L./testing -ltest \
		     -L./testing/lin -llapacktest \
		     -L$(MAGMADIR)/lib -lmagma \
                     -L$(CUDADIR)/lib64  \
                     -L$(MKLDIR)/lib/intel64 $(LIB)

else
 MAGMA_LIBS       := -L./testing -ltest \
		     -L./testing/lin -llapacktest \
		     -L$(MAGMADIR)/lib -lmagma \
                     -L$(CUDADIR)/lib64 -lcublas -lcusparse -lcudart -lcudadevrt \
                     -L$(OPENBLASDIR)/lib -lopenblas  
endif

#gfortran -fPIC                       -fopenmp -Wl,-rpath,/home/cdaa/Lib/magma-2.5.1/lib \
#	-o testing/testing_dgetrf_f testing/testing_dgetrf_f.o \
#	-L./testing -ltest \
#	-L./testing/lin -llapacktest \
#	-L./lib -lmagma \
#	-L/usr/local/cuda-10.0/lib64 -L/home/cdaa/user/OpenBLAS/lib -lopenblas -lcublas -lcusparse -lcudart -lcudadevrt

#g++ -fPIC                       -fopenmp -Wl,-rpath,/home/cdaa/Lib/magma-2.5.1/lib \
#	-o testing/testing_sgemm testing/testing_sgemm.o \
#	-L./testing -ltest \
#	-L./lib -lmagma \
#	-L./testing/lin -llapacktest \
#	-L/usr/local/cuda-10.0/lib64 -L/home/cdaa/user/OpenBLAS/lib -lopenblas -lcublas -lcusparse -lcudart -lcudadevrt

# ----------------------------------------
# Alternatively, using pkg-config (see README.txt):
#MAGMA_CFLAGS   := $(shell pkg-config --cflags magma)

#MAGMA_F90FLAGS := -Dmagma_devptr_t="integer(kind=8)" \
#                  $(shell pkg-config --cflags-only-I magma)

#MAGMA_LIBS     := $(shell pkg-config --libs   magma)

# ----------------------------------------
default:
	@echo "Available make targets are:"
	@echo "  make default    # print this help "
	@echo "  make all        # compiles gemm posv posv_gpu potrf potrf_mgpu "
	@echo "  make gemm       # compiles testing_dgemm"
	@echo "  make gemm_gpu   # compiles testing_sgemm_gpu_f testing_dgemm_gpu_f"
#	@echo "  make posv       # compiles testing_sposv 	testing_dposv"
#	@echo "  make posv_gpu   # compiles testing_sposv_gpu 	testing_dposv_gpu"
	@echo "  make getrf      # compiles testing_sgetrf_f testing_dgetrf_f testing_cgetrf_f testing_zgetrf_f"
	@echo "  make getrf_gpu  # compiles testing_sgetrf_gpu_f testing_dgetrf_gpu_f testing_cgetrf_gpu_f testing_zgetrf_gpu_f"
	@echo "  make potrf      # compiles testing_spotrf_f 	testing_dpotrf_f"
	@echo "  make potrf_gpu  # compiles testing_spotrf_gpu_f testing_dpotrf_gpu_f driver_dpotrf_trs_gpu_f"
#	@echo "  make potrf_mgpu # compiles testing_spotrf_mgpu testing_dpotrf_mgpu "
	@echo "  make clean      # deletes executables and object files"
	@echo "                     "
	@echo "                    testing   (not in all)  "
	@echo "  make testing_BIG_dgemm   # compiles  testing_BIG_dgemm  "
	@echo "  make test_magmaHybrid    # compiles  test_magmaHybrid  "
	@echo "                  PLS  "	
	@echo "  make PLS2019_magma_d_v6       #compiles PLS2019_magma_d_v6 "
	@echo "  make PLS2019_magma_d_v6       #compiles PLS2019_magma_d_v6 "
	@echo "                     "
	@echo "                  driver magma compiler flag "
	@echo "  make driver_dpotrf_trs_f compiler=gnu      # compiles driver_dpotrf_trs_f  compiler=<gnu;intel> "
	@echo "  make test_driver_dpotrf_trs_f compiler=gnu # compiles test_driver_dpotrf_trs_f  compiler=<gnu;intel>"
	@echo "                  driver PLS   compiler flag "
	@echo "  make PLS2019_magma_d_v5 compiler=gnu       # compiles PLS2019_magma_d_v5  compiler=<gnu;intel> "
	@echo "  make test_PLS2019_magma_d_v5 compiler=gnu  # compiles test_PLS2019_magma_d_v5  compiler=<gnu;intel>"
		
all: gemm gemm_gpu getrf potrf_mgpu getrf_gpu
##posv posv_gpu


gemm:  testing_dgemm_f
gemm_gpu: testing_sgemm_gpu_f testing_dgemm_gpu_f
posv: testing_sposv testing_dposv
posv_gpu: testing_sposv_gpu testing_dposv_gpu
getrf: testing_sgetrf_f testing_dgetrf_f testing_cgetrf_f testing_zgetrf_f
getrf_gpu: testing_dgetrf_gpu_f testing_sgetrf_gpu_f testing_cgetrf_gpu_f testing_zgetrf_gpu_f
potrf:  testing_dpotrf_f testing_spotrf_f
potrf_gpu: testing_spotrf_gpu_f testing_dpotrf_gpu_f driver_dpotrf_trs_gpu_f

PLS: PLS2019_magma_d_v5
driver: driver_dpotrf_trs_f

##################################################### objects alpha
INIT = .
SRCDIR = $(INIT)
OBJ=$(INIT)
OBJDIR = $(OBJ)
CALCDIR =$(INIT)


clean:
	-rm -fv *.o *.mod *__genmod.f90

.SUFFIXES:


# ----------------------------------------
# C example
#%.o: %.cpp
#	$(CPP) $(CFLAGS) $(MAGMA_CFLAGS) -c -o $@ $<

# f90 
%.$(o_ext): %.f90
	@echo
	@echo $@ Linking...
	$(F90) $(F90FLAGS) $(CPPFLAGS) $(MAGMA_F90FLAGS) -c -o $@ $<

%.$(o_ext): %.F90 $(PTREXEC)
	@echo
	@echo $@ Linking...
	$(F90) $(F90FLAGS) $(CPPFLAGS) $(PTROPT) -c -o $@ $<

#%.$(o_ext): %.f90
#	$(F90) $(F90FLAGS) $(CPPFLAGS) -c -o $@ $<
#	-mv $(notdir $(basename $@)).mod include/

#%.$(o_ext): %.F90 $(PTREXEC)
#	$(F90) $(F90FLAGS) $(CPPFLAGS) $(PTROPT) -c -o $@ $<
#	-mv $(notdir $(basename $@)).mod include/


#gfortran -O3 -fPIC -DNDEBUG -DADD_ -Wall -Wno-unused-dummy-argument -x f95-cpp-input -I/usr/local/cuda-10.0/include -I./include -I./testing -c -o testing/testing_dgetrf_f.o testing/testing_dgetrf_f.f90

#ifort -O3 -fPIC          -DNDEBUG -DADD_ -warn all -warn nounused -integer-size 64 -I/usr/local/cuda-10.0/include -I/home/cdaa/intel2019/compilers_and_libraries_2019.5.281/linux/mkl/include -I./include -I./testing -c -o testing/testing_dgetrf_f.o testing/testing_dgetrf_f.f90

##   GEMM_gpu_f  

testing_sgemm_gpu_f: printmsgA testing_sgemm_gpu_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.

testing_dgemm_gpu_f: printmsgA testing_dgemm_gpu_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.




testing_dgemm_f: printmsgA testing_dgemm_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.

testing_BIG_dgemm: printmsgA Sub_CPU_Interfaces.o testing_BIG_dgemm.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(word 3,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.

## Xgetrf_f

testing_sgetrf_f: printmsgA testing_sgetrf_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.

#$(testers_f): %: %.$(o_ext)
#	$(F90) $(LDFLAGS) $(RPATH) \
#	-o $@ $< \
#	-L./testing -ltest \
#	-L./testing/lin -llapacktest \
#	-L./lib -lmagma \
#	$(LIBS)

#gfortran -fPIC                       -fopenmp -Wl,-rpath,/home/cdaa/Lib/magma-2.5.1/lib \
#	-o testing/testing_dgetrf_f testing/testing_dgetrf_f.o \
#	-L./testing -ltest \
#	-L./testing/lin -llapacktest \
#	-L./lib -lmagma \
#	-L/usr/local/cuda-10.0/lib64 -L/home/cdaa/user/OpenBLAS/lib -lopenblas -lcublas -lcusparse -lcudart -lcudadevrt

testing_dgetrf_f: printmsgA testing_dgetrf_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.
## asi compila el original
#gfortran -O3 -fPIC -DNDEBUG -DADD_ -Wall -Wno-unused-dummy-argument -x f95-cpp-input -I/usr/local/cuda-10.0/include -I./include -I./testing -c -o testing/testing_dgetrf_f.o testing/testing_dgetrf_f.f90
##mv testing_dgetrf_f.mod include/
##mv: cannot stat `testing_dgetrf_f.mod': No such file or directory
##make[1]: [testing/testing_dgetrf_f.o] Error 1 (ignored)
#gfortran -fPIC                       -fopenmp -Wl,-rpath,/home/cdaa/Lib/magma-2.5.1/lib \
#	-o testing/testing_dgetrf_f testing/testing_dgetrf_f.o \
#	-L./testing -ltest \
#	-L./testing/lin -llapacktest \
#	-L./lib -lmagma \
#	-L/usr/local/cuda-10.0/lib64 -L/home/cdaa/user2/OpenBLAS/lib -lopenblas -lcublas -lcusparse -lcudart -lcudadevrt

#ifort -fPIC -qopenmp -Wl,-rpath,/home/cdaa/Lib/magma-2.5.1/lib \
#	-o testing/testing_dgetrf_f testing/testing_dgetrf_f.o \
#	-L./testing -ltest \
#	-L./testing/lin -llapacktest \
#	-L./lib -lmagma \
#	-L/usr/local/cuda-10.0/lib64 -L/home/cdaa/intel2019/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64 -lmkl_intel_ilp64 -# lmkl_intel_thread -lmkl_core -lpthread -lstdc++ -lm -lcublas -lcusparse -lcudart -lcudadevrt


testing_cgetrf_f: printmsgA testing_cgetrf_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.

testing_zgetrf_f: printmsgA testing_zgetrf_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.



##   Xpotrf

testing_spotrf_f: printmsgA testing_spotrf_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.



testing_dpotrf_f: printmsgA Sub_CPU_Interfaces.o testing_dpotrf_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(word 3,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.

driver_dpotrf_trs_f: printmsgA Sub_CPU_Interfaces.o driver_dpotrf_trs_f.o
	-rm -f $@_${compiler}_${libName}.exe
	@echo
	@echo $@_${compiler}_${libName}.exe Linking...
	$(F90) $(LDFLAGS) -o $@_${compiler}_${libName}.exe $(word 2,$^) $(word 3,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@_${compiler}_${libName}.exe built, enjoy.
	-rm -fv *.o *.mod *__genmod.f90
	
test_driver_dpotrf_trs_f: printmsgA Sub_Utils.o Sub_CPU_Interfaces.o test_driver_dpotrf_trs_f.o
	-rm -f $@_${compiler}_${libName}.exe
	@echo
	@echo $@_${compiler}_${libName}.exe Linking...
	$(F90) $(LDFLAGS) -o $@_${compiler}_${libName}.exe $(word 2,$^) $(word 3,$^) $(word 4,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@_${compiler}_${libName}.exe built, enjoy.
	-rm -fv *.o *.mod *__genmod.f90


###   PLS magma


PLS2019_magma_d_v5: printmsgA Sub_CPU_Interfaces.o mod_IO_S.o SBHSNE2019_magma_d_v5.o PLS2019_magma_d_v5.o
	-rm -f $@_${compiler}_${libName}.exe 
	@echo
	@echo $@_${compiler}_${libName}.exe Linking...
	$(F90) $(LDFLAGS) -o $@_${compiler}_${libName}.exe $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@_${compiler}_${libName}.exe built, enjoy.
	-rm -fv *.o *.mod *__genmod.f90
	
	
test_PLS2019_magma_d_v5: printmsgA Sub_Utils.o Sub_CPU_Interfaces.o mod_IO_S.o SBHSNE2019_magma_d_v5.o test_PLS2019_magma_d_v5.o
	-rm -f $@_${compiler}_${libName}.exe 
	@echo
	@echo $@_${compiler}_${libName}.exe Linking...
	$(F90) $(LDFLAGS) -o $@_${compiler}_${libName}.exe $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) $(word 6,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@_${compiler}_${libName}.exe built, enjoy.
	-rm -fv *.o *.mod *__genmod.f90

PLS2019_magma_d_v6: printmsgA Sub_CPU_Interfaces.o SBHSNE2019_magma_d_v6.o PLS2019_magma_d_v6.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(word 3,$^) $(word 4,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.


test_magmaHybrid: printmsgA  mod_test_magma_hybrid.o test_magmaHybrid.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(word 3,$^)  $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.


##   getrf_gpu_f   funcionando
testing_sgetrf_gpu_f: printmsgA testing_sgetrf_gpu_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.

testing_dgetrf_gpu_f: printmsgA testing_dgetrf_gpu_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.


##  asi compilo el original
#gfortran -O3 -fPIC -DNDEBUG -DADD_ -Wall -Wno-unused-dummy-argument -x f95-cpp-input -I/usr/local/cuda-10.0/include -I./include -I./testing -Dmagma_devptr_t="integer(kind=8)" -c -o testing/testing_dgetrf_gpu_f.o testing/testing_dgetrf_gpu_f.F90
##mv testing_dgetrf_gpu_f.mod include/
##mv: cannot stat `testing_dgetrf_gpu_f.mod': No such file or directory
##make[1]: [testing/testing_dgetrf_gpu_f.o] Error 1 (ignored)

#gfortran -fPIC                       -fopenmp -Wl,-rpath,/home/cdaa/Lib/magma-2.5.1/lib \
#	-o testing/testing_dgetrf_gpu_f testing/testing_dgetrf_gpu_f.o \
#	-L./testing -ltest \
#	-L./testing/lin -llapacktest \
#	-L./lib -lmagma \
#	-L/usr/local/cuda-10.0/lib64 -L/home/cdaa/user2/OpenBLAS/lib -lopenblas -lcublas -lcusparse -lcudart -lcudadevrt


testing_cgetrf_gpu_f: printmsgA testing_cgetrf_gpu_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.

testing_zgetrf_gpu_f: printmsgA testing_zgetrf_gpu_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.



##   POtrf_gpu_f   
testing_spotrf_gpu_f: printmsgA testing_spotrf_gpu_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.

testing_dpotrf_gpu_f: printmsgA testing_dpotrf_gpu_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.

driver_dpotrf_trs_gpu_f: printmsgA driver_dpotrf_trs_gpu_f.o
	-rm -f $@ 
	@echo
	@echo $@ Linking...
	$(F90) $(LDFLAGS) -o $@ $(word 2,$^) $(MAGMA_LIBS)
	@echo "===============>>>>>>>>>>>   "$@ built, enjoy.
	-rm -fv *.o *.mod *__genmod
	





########################################################## message
printmsgA :
	@echo
	@echo Building POST_PROCESSOR for $(SYSTEM)
	@echo Compiler flags  : $(CXXFLAGS)
	@echo Compiler        : $(F90)
	@echo MAGMA   flags   : $(MAGMA_CFLAGS)
	@echo cpp flags       : $(CPPFLAGS)
#
###################################################### end of file

