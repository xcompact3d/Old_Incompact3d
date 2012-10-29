#=======================================================================
# Makefile for Imcompact3D
#=======================================================================

# Choose pre-processing options
#   -DSHM	   - enable shared-memory implementation
#   -DDOUBLE_PREC  - use double-precision (default single)
OPTIONS =

# Choose an FFT engine, available options are:
#   fftw3      - FFTW version 3.x
#   generic    - A general FFT algorithm (no 3rd-party library needed)
FFT=generic

# Paths to FFTW 3
FFTW3_PATH=/home/lining/software/fftw-alpha
FFTW3_INCLUDE = -I$(FFTW3_PATH)/include
FFTW3_LIB = -L$(FFTW3_PATH)/lib -lfftw3

# include PATH 
ifeq ($(FFT),generic)
  INC=
else ifeq ($(FFT),fftw3)
  INC=$(FFTW3_INCLUDE)
endif

# library path
ifeq ($(FFT),generic)
   LIBFFT=
else ifeq ($(FFT),fftw3)
   LIBFFT=$(FFTW3_LIB)
endif

# Compiler 

# PGI
#FC = ftn
#OPTFC = -O3 -fast -Mpreprocess
#CC = cc
#CFLAGS = -O3

# PathScale
#FC = ftn
#OPTFC = -Ofast -cpp
#CC = cc
#CFLAGS = -O3

# GNU
#PREP=/home/lining/build/scalasca/bin/scalasca -inst
#PREP=
FC = $(PREP) mpif90
OPTFC = -g -cpp
CC = mpicc
CFLAGS = -g

#Cray
#FC = ftn
#OPTFC = -O3 -F
#CC = cc
#CFLAGS = -O3

# NAG
#FC = f95
#OPTFC = -O3 -kind=byte -I/usr/local/packages/nag/NAGWare5.1_496/NAGWare_f95-amd64/lib

# List of source files
SRC = decomp_2d.f90 glassman.f90 fft_$(FFT).f90 module_param.f90 io.f90 variables.f90 poisson.f90 schemes.f90 convdiff.f90 incompact3d.f90 navier.f90 derive.f90 parameters.f90 tools.f90 visu.f90

#-----------------------------------------------------------------------
# Normally no need to change anything below

ifneq (,$(findstring DSHM,$(OPTIONS)))   
OBJ =	$(SRC:.f90=.o) alloc_shm.o
else
OBJ =	$(SRC:.f90=.o)
endif	

all: incompact3d

alloc_shm.o: alloc_shm.c
	$(CC) $(CFLAGS) -c $<

incompact3d : $(OBJ)
	$(FC) -O3 -o $@ $(OBJ) $(LIBFFT)

%.o : %.f90
	$(FC) $(OPTFC) $(OPTIONS) $(INC) -c $<

.PHONY: clean 
clean:
	rm -f *.o *.mod incompact3d

.PHONY: realclean
realclean: clean
	rm -f *~ \#*\#
