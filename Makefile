FC     = gfortran-4.2
FFLAGS = -pedantic -O4

#PCLIBS = $(ALGSRC)/libalgencan.a $(DFOSRC)/libdfo.a
PCLIBS = 
#LIBBLAS = /usr/share/lapack-3.2.1/blas_LINUX.a
#LIBLAPACK = /usr/share/lapack-3.2.1/lapack_LINUX.a

all:

eng-toyprob: toyprob.o runeng.o skinny.o engdata.o extremebarrier.o \
	     restoration.o sds.o drandsc.o bobyqadata.o

	$(FC) $(FFLAGS) $^ $(PCLIBS) -o $@

runeng.o: runeng.f90 skinny.mod
	$(FC) $(FFLAGS) -c runeng.f90

skinny.o: skinny.f90 engdata.mod #nelder_mead.mod
	$(FC) $(FFLAGS)  -c skinny.f90

restoration.o: restoration.f90 engdata.mod bobyqadata.mod
	$(FC) $(FFLAGS) -c restoration.f90

extremebarrier.o: extremebarrier.f90 engdata.mod
	$(FC) $(FFLAGS) -c extremebarrier.f90

#mintr_algencan.o: mintr_algencan.f
#	$(FC) $(FFLAGS) -c mintr_algencan.f

%.mod: %.o
	

%.o: %.f90
	$(FC) $(FFLAGS) -c $^

clean:
	rm *.o *.mod -vf
