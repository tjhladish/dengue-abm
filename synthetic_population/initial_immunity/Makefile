CPP      = g++
#CPP      = mpicxx
MAKE     = make --no-print-directory
CFLAGS   = -Wall -pedantic -std=c++11
OPTI     = -O2
LDFLAGS	 = -L $(HPC_GSL_LIB) $(TACC_GSL_LIB)
INCLUDES = -I.

ifdef TACC_GSL_INC
INCLUDES += -I$$TACC_GSL_INC
endif
ifdef HPC_GSL_INC
INCLUDES += -I$$HPC_GSL_INC
endif

LIBS     = -lm -lgsl -lgslcblas

all_mpi: ImmunityGenerator.cpp ImmunityGenerator.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) -c ImmunityGenerator.cpp

all_no_mpi: ImmunityGenerator.cpp ImmunityGenerator.h Makefile
	g++ $(CFLAGS) $(OPTI) $(INCLUDES) -c ImmunityGenerator.cpp

icc_immgen: main.cpp immunity_generator.cpp
	icc $(CFLAGS) $(OPTI) main.cpp immunity_generator.cpp -o immgen $(INCLUDES) $(LDFLAGS) -lm -lgsl -lgslcblas

immgen: yucatan_immunity_generator.cpp
	g++ $(CFLAGS) $(OPTI) yucatan_immunity_generator.cpp -o immgen -lm -lgsl -lgslcblas

clean:
	rm -f *.o 
