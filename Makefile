SHELL=/bin/bash
G++VER := $(shell command -v g++-4.9)

ifndef G++VER
CPP:=g++
else
CPP:=g++-4.9
endif

MAKE     	= make --no-print-directory
CFLAGS   	= -Wall -Wextra -pedantic -std=c++11
#OPTI     	= -g
OPTI     	= -O2
LDFLAGS	 	= -L $(HPC_GSL_LIB) $(TACC_GSL_LIB)
INCLUDES 	= -I $(HPC_GSL_INC) $(TACC_GSL_INC)
LIBS     	= -lm -lgsl -lgslcblas
DEFINES  	= -DVERBOSE 

default: model

model: $(OBJS) Makefile simulator.h Person.o Location.o Mosquito.o Community.o driver.o Parameters.o Utility.o
	$(CPP) $(CFLAGS) $(OPTI) -o model Person.o Location.o Mosquito.o Community.o driver.o Parameters.o Utility.o $(OBJS) $(LDFLAGS) $(LIBS)

%.o: %.cpp Community.h Location.h Mosquito.h Utility.h Parameters.h Person.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c $<

clean:
	rm -f *.o model *~
