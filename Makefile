### make changes accordingly ###
CC       = gcc
CPP      = g++
CLINKER  = gcc
CCLINKER = g++
MAKE     = make --no-print-directory
SHELL    = /bin/sh
CFLAGS   = -Wall -pedantic 
#OPTI     = -p
OPTI     = -O2
#OPTI     = -g
LDFLAGS	= 
INCLUDES	= 
LIBS	= -lm -lgsl -lgslcblas
DEFINES = -DVERBOSE 

default: model

model: $(OBJS) Makefile Person.o Location.o Mosquito.o Community.o driver.o Parameters.o 
	$(CCLINKER) $(OPTI) -o model Person.o Location.o Mosquito.o Community.o driver.o Parameters.o $(OBJS) $(LDFLAGS) $(LIBS)

%.o: %.cpp Parameters.h Person.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c $<

zip: *.cpp *.h Makefile README LICENSE.txt HISTORY.txt *-bangphae.txt
	cd ..; zip denguemodelcode/denguemodel.zip denguemodelcode/README denguemodelcode/LICENSE.txt denguemodelcode/HISTORY.txt denguemodelcode/gpl.txt denguemodelcode/Makefile denguemodelcode/*.cpp denguemodelcode/*.h denguemodelcode/*-bangphae.txt

emacs:
	emacs Makefile *.h *.cpp README LICENSE.txt HISTORY.txt &

clean:
	rm -f *.o model *~
