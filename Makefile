### make changes accordingly ###
CC       = gcc
CPP      = g++
CLINKER  = gcc
CCLINKER = g++
MAKE     = make --no-print-directory
SHELL    = /bin/sh
CFLAGS		= -Wall -pedantic 
OPTI            = -O3
LDFLAGS	= 
INCLUDES	= 
LIBS	= -lm -lgsl -lgslcblas
DEFINES = -DVERBOSE 

default: model

model: $(OBJS) Makefile Person.o Location.o Mosquito.o Community.o driver.o
	$(CCLINKER) -o model Person.o Location.o Mosquito.o Community.o driver.o $(OBJS) $(LDFLAGS) $(LIBS)

%.o: %.cpp Person.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c $<

zip: *.cpp *.h Makefile README LICENSE.txt HISTORY.txt *-bangphae.txt
	cd ..; zip denguemodelcode/denguemodel.zip denguemodelcode/README denguemodelcode/LICENSE.txt denguemodelcode/HISTORY.txt denguemodelcode/gpl.txt denguemodelcode/Makefile denguemodelcode/*.cpp denguemodelcode/*.h denguemodelcode/*-bangphae.txt

emacs:
	emacs Makefile *.h *.cpp README LICENSE.txt HISTORY.txt &

clean:
	rm -f *.o model *~
