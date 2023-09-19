MAKE     	= make --no-print-directory
CFLAGS   	= -Wall -Wextra -pedantic -std=c++17
#OPTI     	= -g
OPTI     	= -O2
LIBS     	= -lm -lgsl -lgslcblas
DEFINES  	= -DVERBOSE 

model: Makefile simulator.h Person.o Location.o Mosquito.o Community.o Parameters.o Utility.o

%.o: %.cpp Community.h Location.h Mosquito.h Utility.h Parameters.h Person.h Makefile
	$(CXX) $(CFLAGS) $(OPTI) $(DEFINES) -c $<

clean:
	rm -f *.o
