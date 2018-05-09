CXXFLAGS      = -O2 
#LDFLAGS       =  -lm -ldl
LINK = -lstdc++

CC=g++
CFLAGS=-c -Wall -O3        

SOURCES= DM.o crossection.o main.o


csect: $(SOURCES)
	$(CC) -o csect $(SOURCES) $(LDFLAGS)

DM.o: DM.h DM.cc
	$(CC) $(CFLAGS)  DM.cc   -o $@

crossection.o:   crossection.h constants.h  util2.h DM.h ../jednostki.h crossection.cc
	$(CC) $(CFLAGS) crossection.cc   -o $@
	

main.o:  crossection.h  main.cpp
	$(CC) $(CFLAGS)  main.cpp -o $@
	
clean:
	rm -rf *.o  csect  


