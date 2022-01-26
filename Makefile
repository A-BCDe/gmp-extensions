CC = g++
CFLAGS = -Wall
LDFLAGS = -lgmp -lgmpxx
LIB = extensions
INCLUDE = -I$(LIB)/include/

all:	polynomial matrix

polynomial:	app/polynomial.cc lib
	$(CC) app/polynomial.cc $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o polynomial polynomial.o matrix.o

lib:	lib-polynomial lib-matrix

matrix:	app/matrix.cc lib-matrix
	$(CC) app/matrix.cc $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o matrix matrix.o

lib-matrix:	$(LIB)/src/matrix.cc $(LIB)/include/matrix.h
	$(CC) $(LIB)/src/matrix.cc $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o matrix.o

lib-polynomial:	$(LIB)/src/polynomial.cc $(LIB)/include/polynomial.h
	$(CC) $(LIB)/src/polynomial.cc $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o polynomial.o

clean:
	rm -rf polynomial polynomial.o matrix matrix.o

.PHONY:	all clean main
