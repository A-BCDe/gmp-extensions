CC = g++
CFLAGS = -Wall -O0
LDFLAGS = -lgmp -lgmpxx
LIB = extensions
INCLUDE = -I$(LIB)/include/

all:	polynomial matrix

polynomial:	app/polynomial.cc lib
	$(CC) app/polynomial.cc $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o polynomial polynomial.o matrix.o prime_generator.o

lib:	lib-polynomial lib-matrix lib-prime-generator

matrix:	app/matrix.cc lib-matrix
	$(CC) app/matrix.cc $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o matrix matrix.o

lib-matrix:	$(LIB)/src/matrix.cc $(LIB)/include/matrix.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o matrix.o

lib-polynomial:	$(LIB)/src/polynomial.cc $(LIB)/include/polynomial.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o polynomial.o

lib-prime-generator:	$(LIB)/src/prime_generator.cc $(LIB)/include/prime_generator.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o prime_generator.o

clean:
	rm -rf polynomial polynomial.o matrix matrix.o

.PHONY:	all clean main
