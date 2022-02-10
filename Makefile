CC = g++
CFLAGS = -Wall -O0
LDFLAGS = -lgmp -lgmpxx
LIB = extensions
INCLUDE = -I$(LIB)/include/

all:	polynomial matrix

polynomial:	app/polynomial.cc lib
	$(CC) app/polynomial.cc $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o polynomial polynomial.o matrix.o prime_generator.o vector.o

lib:	lib-integer-vector lib-polynomial lib-integer-matrix lib-prime-generator lib-polynomial-matrix

matrix:	app/matrix.cc lib-integer-matrix
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o matrix matrix.o

lib-integer-vector:	$(LIB)/src/vector.cc $(LIB)/include/vector.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o vector.o

lib-integer-matrix:	$(LIB)/src/matrix.cc $(LIB)/include/matrix.h lib-integer-vector
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o matrix.o

lib-polynomial:	$(LIB)/src/polynomial.cc $(LIB)/include/polynomial.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o polynomial.o

lib-polynomial-matrix:	$(LIB)/src/polynomial_matrix.cc $(LIB)/include/polynomial_matrix.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o polynomial_matrix.o

lib-prime-generator:	$(LIB)/src/prime_generator.cc $(LIB)/include/prime_generator.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o prime_generator.o

clean:
	rm -rf polynomial polynomial.o matrix matrix.o

.PHONY:	all clean main
