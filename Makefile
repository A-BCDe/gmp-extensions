CC = g++
CFLAGS = --std=c++14 -Wall -pthread -lpthread -O2 -DNDEBUG
LDFLAGS = -lgmp -lgmpxx
LIB = extensions
INCLUDE = -I$(LIB)/include/

all:	polynomial matrix

polynomial:	examples/polynomial.cc lib
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o polynomial polynomial.o matrix.o prime_generator.o vector.o number_theoretic.o lattice.o integer_random.o

lib:	lib-integer-vector lib-polynomial lib-integer-matrix lib-prime-generator lib-polynomial-matrix lib-number-theoretic lib-integer-lattice lib-integer-random

matrix:	examples/matrix.cc lib-integer-matrix lib-number-theoretic lib-prime-generator lib-integer-random
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o matrix matrix.o vector.o number_theoretic.o prime_generator.o integer_random.o

lib-number-theoretic:	$(LIB)/src/number_theoretic.cc $(LIB)/include/number_theoretic.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o number_theoretic.o

lib-integer-vector:	$(LIB)/src/vector.cc $(LIB)/include/vector.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o vector.o

lib-integer-matrix:	$(LIB)/src/matrix.cc $(LIB)/include/matrix.h lib-integer-vector lib-number-theoretic
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o matrix.o

lib-integer-lattice:	$(LIB)/src/lattice.cc $(LIB)/include/lattice.h lib-integer-vector lib-integer-matrix lib-number-theoretic
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o lattice.o

lib-polynomial:	$(LIB)/src/polynomial.cc $(LIB)/include/polynomial.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o polynomial.o

lib-polynomial-matrix:	$(LIB)/src/polynomial_matrix.cc $(LIB)/include/polynomial_matrix.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o polynomial_matrix.o

lib-prime-generator:	$(LIB)/src/prime_generator.cc $(LIB)/include/prime_generator.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o prime_generator.o

lib-integer-random:		$(LIB)/src/integer_random.cc $(LIB)/include/integer_random.h
	$(CC) $< $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o integer_random.o

clean:
	rm -rf polynomial polynomial.o matrix matrix.o

.PHONY:	all clean main
