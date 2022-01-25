CC = g++
CFLAGS = -Wall
LDFLAGS = -lgmp -lgmpxx
INCLUDE = -Iinclude/

all:	main

main:	src/main.cc integer-polynomial
	$(CC) src/main.cc $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o main

integer-polynomial:	src/integer_polynomial.cc include/integer_polynomial.h
	$(CC) src/integer_polynomial.cc $(INCLUDE) $(CFLAGS) $(LDFLAGS) -c -o integer-polynomial.o

clean:
	rm -rf main integer-polynomial.o

.PHONY:	all clean