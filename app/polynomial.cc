#include <gmpxx.h>

#include <iostream>
#include <prime_generator.h>
#include <polynomial.h>

int main() {
	std::cout << "Hello, world!" << std::endl;
	project::integer_polynomial poly({"-2", "0", "1", "1"});
	std::cout << poly << '\n';
	poly.factorize(3);
}
