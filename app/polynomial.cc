#include <gmpxx.h>

#include <iostream>
#include <prime_generator.h>
#include <polynomial.h>

int main() {
	std::cout << "Hello, world!" << std::endl;
	project::integer_polynomial poly(std::vector<std::string>{ "1", "-3", "-1", "-3", "1", "-3", "1" });

	mpz_class p = 11;
	std::cout << poly << '\n';
	auto factors = poly.factorize(p);
	std::cout << "Irreducible factors of " << poly << " mod " << p << " are:\n";
	for(auto const &v : factors) {
		std::cout << v << '\n';
	}
}
