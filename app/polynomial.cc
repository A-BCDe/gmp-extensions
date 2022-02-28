#include <gmpxx.h>

#include <iostream>
#include <prime_generator.h>
#include <polynomial.h>

int main() {
	std::cout << "Hello, world!" << std::endl;
	//project::integer_polynomial poly(std::vector<std::string>{ "1", "-3", "-1", "-3", "1", "-3", "1" });
	// (x + 1)(x^2 + 5x + 3)(x^3 + 2x^2 + 3x + 4)
	//project::integer_polynomial poly(std::vector<std::string>{ "12", "41", "54", "41", "23", "8", "1" });
	// (x - 1)(x - 2)
	//project::integer_polynomial poly(std::vector<std::string>{ "2", "-3", "1" });
	// (x - 1)(x - 2)(x + 3)
	//project::integer_polynomial poly(std::vector<std::string>{ "6", "-7", "0", "1" });
	// (x^2 + 1)(x^2 + 2)
	//project::integer_polynomial poly(std::vector<std::string>{ "2", "0", "3", "0", "1" });
	// x^2 + 2
	//project::integer_polynomial poly(std::vector<std::string>{ "2", "0", "1" });
	// (x^2 + 2)(x^3 + 2)
	//project::integer_polynomial poly(std::vector<std::string>{ "4", "0", "2", "2", "0", "1" });
	// (x^2 + 2)(x^3 + 3)
	//project::integer_polynomial poly(std::vector<std::string>{ "6", "0", "3", "2", "0", "1" });
	// (x^2 + 2)(x^2 - 2)
	project::integer_polynomial poly(std::vector<std::string>{ "-4", "0", "0", "0", "1" });

	mpz_class const p = 3;
	std::cout << poly << '\n';
	auto factors = poly.factorize(p);
	std::cout << "Irreducible factors of " << poly << " mod " << p << " are:\n";
	for(auto const &v : factors) {
		std::cout << v << '\n';
	}

	factors = poly.factorize();
	std::cout << "Irreducible factors of " << poly << " are:\n";
	for(auto const &v : factors) {
		std::cout << v << '\n';
	}
}
