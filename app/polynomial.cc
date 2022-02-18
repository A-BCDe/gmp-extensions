#include <gmpxx.h>

#include <iostream>
#include <prime_generator.h>
#include <polynomial.h>

int main() {
	std::cout << "Hello, world!" << std::endl;
	//project::integer_polynomial poly(std::vector<std::string>{ "1", "-3", "-1", "-3", "1", "-3", "1" });
	project::integer_polynomial poly(std::vector<std::string>{ "12", "41", "54", "41", "23", "8", "1" });

	mpz_class p = 11;
	std::cout << poly << '\n';
	auto factors = poly.factorize(p);
	std::cout << "Irreducible factors of " << poly << " mod " << p << " are:\n";
	for(auto const &v : factors) {
		std::cout << v << '\n';
	}

	auto pair = poly.hensel_lifting(factors[0], factors[1] * factors[2], p, 0, 11);
	auto g = std::move(pair.first);
	auto h = std::move(pair.second);
	std::cout << "g = " << g << '\n';
	std::cout << "h = " << h << '\n';
}
