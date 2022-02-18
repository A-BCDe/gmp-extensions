#ifndef POLYNOMIAL_PROJECT_NUMBER_THEORETIC_H
#define POLYNOMIAL_PROJECT_NUMBER_THEORETIC_H

#include <gmpxx.h>

#include <prime_generator.h>

#include <cassert>

namespace project {

	void invert(mpz_class &res, mpz_class const &a, mpz_class const &p);

	std::pair<mpz_class, mpz_class> extended_euclidean_algorithm(mpz_class const &n, mpz_class const &m);

	[[nodiscard]] bool is_prime(mpz_class const&);

}

#endif