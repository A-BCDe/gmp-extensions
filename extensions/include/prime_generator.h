#ifndef POLYNOMIAL_PROJECT_PRIME_GENERATOR_H
#define POLYNOMIAL_PROJECT_PRIME_GENERATOR_H

#include <gmpxx.h>

#include <functional>
#include <mutex>
#include <ostream>
#include <vector>

namespace gmp_extensions {

	class prime_generator {
	public:
		explicit prime_generator(size_t sieve_size = 10000000, size_t bound = 1);

		void print_primes(std::ostream &os) const;
		mpz_class next_prime();

		std::vector<mpz_class> const &found_primes() const { return primes; }

		mpz_class smallest_prime(std::function<bool(mpz_class const&)> const&);

	private:
		size_t const sieve_size;
		mpz_class bound;
		std::vector<mpz_class> primes; // found primes up to now
		std::vector<bool> sieve;
		mutable std::mutex mutex;

		mpz_class next_prime_inner();
	};

}

#endif
