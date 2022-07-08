#include "prime_generator.h"

#include <algorithm>
#include <iostream>

namespace gmp_extensions {
	// assume unique call
	/**
	 * Constructor.
	 * Constructs a prime_generator.
	 * @param sieve_size Sieve size
	 * @param bound The last value the prime_generator has checked
	 */
	prime_generator::prime_generator(size_t sieve_size, size_t bound)
	  : sieve_size(sieve_size), bound(bound), sieve(sieve_size) {
		for(size_t pivot = 0; pivot < bound; pivot += sieve_size) {
			std::fill(sieve.begin(), sieve.end(), false);

			// fill sieve with precalculated primes
			for(auto const &p : primes) {
				auto const pui = p.get_si();
				mpz_class idx = (pivot / p) * p - pivot;
				if(idx < 0) idx += p;
				if(idx >= sieve_size) continue;
				auto const idxui = idx.get_ui();
				for(auto i = idxui; i < sieve_size; i += pui) {
					sieve[i] = true;
				}
			}

			// manual sieve filling
			for(size_t i = 0; i < sieve_size && i < bound - pivot; i++) {
				auto const cur = i + pivot;
				if(cur < 2) continue;
				if(sieve[i]) continue;
				primes.emplace_back(cur);
				for(size_t j = i + cur; j < sieve_size; j += cur) {
					sieve[j] = true;
				}
			}
		}
	}

	/**
	 * Prints primes calculated up to now to os.
	 * @param os ostream to print primes.
	 */
	void prime_generator::print_primes(std::ostream &os) const {
		std::lock_guard<std::mutex> guard(mutex);
		for(auto const &p : primes) {
			os << p.get_str() << ' ';
		}
	}

	/**
	 * Calculate next prime the prime_generator has not calculated.
	 * @return The next prime the prime_generator has not calculated
	 */
	mpz_class prime_generator::next_prime() {
		std::lock_guard<std::mutex> guard(mutex);
		return next_prime_inner();
	}

	/**
	 * Calculate the smallest prime p which makes f(p) true
	 * @param f Some function which gets a prime as an input and returns a boolean value.
	 * @return the smallest prime p which makes f(p) true
	 */
	mpz_class prime_generator::smallest_prime(std::function<bool(mpz_class const&)> const &f) {
		std::lock_guard<std::mutex> guard(mutex);
		for(auto const &p : primes) {
			if(f(p)) return p;
		}
		while(true) {
			auto p = next_prime_inner();
			if(f(p)) return p;
		}
	}

	/**
	 * Helper function to manage multithreading.
	 * @return A prime.
	 */
	mpz_class prime_generator::next_prime_inner() {
		if(bound < 2) bound = 2;
		while(true) {
			mpz_class pivot = (bound / sieve_size) * sieve_size;
			mpz_class index = bound - pivot;
			auto const indexui = index.get_ui();
			if(index == 0) {
				std::fill(sieve.begin(), sieve.end(), false);

				// fill sieve with precalculated primes
				for(auto const &p : primes) {
					auto const pui = p.get_si();
					mpz_class idx = (pivot / p) * p - pivot;
					if(idx < 0) idx += p;
					if(idx >= sieve_size) continue;
					auto const idxui = idx.get_ui();
					for(auto i = idxui; i < sieve_size; i += pui) {
						sieve[i] = true;
					}
				}
			}

			if(!sieve[indexui]) {
				primes.emplace_back(bound);
				if(bound < sieve_size) {
					auto const boundui = bound.get_ui();
					for(auto i = indexui + boundui; i < sieve_size; i += boundui) {
						sieve[i] = true;
					}
				}
				return bound++;
			}

			++bound;
		}
	}
}