#include "number_theoretic.h"

#include <mutex>
#include <random>

#include <prime_generator.h>
#include <integer_random.h>

namespace project {

	/**
	 * A single prime generator for less computation.
	 * @return The prime generator
	 */
	prime_generator &the_prime_generator() {
		static std::mutex mutex;
		mutex.lock();
		static prime_generator pg;
		mutex.unlock();
		return pg;
	}

	/**
	 * Inverse of a modulo p.
	 * @param res inverse of a modulo p
	 * @param a An integer
	 * @param p A non-zero integer
	 */
	void invert(mpz_class &res, const mpz_class &a, const mpz_class &p) {
#ifndef NDEBUG
		int result =
#endif
		mpz_invert(res.get_mpz_t(), a.get_mpz_t(), p.get_mpz_t());
		assert(result);
	}

	/**
	 * Calculates a and b which n * a + m * b = gcd(n, m).
	 * @param n An integer
	 * @param m An integer
	 * @return { a, b } which n * a + m * b = gcd(n, m).
	 */
	std::pair<mpz_class, mpz_class> extended_euclidean_algorithm(mpz_class const &n, mpz_class const &m) {
		if(m == 0) return { 1, 0 };
		auto const p = extended_euclidean_algorithm(m, n % m);
		return { p.second, p.first - (n / m) * p.second };
	}

	/**
	 * Calculates an euler totient of n.
	 * @param n A positive integer
	 * @return An euler totient of n
	 */
	mpz_class euler_totient(mpz_class n) {
		auto &pg = the_prime_generator();
		mpz_class res(n), tmp;
		bool flag = true;
		for(auto &p : pg.found_primes()) {
			if(p * p > n) {
				flag = false;
				break;
			}
			if(n % p == 0) {
				mpz_remove(tmp.get_mpz_t(), n.get_mpz_t(), p.get_mpz_t());
				n = tmp;
				res /= p;
				res *= p - 1;
			}
		}
		if(n == 1) {
			return res;
		}
		if(!flag) {
			while(true) {
				mpz_class const p = pg.next_prime();
				if(p * p > n) {
					break;
				}
				if(n % p == 0) {
					mpz_remove(tmp.get_mpz_t(), n.get_mpz_t(), p.get_mpz_t());
					n = tmp;
					res /= p;
					res *= p - 1;
				}
			}
		}
		if(n != 1) {
			res /= n;
			res *= n - 1;
		}
		return res;
	}

	/**
	 * Checks if n is prime.
	 * @param n An integer
	 * @return true if n is prime, false if not.
	 */
	bool is_prime(mpz_class const &n) {
		for(mpz_class i = 2; i * i <= n; i++) {
			if(n % i == 0) return false;
		} return true;
	}

	/**
	 * Prime factorization of n.
	 * @param n An integer
	 * @return A vector of primes.
	 */
	std::vector<mpz_class> prime_factorization(mpz_class n) {
		std::vector<mpz_class> fac;
		for(mpz_class i = 2; i * i <= n; ++i) {
			while(n % i == 0) {
				n /= i;
				fac.push_back(i);
			}
		}
		if(n != 1) {
			fac.push_back(n);
		}
		return fac;
	}

	/**
	 * Check if n is prime using Fermat's little theorem.
	 * There is a possibility that n is not prime even if the function returns true.
	 * @param n A positive integer to check primality.
	 * @param iter Maximum number of iterations.
	 * @return The result of the primality test.
	 */
	bool primality_test_fermat_little_thm(mpz_class const &n, mpz_class const &iter) {
		assert(n > 0);
		mpz_class tmp;
		if(iter >= n) {
			for(mpz_class i = 2; i < n; i++) {
				mpz_powm(tmp.get_mpz_t(), i.get_mpz_t(), n.get_mpz_t(), n.get_mpz_t());
				if(tmp != i) {
					return false;
				}
			}
			return true;
		}

		integer_random ir;
		for(mpz_class i = 0; i < iter; i++) {
			mpz_class const rnd = ir(2, n);
			mpz_powm(tmp.get_mpz_t(), rnd.get_mpz_t(), n.get_mpz_t(), n.get_mpz_t());
			if(tmp != i) {
				return false;
			}
		}
		return true;
	}
}