#include <number_theoretic.h>

namespace project {

	void invert(mpz_class &res, const mpz_class &a, const mpz_class &p) {
#ifndef NDEBUG
		int result =
#endif
		mpz_invert(res.get_mpz_t(), a.get_mpz_t(), p.get_mpz_t());
		assert(result);
	}

	bool is_prime(mpz_class const &n) {
		for(mpz_class i = 2; i * i <= n; i++) {
			if(n % i == 0) return false;
		} return true;
	}

	std::pair<mpz_class, mpz_class> extended_euclidean_algorithm(mpz_class const &n, mpz_class const &m) {
		if(m == 0) return { 1, 0 };
		auto const p = extended_euclidean_algorithm(m, n % m);
		return { p.second, p.first - (n / m) * p.second };
	}

}