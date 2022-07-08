#include <integer_random.h>

#include <cassert>

namespace gmp_extensions {

	/**
	 * Picks a random integer from a uniform distribution [lo, hi).
	 * @param lo Lowest possible integer
	 * @param hi Highest possible integer - 1
	 * @return A random integer from a uniform distribution [lo, hi).
	 */
	mpz_class integer_random::operator()(const mpz_class &lo, const mpz_class &hi) {
		assert(hi > lo);
		mpz_class m, n = hi - lo;
		mpz_urandomm(m.get_mpz_t(), randstate, n.get_mpz_t());
		m += lo;
		return m;
	}

}