#ifndef POLYNOMIAL_PROJECT_NUMBER_THEORETIC_H
#define POLYNOMIAL_PROJECT_NUMBER_THEORETIC_H

#include <gmpxx.h>

#include <prime_generator.h>

#include <cassert>

namespace gmp_extensions {

    prime_generator &the_prime_generator();

    void invert(mpz_class &res, mpz_class const &a, mpz_class const &p);

    std::pair<mpz_class, mpz_class> extended_euclidean_algorithm(mpz_class const &n, mpz_class const &m);

    mpz_class euler_totient(mpz_class const &);

    [[nodiscard]]
    bool is_prime(mpz_class const &);

    [[nodiscard]]
    std::vector<mpz_class> prime_factorization(mpz_class n);

    [[nodiscard]]
    bool primality_test_fermat_little_thm(mpz_class const &n, size_t iter = 50);
}  // namespace gmp_extensions

#endif