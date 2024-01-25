#include <lattice.h>

#include <vector.h>

#include <cassert>
#include <iostream>

namespace gmp_extensions {

    lattice::lattice(std::vector<vector> bases)
     : bases(std::move(bases)) {
        LLL(3, 2);
    }

    void lattice::LLL(mpz_class const &alpha_numerator, mpz_class const &alpha_denominator) {
        LLL(bases, alpha_numerator, alpha_denominator);
    }

    // Thiemann2020_Article_FormalizingTheLLLBasisReductio.pdf
    // Algorithm 5
    void lattice::LLL(std::vector<vector> &f, mpz_class const &num, mpz_class const &denom) {
        assert(is_linearly_independent(f));

        size_t i = 0;
        bool upw = true;
        auto const m = f.size();
        mpz_class c, tmp1, tmp2, tmp3, tmp4;

        // tilde{mu}
        // Algorithm 2
        auto tilde_mu = GSO_computation(f);
        std::vector<mpz_class> d;
        d.reserve(m + 1);
        d.emplace_back(1);
        for(size_t j = 0; j < m; j++) { d.emplace_back(tilde_mu(j, j)); }
        while(i < m) {
            if(upw && i) {
                size_t j = i;
                do {
                    j--;
                    tmp1 = 2 * tilde_mu(i, j) + d[j + 1];
                    tmp2 = 2 * d[j + 1];
                    mpz_fdiv_q(c.get_mpz_t(), tmp1.get_mpz_t(), tmp2.get_mpz_t());
                    if(c == 0) continue;
                    f[i] -= c * f[j];
                    tilde_mu(i, j) -= c * d[j + 1];
                    for(size_t k = 0; k < j; k++) { tilde_mu(i, k) -= c * tilde_mu(j, k); }
                }
                while(j);
            }
            if(!i || d[i] * d[i] * denom <= d[i - 1] * d[i + 1] * num) {
                i++;
                upw = true;
                continue;
            }
            for(size_t j = 0; j < i - 1; j++) { std::swap(tilde_mu(i - 1, j), tilde_mu(i, j)); }
            for(size_t j = i + 1; j < m; j++) {
                tmp1 = tilde_mu(i, i - 1) * tilde_mu(j, i - 1) + tilde_mu(j, i) * d[i - 1];
                tmp2 = d[i + 1] * tilde_mu(j, i - 1) - tilde_mu(i, i - 1) * tilde_mu(j, i);

                mpz_fdiv_q(tmp3.get_mpz_t(), tmp1.get_mpz_t(), d[i].get_mpz_t());
                mpz_fdiv_q(tmp4.get_mpz_t(), tmp2.get_mpz_t(), d[i].get_mpz_t());
                tilde_mu(j, i - 1) = tmp3;
                tilde_mu(j, i) = tmp4;
            }
            tmp1 = d[i + 1] * d[i - 1] + tilde_mu(i, i - 1) * tilde_mu(i, i - 1);
            tmp2 = tmp1 / d[i];
            mpz_fdiv_q(tmp2.get_mpz_t(), tmp1.get_mpz_t(), d[i].get_mpz_t());
            d[i] = tmp2;
            std::swap(f[i], f[i - 1]);
            upw = false;
            i--;
        }
    }

    // Thiemann2020_Article_FormalizingTheLLLBasisReductio.pdf
    // Algorithm 2
    integer_matrix lattice::GSO_computation(std::vector<vector> const &f) {
        auto const m = f.size();
        integer_matrix mu(m, m);
        mpz_class sigma, tmp1;
        for(size_t i = 0; i < m; i++) {
            mu(i, 0) = inner_product(f[i], f[0]);
            for(size_t j = 1; j <= i; j++) {
                sigma = mu(i, 0) * mu(j, 0);
                for(size_t l = 1; l < j; l++) {
                    tmp1 = mu(l, l) * sigma + mu(i, l) * mu(j, l);
                    mpz_fdiv_q(sigma.get_mpz_t(), tmp1.get_mpz_t(), mu(l - 1, l - 1).get_mpz_t());
                }
                mu(i, j) = mu(j - 1, j - 1) * inner_product(f[i], f[j]) - sigma;
            }
        }
        return mu;
    }
}  // namespace gmp_extensions