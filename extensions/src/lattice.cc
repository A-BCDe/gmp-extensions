#include <lattice.h>

#include <vector.h>

#include <cassert>
#include <iostream>

namespace project {

	lattice::lattice(std::vector<vector> bases) : bases(std::move(bases)) {
		LLL(3, 2);
		std::cout << "bases:\n";
		for(auto const &now : this->bases) {
			std::cout << now << '\n';
		}
	}

	void lattice::LLL(mpz_class const &alpha_numerator, mpz_class const &alpha_denominator) {
		LLL(bases, alpha_numerator, alpha_denominator);
	}

	// Thiemann2020_Article_FormalizingTheLLLBasisReductio.pdf
	// Algorithm 5
	void lattice::LLL(std::vector<vector> &f,
	                  mpz_class const &alpha_numerator,
	                  mpz_class const &alpha_denominator) {
		assert(is_linearly_independent(f));
		size_t i = 0;
		bool upw = true;
		auto const m = f.size();
		mpz_class tmp1, tmp2;

		// tilde{mu}
		// Algorithm 2
		auto tilde_mu = GSO_computation(f);
		std::cout << "tilde_mu:\n" << tilde_mu << "\n\n";
		std::vector<mpz_class> d;
		d.reserve(m + 1);
		d.emplace_back(1);
		for(size_t j = 0; j < f.size(); j++) {
			//std::cout << "!j = " << j << std::endl;
			d.emplace_back(tilde_mu(j, j));
		}
		while(i < m) {
			if(upw && i) {
				size_t j = i;
				do {
					j--;
					tmp1 = 2 * tilde_mu(i, j) + d[j + 1];
					tmp2 = 2 * d[j + 1];
					mpz_class c;
					//c = tmp1 / tmp2;
					mpz_fdiv_q(c.get_mpz_t(), tmp1.get_mpz_t(), tmp2.get_mpz_t());
					if(c == 0) continue;
					f[i] -= c * f[j];
					tilde_mu(i, j) -= c * d[j + 1];
					for(size_t k = 0; k < j; k++) {
						tilde_mu(i, k) -= c * tilde_mu(j, k);
					}
				} while(j);
			}
			if(!i || d[i] * d[i] * alpha_denominator <= d[i - 1] * d[i + 1] * alpha_numerator) {
				i++;
				upw = true;
				continue;
			}
			for(size_t j = 0; j < i - 1; j++) {
				std::swap(tilde_mu(i - 1, j), tilde_mu(i, j));
			}
			for(size_t j = i + 1; j < m; j++) {
				tmp1 = tilde_mu(i, i - 1) * tilde_mu(j, i - 1) + tilde_mu(j, i) * d[i - 1];
				tmp2 = d[i + 1] * tilde_mu(j, i - 1) - tilde_mu(i, i - 1) * tilde_mu(j, i);

				//tilde_mu(j, i - 1) = tmp1 / d[i];
				//tilde_mu(j, i) = tmp2 / d[i];
				mpz_fdiv_q(tilde_mu(j, i - 1).get_mpz_t(), tmp1.get_mpz_t(), d[i].get_mpz_t());
				mpz_fdiv_q(tilde_mu(j, i).get_mpz_t(), tmp2.get_mpz_t(), d[i].get_mpz_t());
			}
			tmp1 = d[i + 1] * d[i - 1] + tilde_mu(i, i - 1) * tilde_mu(i, i - 1);
			//tmp2 = tmp1 / d[i];
			mpz_fdiv_q(tmp2.get_mpz_t(), tmp1.get_mpz_t(), d[i].get_mpz_t());
			d[i] = tmp2;
			std::swap(f[i], f[i - 1]);
			upw = false;
			i--;
		}
	}

	// Thiemann2020_Article_FormalizingTheLLLBasisReductio.pdf
	// Algorithm 2
	integer_matrix lattice::GSO_computation(std::vector<vector> const &bases) {
		auto const m = bases.size();
		integer_matrix mu(m, m);
		mpz_class sigma, tmp1;
		for (size_t i = 0; i < m; i++) {
			mu(i, 0) = inner_product(bases[i], bases[0]);
			for (size_t j = 1; j <= i; j++) {
				sigma = mu(i, 0) * mu(j, 0);
				for (size_t k = 1; k < j; k++) {
					tmp1 = mu(k, k) * sigma + mu(i, k) * mu(j, k);
					//sigma = tmp1 / mu(k - 1, k - 1);
					mpz_fdiv_q(sigma.get_mpz_t(), tmp1.get_mpz_t(), mu(k - 1, k - 1).get_mpz_t());
				}
				mu(i, j) = mu(j - 1, j - 1) * inner_product(bases[i], bases[j]) - sigma;
			}
		}
		return mu;
	}
}