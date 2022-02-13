#include "polynomial.h"


#include <matrix.h>
#include <prime_generator.h>
#include <polynomial_matrix.h>

#include <cmath>
#include <iostream>

namespace project {

	// polynomial.cc specific functions

	static integer_polynomial make_mono(
			mpz_class const &coef,
			mpz_class deg,
			integer_polynomial const &poly,
			mpz_class const &p) {
		assert(!poly.is_zero());
		assert(coef % p != 0);
		integer_polynomial result({ (coef % p + p) % p });
		integer_polynomial mul = integer_polynomial::monomial(1, 1);
		while(deg > 0) {
			if(mpz_odd_p(deg.get_mpz_t())) {
				result *= mul;
				result.modulo_eq(poly, p);
			}
			mul *= mul;
			mul.modulo_eq(poly, p);
			deg /= 2;
		}
		return result;
	}

	// class member functions

	integer_polynomial::integer_polynomial(std::vector<mpz_class> coef)
	  : coef(std::move(coef)) { clean(); }

	integer_polynomial::integer_polynomial(std::vector<std::string> const &coef) {
		this->coef.reserve(coef.size());
		for(std::string const &now : coef) {
			this->coef.emplace_back(now.c_str());
		}
	}

	integer_polynomial::integer_polynomial(vector const &v) {
		coef.reserve(v.length());
		for(size_t i = 0; i < v.length(); i++) {
			coef.emplace_back(v[i]);
		}
		clean();
	}

	integer_polynomial integer_polynomial::monomial(mpz_class const &coefficient, size_t deg) {
		if(coefficient == 0) return integer_polynomial();
		integer_polynomial polynomial;
		polynomial.coef.reserve(deg + 1);
		for(size_t i = 0; i < deg; i++) {
			polynomial.coef.emplace_back(0);
		} polynomial.coef.emplace_back(coefficient);
		return polynomial;
	}

	integer_polynomial integer_polynomial::constant(mpz_class const &c) {
		if(c == 0) return integer_polynomial();
		return integer_polynomial({ c });
	}



	void integer_polynomial::clean() noexcept {
		size_t i = coef.size() - 1;
		while(true) {
			if(coef[i] != 0) {
				coef.resize(i + 1);
				return;
			}
			if(!i) {
				coef.clear();
				return;
			}
			i--;
		}
	}

	integer_polynomial &integer_polynomial::modulo_eq(integer_polynomial poly, mpz_class const &p) {
		assert(p != 0);
		mpz_class tmp;
		while(poly.leading() % p == 0) poly.coef.pop_back();
		assert(!poly.is_zero());
		if(poly.leading() < 0) poly = -poly;
		if(!poly.is_monic()) {
			int result = mpz_invert(tmp.get_mpz_t(), poly.leading().get_mpz_t(), p.get_mpz_t());
			assert(result);
			for(size_t i = 0; i < poly.coef.size() - 1; i++) {
				poly.coef[i] = (poly.coef[i] * tmp) % p;
			} poly.coef.back() = 1;
		}
		else {
			for(size_t i = 0; i < poly.coef.size() - 1; i++) {
				poly.coef[i] %= p;
			}
		}

		while(coef.size() >= poly.coef.size()) {
			if((coef.back() %= p) == 0) {
				clean();
				continue;
			}
			for(size_t i = 0; i < poly.coef.size(); i++) {
				coef[coef.size() - poly.coef.size() + i] -= poly.coef[i] * leading();
			}
			clean();
		}
		for(auto &now : coef) {
			now = ((now % p) + p) % p;
		}
		return *this;
	}

	integer_polynomial integer_polynomial::modulo(integer_polynomial poly, mpz_class const &p) const {
		integer_polynomial polynomial(*this);
		return polynomial.modulo_eq(std::move(poly), p);
	}

	integer_polynomial &integer_polynomial::derivative_eq() {
		for(size_t i = 0; i < coef.size() - 1; i++) {
			coef[i] = coef[i] * (i + 1);
		}
		clean();
		return *this;
	}

	integer_polynomial integer_polynomial::derivative() const {
		integer_polynomial polynomial(*this);
		return polynomial.derivative_eq();
	}

	// TODO: Fix it.
	// Error when:
	// std::vector<std::string>{ "2", "1", "3", "1", "1" }
	// p = 11
	integer_polynomial &integer_polynomial::divexact_modulo_eq(integer_polynomial const &poly, mpz_class const &p) {
		assert(!poly.is_zero());
		assert(is_prime(p));

		mpz_class inv, tmp;
		int result = mpz_invert(inv.get_mpz_t(), poly.leading().get_mpz_t(), p.get_mpz_t());
		assert(result);

		std::vector<mpz_class> v;
		while(coef.size() >= poly.coef.size()) {
			coef.back() %= p;
			tmp = inv * coef.back();
			for(size_t i = 0; i < poly.coef.size() - 1; i++) {
				coef[coef.size() - poly.coef.size() + i] -= poly.coef[i] * tmp;
			}
			v.emplace_back((tmp % p + p) % p);
			coef.pop_back();
		}
#ifndef NDEBUG
		for(auto & i : coef) assert(i % p == 0);
#endif
		std::reverse(v.begin(), v.end());
		coef = std::move(v);
		clean();

		return *this;
	}

	integer_polynomial integer_polynomial::divexact_modulo(integer_polynomial const &poly, mpz_class const &p) const {
		integer_polynomial polynomial(*this);
		return polynomial.divexact_modulo_eq(poly, p);
	}

	bool integer_polynomial::is_square_free() const {
		return gcd(*this, derivative()).is_constant();
	}

	bool integer_polynomial::is_square_free(mpz_class const &p) const {
		assert(p == 0 || is_prime(p));
		if(p == 0) return is_square_free();
		return gcd(*this, derivative(), p).is_constant();
	}

	integer_polynomial &integer_polynomial::negate() noexcept {
		for(auto &now : coef) now = -now;
		return *this;
	}


	// Based on Thiemann2020_Article_FormalizingTheLLLBasisReductio
	std::vector<integer_polynomial> integer_polynomial::factorize() const {

		assert(is_square_free());
		assert(degree());

		// 1
		mpz_class const two(2);
		mpz_class tmp, f(0);
		mpz_class b = leading();
		mpz_class B;
		mpz_pow_ui(B.get_mpz_t(), two.get_mpz_t(), 5 * degree() * degree());
		for(auto const &now : coef) {
			f += now * now;
		}
		mpz_pow_ui(tmp.get_mpz_t(), f.get_mpz_t(), degree());
		B *= tmp;
		B = sqrt(B) + 1;

		// 2
		prime_generator pg;
		mpz_class p;
		do {
			p = pg.next_prime();
		} while(mpz_divisible_p(leading().get_mpz_t(), p.get_mpz_t()) || is_square_free(p));

		mpz_class l(0);
		tmp = 1;
		while(B >= tmp) {
			tmp *= p;
			++l;
		}

		// 3
		auto factor_p = factorize(p);
		for(auto &now : factor_p) {
			if(now.leading() == 1) continue;
			mpz_invert(tmp.get_mpz_t(), now.leading().get_mpz_t(), p.get_mpz_t());
			for(size_t i = 0; i < now.coef.size() - 1; i++) {
				now.coef[i] = (now.coef[i] * tmp) % p;
			} now.coef.back() = 1;
		}

		// 4


		return {};
	}

	std::vector<integer_polynomial> integer_polynomial::factorize(mpz_class const &p) const {
		// Berlekamp
		integer_matrix L(degree(), degree());
		mpz_class tmp;

		for(size_t i = 0; i < degree(); i++) {
			tmp = p * i;
			auto mono = make_mono(1, tmp, *this, p);
			if(!mono.is_zero()) {
				for (size_t j = 0; j <= mono.degree(); j++) L(j, i) = mono[j];
			}
		}
		for(size_t i = 0; i < degree(); i++) L(i, i) = (L(i, i) + p - 1) % p;

		auto space = kernel(L, p);
		for(auto &v : space) {
			for(size_t i = 0; i < degree(); i++) {
				v[i] = ((v[i] % p) + p) % p;
			}
		}

#ifndef NDEBUG
		for(size_t i = 0; i < degree(); i++) assert(space[0][i] == (!i));
#endif

		std::vector<integer_polynomial> factors;
		factors.reserve(space.size());
		factors.emplace_back(*this);

		// http://newweb.cecm.sfu.ca/CAG/theses/chelsea.pdf
		size_t r = 1;
		while(factors.size() < space.size()) {
			for(auto &poly : factors) {
				assert(r < space.size());
				for(mpz_class i = 0; i < p; i++) {
					auto polynomial = integer_polynomial(space[r]) - constant(i);
					auto g = gcd(polynomial, poly, p);
					assert(!g.is_zero());
					if(!(g.is_constant() && g[0] == 1) && g != poly) {
						poly.divexact_modulo_eq(g, p);
						factors.emplace_back(g);
					}
					if(factors.size() == space.size()) goto END;
				}
			}
			r++;
		}
END:
		return factors;
	}

	mpz_class integer_polynomial::operator()(mpz_class const &x) const {
		mpz_class res(coef[0]), mul(x);
		for(size_t i = 1; i < coef.size(); i++) {
			res += mul * coef[i];
			mul *= x;
		} return res;
	}

	integer_polynomial &integer_polynomial::add_eq(integer_polynomial const &poly) {
		size_t i;
		if(poly.coef.size() > coef.size()) {
			coef.resize(poly.coef.size());
			for(i = 0; i < coef.size(); i++) coef[i] += poly.coef[i];
			for(; i < poly.coef.size(); i++) coef[i] = poly.coef[i];
		}
		else {
			for(i = 0; i < poly.coef.size(); i++) coef[i] += poly.coef[i];
		}
		clean();
		return *this;
	}

	integer_polynomial &integer_polynomial::sub_eq(integer_polynomial const &poly) {
		if(this == &poly) return *this = integer_polynomial();
		size_t i = 0;
		if(poly.coef.size() > coef.size()) {
			coef.resize(poly.coef.size());
			for(i = 0; i < coef.size(); i++) coef[i] -= poly.coef[i];
			for(; i < poly.coef.size(); i++) coef[i] = -poly.coef[i];
		}
		else {
			for(i = 0; i < poly.coef.size(); i++) coef[i] -= poly.coef[i];
		}
		clean();
		return *this;
	}

	integer_polynomial &integer_polynomial::mul_eq(integer_polynomial const &poly) {
		if(is_zero()) return *this;
		std::vector<mpz_class> tmp(coef.size() + poly.coef.size() - 1);
		for(size_t i = 0; i < coef.size(); i++) {
			for(size_t j = 0; j < poly.coef.size(); j++) {
				tmp[i + j] += coef[i] * poly.coef[j];
			}
		}
		coef = std::move(tmp);
		clean();
		return *this;
	}

	integer_polynomial &integer_polynomial::mul_scalar_eq(mpz_class const &n) {
		if(n == 0) return *this = integer_polynomial();
		for(auto &now : coef) now *= n;
		return *this;
	}

	integer_polynomial integer_polynomial::neg() const {
		integer_polynomial polynomial(*this);
		return polynomial.negate();
	}

	integer_polynomial integer_polynomial::add(integer_polynomial const &poly) const {
		integer_polynomial polynomial(*this);
		return polynomial.add_eq(poly);
	}

	integer_polynomial integer_polynomial::sub(integer_polynomial const &poly) const {
		integer_polynomial polynomial(*this);
		return polynomial.sub_eq(poly);
	}

	integer_polynomial integer_polynomial::mul(integer_polynomial const &poly) const {
		integer_polynomial polynomial(*this);
		return polynomial.mul_eq(poly);
	}

	integer_polynomial &integer_polynomial::mul_scalar(mpz_class const &n) const {
		integer_polynomial polynomial(*this);
		return polynomial.mul_scalar_eq(n);
	}

	bool integer_polynomial::equal(integer_polynomial const &poly) const {
		if(coef.size() != poly.coef.size()) return false;
		for(size_t i = 0; i < coef.size(); i++) {
			if(coef[i] != poly.coef[i]) return false;
		} return true;
	}

	std::string integer_polynomial::get_str() const {
		std::string res;
		size_t pow = 0;
		bool start = true;
		if(is_zero()) return res = "0";
		for(auto const &now : coef) {
			if(now == 0) {
				pow++;
				continue;
			}
			if(pow) {
				if(!start && now > 0) res += '+';
				if(now == -1) res += '-';
				else if(now != 1) res += now.get_str();
				if(pow > 1) res += "x^" + std::to_string(pow);
				else res += 'x';
				start = false;
			}
			else { res += now.get_str(); start = false; }
			pow++;
		}
		return res;
	}


	// friend functions

	[[nodiscard]] integer_polynomial gcd(integer_polynomial a, integer_polynomial b) {
		if(a.coef.size() < b.coef.size()) return gcd(std::move(b), std::move(a));
		if(b.coef.empty()) {
			mpz_class g = a.coef[0];
			mpz_class q;
			for(size_t i = 1; i < a.coef.size(); i++) {
				g = gcd(g, a.coef[i]);
			} for(auto &now : a.coef) {
				mpz_divexact(q.get_mpz_t(), now.get_mpz_t(), g.get_mpz_t());
				now = q;
			} return a;
		}
		mpz_class leading_gcd = gcd(a.leading(), b.leading());
		mpz_class q;
		mpz_divexact(q.get_mpz_t(), b.leading().get_mpz_t(), leading_gcd.get_mpz_t());
		a *= q;
		mpz_divexact(q.get_mpz_t(), a.leading().get_mpz_t(), leading_gcd.get_mpz_t());
		b *= q;
		for(size_t i = 0; i < b.coef.size(); i++) {
			a.coef[a.coef.size() - b.coef.size() + i] -= b.coef[i];
		} a.clean();
		return gcd(std::move(b), std::move(a));
	}

	integer_polynomial gcd(integer_polynomial a, integer_polynomial b, const mpz_class &p) {
		assert(is_prime(p));
		if(a.is_zero()) return gcd(std::move(b), std::move(a), p);
		if(b.is_zero()) {
			// convert all coefficients of `a` to non-negative
			// and modulo p
			for(auto &now : a.coef) {
				now = ((now % p) + p) % p;
			}

			// multiply inverse of leading coefficient
			// and modulo p
			mpz_class inv;
			int result = mpz_invert(inv.get_mpz_t(), a.leading().get_mpz_t(), p.get_mpz_t());
			assert(result);
			for(size_t i = 0; i < a.coef.size() - 1; i++) {
				a.coef[i] = (a.coef[i] * inv) % p;
			} a.coef.back() = 1;
			return a;
		}
		if(a.leading() < 0) return gcd(-a, std::move(b), p);
		if(b.leading() < 0) return gcd(std::move(a), -b, p);
		if(a.coef.size() < b.coef.size()) return gcd(std::move(b), std::move(a), p);

		mpz_class inv;
		int result = mpz_invert(inv.get_mpz_t(), b.leading().get_mpz_t(), p.get_mpz_t());
		assert(result);
		for(size_t i = 0; i < b.coef.size() - 1; i++) {
			b.coef[i] = (b.coef[i] * inv * a.leading()) % p;
		} b.coef.back() = a.leading();

		for(size_t i = 0; i < b.coef.size(); i++) {
			auto const off = a.coef.size() - b.coef.size() + i;
			a.coef[off] -= b.coef[i];
			a.coef[off] %= p;
		} a.clean();
		return gcd(std::move(b), std::move(a), p);
	}

	std::ostream &operator<<(std::ostream &os, integer_polynomial const &poly) {
		os << poly.get_str();
		return os;
	}

	// outside class

	bool is_prime(mpz_class const &n) {
		for(mpz_class i = 2; i * i <= n; i++) {
			if(n % i == 0) return false;
		} return true;
	}
}
