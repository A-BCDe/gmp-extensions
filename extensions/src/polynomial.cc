#include "polynomial.h"

#include <lattice.h>

#include <matrix.h>
#include <prime_generator.h>
#include <number_theoretic.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>

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

	inline static size_t find_power(size_t l, size_t r, mpz_class const &p, mpz_class const &n) {
		mpz_class tmp;
		while(l < r) {
			auto const m = ((l + r) >> 1);
			mpz_pow_ui(tmp.get_mpz_t(), p.get_mpz_t(), m);
			if(n >= tmp) l = m + 1;
			else r = m;
		}
		return l;
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
		if(is_zero()) return;
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

	bool integer_polynomial::is_divisible_modulo(integer_polynomial const &poly, mpz_class const &p) const {
		assert(!poly.is_zero());
		assert(is_prime(p));

		integer_polynomial polynomial(*this);
		mpz_class inv, tmp;
		invert(inv, poly.leading(), p);

		integer_polynomial monic_poly(inv * poly);
		for(auto &now : monic_poly.coef) {
			now %= p;
		}

		while(polynomial.coef.size() >= poly.coef.size()) {
			polynomial.coef.back() %= p;
			if(polynomial.leading() != 0) {
				for(size_t i = 1; i < poly.coef.size(); i++) {
					polynomial.coef[polynomial.degree() - i]
						-= polynomial.leading() * monic_poly[monic_poly.degree() - i];
				}
				polynomial.coef.back() = 0;
			}
			polynomial.clean();
		}

		return std::all_of(polynomial.coef.begin(), polynomial.coef.end(), [&p](mpz_class const &n) {
			return n % p == 0;
		});
	}

	bool integer_polynomial::is_divisible(integer_polynomial const &poly) const {
		assert(!poly.is_zero());
		auto *polynomial_p = divexact(poly);
		if(!polynomial_p) {
			return false;
		}
		delete[] polynomial_p;
		return true;
	}

	integer_polynomial &integer_polynomial::modulo_eq(integer_polynomial poly, mpz_class const &p) {
		assert(p != 0);
		mpz_class tmp;
		while(poly.leading() % p == 0) poly.coef.pop_back();
		assert(!poly.is_zero());
		if(poly.leading() < 0) poly = -poly;
		if(!poly.is_monic()) {
			invert(tmp, poly.leading(), p);
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

	integer_polynomial &integer_polynomial::primitive_eq() {
		if(is_zero()) return *this;
		mpz_class g = coef[0];
		mpz_class tmp;
		for(size_t i = 1; i < coef.size(); i++) {
			g = gcd(g, coef[i]);
		}
		for(auto &now : coef) {
			mpz_divexact(tmp.get_mpz_t(), now.get_mpz_t(), g.get_mpz_t());
			now = tmp;
		}
		return *this;
	}

	integer_polynomial integer_polynomial::primitive() const {
		integer_polynomial polynomial(*this);
		return polynomial.primitive_eq();
	}

	integer_polynomial &integer_polynomial::derivative_eq() {
		for(size_t i = 0; i < coef.size() - 1; i++) {
			coef[i] = coef[i + 1] * (i + 1);
		} coef.pop_back();
		clean();
		return *this;
	}

	integer_polynomial integer_polynomial::derivative() const {
		integer_polynomial polynomial(*this);
		return polynomial.derivative_eq();
	}

	integer_polynomial *integer_polynomial::divexact(integer_polynomial const &poly) const {
		assert(!poly.is_zero());
		if(is_zero()) {
			return new(std::nothrow) integer_polynomial;
		}
		if(coef.size() < poly.coef.size()) {
			return nullptr;
		}
		integer_polynomial polynomial(*this);
		integer_polynomial quotient;
		mpz_class tmp;

		size_t deg = polynomial.coef.size() - poly.coef.size() + 1;

		do {
			deg--;
			if(!mpz_divisible_p(polynomial.coef.back().get_mpz_t(), poly.leading().get_mpz_t())) {
				return nullptr;
			}
			mpz_divexact(tmp.get_mpz_t(), polynomial.coef.back().get_mpz_t(), poly.leading().get_mpz_t());
			for(size_t i = 0; i < poly.coef.size() - 1; i++) {
				polynomial.coef[deg + i] -= tmp * poly[i];
			} polynomial.coef.pop_back();
			quotient.coef.emplace_back(tmp);
		} while(deg);

		for(auto const &now : polynomial.coef) {
			if(now != 0) {
				return nullptr;
			}
		}

		std::reverse(quotient.coef.begin(), quotient.coef.end());

		return new(std::nothrow) integer_polynomial(std::move(quotient));
	}

	integer_polynomial &integer_polynomial::divexact_modulo_eq(integer_polynomial const &poly, mpz_class const &p) {
		assert(!poly.is_zero());
		assert(is_prime(p));

		mpz_class inv, tmp;
		invert(inv, poly.leading(), p);

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
		mpz_class tmp, f;
		mpz_class b = leading();
		mpz_class B, Bp;
		mpz_ui_pow_ui(B.get_mpz_t(), 2, 5 * degree() * degree());

		for(auto const &now : coef) {
			f += now * now;
		}
		mpz_pow_ui(tmp.get_mpz_t(), f.get_mpz_t(), degree());
		B *= tmp;
		Bp = sqrt(B);
		if(Bp * Bp != B) ++Bp;
		B = Bp;

		// 2
		prime_generator pg;
		mpz_class p, p_power_l;
		do {
			p = pg.next_prime();
		} while(mpz_divisible_p(leading().get_mpz_t(), p.get_mpz_t()) || !is_square_free(p));

		size_t l = 1;
		tmp = p;
		while(B > tmp) {
			tmp *= tmp;
			l <<= 1;
		}
		l = find_power((l >> 1) + 1, l + 1, p, B);
		mpz_pow_ui(p_power_l.get_mpz_t(), p.get_mpz_t(), l);

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
		integer_polynomial tmp_poly1;
		integer_polynomial tmp_poly2;
		integer_polynomial tmp_poly3(*this);
		std::vector<std::pair<integer_polynomial, int>> hensel_v;

		size_t idx = 0;
		for(auto const &poly : factor_p) {
			tmp_poly1 = poly;
			tmp_poly2 = tmp_poly3.divexact_modulo(tmp_poly1, p);
			auto pair = tmp_poly3.hensel_lifting(tmp_poly1, tmp_poly2, p, 0, l);
			hensel_v.emplace_back(std::move(pair.first), idx++);
			tmp_poly3 = std::move(pair.second);
		}

		// 5 ~ 14
		std::set<size_t> T; // T
		std::vector<integer_polynomial> result; // G
		for(size_t i = 0; i < hensel_v.size(); i++) T.insert(i);
		std::sort(hensel_v.begin(), hensel_v.end(),
				  [](std::pair<integer_polynomial, size_t> const &a,
						  std::pair<integer_polynomial, size_t> const &b) {
			return a.first.degree() < b.first.degree();
		});
		integer_polynomial f_star(*this);
		size_t offset = 0;
	LP: while(!T.empty()) {
			assert(offset < T.size());
			auto const current_idx = *std::next(T.rbegin(), static_cast<ptrdiff_t>(offset++));
			auto u = hensel_v[current_idx].first;

			std::set<size_t> degree_set; // Deg
			degree_set.insert(0);
			for(auto now_idx : T) { // Dynamic Programming
				if(now_idx == current_idx) continue;
				size_t now_deg = hensel_v[now_idx].first.degree();
				auto degree_set_tmp = degree_set;
				for(auto prev : degree_set) {
					degree_set_tmp.insert(prev + now_deg);
				}
				degree_set = std::move(degree_set_tmp);
			}

			for(size_t now_deg : degree_set) {
				size_t const k = now_deg + 1;
				size_t const j = u.coef.size() + now_deg;
				mpz_ui_pow_ui(B.get_mpz_t(), 2, 5 * j * j);

				f = 0; // ||f*||
				for(auto const &now : f_star.coef) f += now * now;
				mpz_pow_ui(tmp.get_mpz_t(), f.get_mpz_t(), j << 1);
				Bp = sqrt(B * tmp);
				if(Bp * Bp < B) ++Bp;
				B = Bp;

				size_t lp = 1; // l'
				tmp = p;
				while(B > tmp) {
					tmp *= tmp;
					lp <<= 1;
				}
				lp = find_power((lp >> 1) + 1, lp + 1, p, B);
				if(l > lp) {
					l = lp;
					mpz_pow_ui(p_power_l.get_mpz_t(), p.get_mpz_t(), l);
				}

				// 12: LLL
				// Make a lattice L_{u,k}
				std::vector<vector> vectors;
				vectors.reserve(j);
				for(size_t i = 0; i < k; i++) {
					vector v(j);
					for(size_t ii = 0; ii < u.coef.size(); ii++) {
						v[ii + i] = u[u.degree() - ii];
					}
					vectors.emplace_back(std::move(v));
				}
				for(size_t i = k; i < j; i++) {
					vectors.emplace_back(j);
					vectors.back()[i] = p_power_l;
				}

				lattice lat(vectors);
				auto const g = integer_polynomial(lat.short_vector().reverse());
				integer_polynomial ppg(g.primitive());

				// 14
				if(abs(g.leading()) < p_power_l) {
					auto *poly_p = divexact(ppg);
					if(!poly_p) {
						continue;
					}

					// 13, calculate T <- T - S
					for(auto it = T.begin(); it != T.end(); ) {
						if(ppg.is_divisible_modulo(factor_p[hensel_v[*it].second], p)) {
							it = T.erase(it);
							continue;
						}
						++it;
					}

					f_star = std::move(*poly_p);
					poly_p = nullptr;
					delete poly_p;
					b = f_star.leading();

					result.emplace_back(std::move(ppg));
					offset = 0;
					goto LP;
				}
			}
		}

		// Removed 15th line as it seems incorrect
		return result;
	}

	std::vector<integer_polynomial> integer_polynomial::factorize(mpz_class const &p) const {
		assert(is_prime(p));

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

	// Mignotte1992_Book_MathematicsForComputerAlgebra.pdf
	std::pair<integer_polynomial, integer_polynomial> integer_polynomial::hensel_lifting(
			integer_polynomial g,
			integer_polynomial h,
			mpz_class const &p,
			size_t k,
			size_t desired_power) const {
		mpz_class tmp;
		mpz_class pk2i; // p^(k + 2^i)
		mpz_class current_modulo;
		mpz_pow_ui(tmp.get_mpz_t(), p.get_mpz_t(), k); // tmp = p^k
		current_modulo = tmp * tmp * p; // current_modulo = p^(2k + 1)
		pk2i = tmp * p; // pk2i = p^(k+1)
		integer_polynomial diff(*this - g * h);
#ifndef NDEBUG
		mpz_class tmp2 = resultant(g, h);
		assert(mpz_divisible_p(tmp2.get_mpz_t(), tmp.get_mpz_t()));
		tmp *= p;
		assert(!mpz_divisible_p(tmp2.get_mpz_t(), tmp.get_mpz_t()));
		mpz_divexact(tmp2.get_mpz_t(), tmp.get_mpz_t(), p.get_mpz_t());

		for(size_t i = 0; i < diff.coef.size(); i++) {
			assert(mpz_divisible_p(diff[i].get_mpz_t(), current_modulo.get_mpz_t()));
		}
#endif
		mpz_pow_ui(tmp.get_mpz_t(), p.get_mpz_t(), 1);
		// tmp = p
		// current_modulo = p^(2k + 1)
		// pk2i = p^(k + 1)
		size_t current_power = ((k << 1) | 1); // current_power = 2k + 1
		size_t current_2i = 1; // current_2i = 1

		while(current_power <= desired_power) {
			// TODO: Fixit.
			// current_modulo = p^(2k + 2^i)
			// pk2i = p^(k + 2^i)
			integer_matrix matrix(degree() + 1, degree() + 1);
			for(size_t i = 0; i < g.degree(); i++) {
				for(size_t j = 0; j <= h.degree(); j++) {
					matrix(i + j, i) = h[j];
				}
			}
			for(size_t i = 0; i <= h.degree(); i++) {
				for(size_t j = 0; j <= g.degree(); j++) {
					matrix(i + j, i + g.degree()) = g[j];
				}
			}

			vector b(degree() + 1);
			for(size_t i = 0; i < diff.coef.size(); i++) {
				mpz_divexact(b[i].get_mpz_t(), diff[i].get_mpz_t(), pk2i.get_mpz_t());
			}
			auto result = solve(std::move(matrix), std::move(b), current_modulo);

			std::vector<mpz_class> g_star_v(g.degree()), h_star_v(h.degree() + 1);
			for(size_t i = 0; i < g.degree(); i++) g_star_v[i] = result[i];
			for(size_t i = 0; i <= h.degree(); i++) h_star_v[i] = result[i + g.degree()];
			integer_polynomial g_star(std::move(g_star_v)), h_star(std::move(h_star_v));
			g += g_star * current_modulo; h += h_star * current_modulo;

			current_modulo *= tmp;
			pk2i *= tmp;
			tmp *= tmp;
			diff = *this - g * h;
			current_power += current_2i;
			current_2i <<= 1;
		}
		mpz_pow_ui(tmp.get_mpz_t(), p.get_mpz_t(), desired_power);
		for(size_t i = 0; i <= g.degree(); i++) {
			g[i] %= tmp;
			if(g[i] > ((tmp - 1) >> 1)) g[i] -= tmp;
			else if(g[i] < -((tmp - 1) >> 1)) g[i] += tmp;
		}
		for(size_t i = 0; i < h.degree(); i++) h[i] %= tmp;
		return { std::move(g), std::move(h) };
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

	integer_polynomial integer_polynomial::mul_scalar(mpz_class const &n) const {
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
			if(a.coef.empty()) return a;
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
		mpz_divexact(q.get_mpz_t(), a.leading().get_mpz_t(), b.leading().get_mpz_t());

		for(size_t i = 0; i < b.coef.size(); i++) {
			a.coef[a.coef.size() - b.coef.size() + i] -= b.coef[i] * q;
		}
		assert(a.coef.back() == 0);
		a.clean();

		if(!a.coef.empty()) {
			leading_gcd = a.leading();
			for(size_t i = 0; i < a.coef.size() - 1; i++) {
				leading_gcd = gcd(leading_gcd, a.coef[i]);
			} for(auto &now : a.coef) {
				mpz_divexact(q.get_mpz_t(), now.get_mpz_t(), leading_gcd.get_mpz_t());
				now = q;
			}
		}

		return gcd(std::move(b), std::move(a));
	}

	integer_polynomial gcd(integer_polynomial a, integer_polynomial b, const mpz_class &p) {
		assert(is_prime(p));
		for(auto &now : a.coef) now %= p;
		for(auto &now : b.coef) now %= p;
		a.clean(); b.clean();
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
			invert(inv, a.leading(), p);
			for(size_t i = 0; i < a.coef.size() - 1; i++) {
				a.coef[i] = (a.coef[i] * inv) % p;
			} a.coef.back() = 1;
			return a;
		}
		if(mpz_divisible_p(a.leading().get_mpz_t(), p.get_mpz_t())) {
			a.coef.pop_back();
			return gcd(std::move(a), std::move(b));
		}
		if(mpz_divisible_p(b.leading().get_mpz_t(), p.get_mpz_t())) {
			b.coef.pop_back();
			return gcd(std::move(a), std::move(b));
		}
		if(a.leading() < 0) return gcd(-a, std::move(b), p);
		if(b.leading() < 0) return gcd(std::move(a), -b, p);
		if(a.coef.size() < b.coef.size()) return gcd(std::move(b), std::move(a), p);

		mpz_class inv;
		invert(inv, b.leading(), p);
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

	mpz_class resultant(integer_polynomial const &f, integer_polynomial const &g) {
		return determinant(integer_matrix::sylvester_matrix(f, g));
	}
}
