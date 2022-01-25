#include <integer_polynomial.h>

namespace project {

	integer_polynomial::integer_polynomial(std::vector<mpz_class> coef)
	  : coef(std::move(coef)) {}

	integer_polynomial::integer_polynomial(std::vector<std::string> const &coef) {
		this->coef.reserve(coef.size());
		for(std::string const &now : coef) {
			this->coef.emplace_back(now.c_str());
		}
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

	mpz_class integer_polynomial::operator()(const mpz_class &x) const {
		mpz_class res(coef[0]), mul(x);
		for(size_t i = 1; i < coef.size(); i++) {
			res += mul * coef[i];
			mul *= x;
		} return res;
	}

	integer_polynomial &integer_polynomial::negate() noexcept {
		for(auto &now : coef) now = -now;
		return *this;
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
		std::vector<mpz_class> tmp(std::move(coef));
		coef.clear();
		coef.resize(tmp.size() + poly.degree());
		for(size_t i = 0; i < tmp.size(); i++) {
			for(size_t j = 0; j <= poly.degree(); j++) {
				coef[i + j] += tmp[i] * poly.coef[j];
			}
		}
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
}