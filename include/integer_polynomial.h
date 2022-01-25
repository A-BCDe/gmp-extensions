#ifndef POLYNOMIAL_PROJECT_INTEGER_POLYNOMIAL_H
#define POLYNOMIAL_PROJECT_INTEGER_POLYNOMIAL_H

#include <gmpxx.h>

#include <cassert>
#include <string>
#include <vector>

namespace project {
	class integer_polynomial {
	public:
		explicit integer_polynomial() = default;
		explicit integer_polynomial(std::vector<mpz_class>);
		explicit integer_polynomial(std::vector<std::string> const&);

		[[nodiscard]] size_t degree() const { assert(!is_zero()); return coef.size() - 1; }
		[[nodiscard]] bool is_zero() const noexcept { return coef.empty(); }
		integer_polynomial &negate() noexcept;

		// https://mathworld.wolfram.com/SylvesterMatrix.html
		// https://mathworld.wolfram.com/Resultant.html
		// https://mathworld.wolfram.com/PolynomialDiscriminant.html
		mpz_class discriminant() const;
		std::vector<integer_polynomial> factorize() const;

		mpz_class &operator[](size_t idx) { assert(idx < coef.size()); return coef[idx]; }
		mpz_class const &operator[](size_t idx) const { assert(idx < coef.size()); return coef[idx]; }

		mpz_class operator()(mpz_class const&) const;

		[[nodiscard]] integer_polynomial operator-() const { return neg(); }

		integer_polynomial &operator+=(integer_polynomial const &poly) { return add_eq(poly); }
		integer_polynomial &operator-=(integer_polynomial const &poly) { return sub_eq(poly); }
		integer_polynomial &operator*=(integer_polynomial const &poly) { return mul_eq(poly); }

		[[nodiscard]] integer_polynomial operator+(integer_polynomial const &poly) const { return add(poly); }
		[[nodiscard]] integer_polynomial operator-(integer_polynomial const &poly) const { return sub(poly); }
		[[nodiscard]] integer_polynomial operator*(integer_polynomial const &poly) const { return mul(poly); }

	private:
		std::vector<mpz_class> coef;

		void clean() noexcept;

		integer_polynomial &add_eq(integer_polynomial const&);
		integer_polynomial &sub_eq(integer_polynomial const&);
		integer_polynomial &mul_eq(integer_polynomial const&);

		[[nodiscard]] integer_polynomial neg() const;
		[[nodiscard]] integer_polynomial add(integer_polynomial const&) const;
		[[nodiscard]] integer_polynomial sub(integer_polynomial const&) const;
		[[nodiscard]] integer_polynomial mul(integer_polynomial const&) const;
	};
}

#endif