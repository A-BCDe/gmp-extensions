#ifndef POLYNOMIAL_PROJECT_INTEGER_POLYNOMIAL_H
#define POLYNOMIAL_PROJECT_INTEGER_POLYNOMIAL_H

#include <gmpxx.h>

#include <vector.h>

#include <cassert>
#include <string>
#include <vector>

namespace project {
	class integer_polynomial {
	public:
		explicit integer_polynomial() = default;
		explicit integer_polynomial(std::vector<mpz_class>);
		explicit integer_polynomial(std::vector<std::string> const&);
		explicit integer_polynomial(vector const&);

		[[nodiscard]] static integer_polynomial monomial(mpz_class const&, size_t);
		[[nodiscard]] static integer_polynomial constant(mpz_class const&);

		[[nodiscard]] size_t degree() const noexcept { assert(!is_zero()); return coef.size() - 1; }
		[[nodiscard]] mpz_class const &leading() const { assert(!is_zero()); return coef.back(); }

		[[nodiscard]] bool is_zero() const noexcept { return coef.empty(); }
		[[nodiscard]] bool is_constant() const noexcept { return coef.size() < 2; }
		[[nodiscard]] bool is_square_free() const;
		[[nodiscard]] bool is_square_free(mpz_class const &p) const;
		[[nodiscard]] bool is_monic() const noexcept { return !coef.empty() && coef.back() == 1; }
		[[nodiscard]] bool is_divisible_modulo(integer_polynomial const&, mpz_class const&) const;
		[[deprecated("You might want to use divexact instead.\n"), nodiscard]]
		bool is_divisible(integer_polynomial const&) const;

		integer_polynomial &negate() noexcept;

		integer_polynomial &modulo_eq(integer_polynomial, mpz_class const&);
		[[nodiscard]] integer_polynomial modulo(integer_polynomial, mpz_class const&) const;

		integer_polynomial &primitive_eq();
		[[nodiscard]] integer_polynomial primitive() const;

		integer_polynomial &derivative_eq();
		[[nodiscard]] integer_polynomial derivative() const;

		[[nodiscard]] integer_polynomial *divexact(integer_polynomial const&) const;

		integer_polynomial &divexact_modulo_eq(integer_polynomial const&, mpz_class const&);
		[[nodiscard]] integer_polynomial divexact_modulo(integer_polynomial const &poly, mpz_class const &p) const;

		std::vector<integer_polynomial> factorize() const;
		std::vector<integer_polynomial> factorize(mpz_class const&) const;

		std::pair<integer_polynomial, integer_polynomial> hensel_lifting(integer_polynomial, integer_polynomial, mpz_class const&, size_t, size_t) const;

		mpz_class &operator[](size_t idx) { assert(idx < coef.size()); return coef[idx]; }
		mpz_class const &operator[](size_t idx) const { assert(idx < coef.size()); return coef[idx]; }

		mpz_class operator()(mpz_class const&) const;

		[[nodiscard]] integer_polynomial operator-() const { return neg(); }

		integer_polynomial &operator+=(integer_polynomial const &poly) { return add_eq(poly); }
		integer_polynomial &operator-=(integer_polynomial const &poly) { return sub_eq(poly); }
		integer_polynomial &operator*=(integer_polynomial const &poly) { return mul_eq(poly); }
		integer_polynomial &operator*=(mpz_class const &n) { return mul_scalar_eq(n); }

		[[nodiscard]] integer_polynomial operator+(integer_polynomial const &poly) const { return add(poly); }
		[[nodiscard]] integer_polynomial operator-(integer_polynomial const &poly) const { return sub(poly); }
		[[nodiscard]] integer_polynomial operator*(integer_polynomial const &poly) const { return mul(poly); }
		[[nodiscard]] integer_polynomial operator*(mpz_class const &n) const { return mul_scalar(n); }
		[[nodiscard]] friend integer_polynomial operator*(mpz_class const &n, integer_polynomial const &poly) { return poly * n; }

		[[nodiscard]] bool operator==(integer_polynomial const &poly) const { return equal(poly); }
		[[nodiscard]] bool operator!=(integer_polynomial const &poly) const { return !equal(poly); }

		std::string get_str() const;

		friend integer_polynomial gcd(integer_polynomial, integer_polynomial);
		friend integer_polynomial gcd(integer_polynomial, integer_polynomial, mpz_class const&);

		friend std::ostream &operator<<(std::ostream&, integer_polynomial const&);

	private:
		std::vector<mpz_class> coef;

		void clean() noexcept;

		integer_polynomial &add_eq(integer_polynomial const&);
		integer_polynomial &sub_eq(integer_polynomial const&);
		integer_polynomial &mul_eq(integer_polynomial const&);
		integer_polynomial &mul_scalar_eq(mpz_class const&);

		[[nodiscard]] integer_polynomial neg() const;
		[[nodiscard]] integer_polynomial add(integer_polynomial const&) const;
		[[nodiscard]] integer_polynomial sub(integer_polynomial const&) const;
		[[nodiscard]] integer_polynomial mul(integer_polynomial const&) const;
		[[nodiscard]] integer_polynomial mul_scalar(mpz_class const&) const;

		[[nodiscard]] bool equal(integer_polynomial const&) const;
	};

	mpz_class resultant(integer_polynomial const&, integer_polynomial const&);
}

#endif