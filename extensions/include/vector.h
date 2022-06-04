#ifndef POLYNOMIAL_PROJECT_VECTOR_H
#define POLYNOMIAL_PROJECT_VECTOR_H

#include <gmpxx.h>

#include <cassert>
#include <vector>
#include <ostream>

namespace project {

	class vector {
	public:
		explicit vector(size_t len);
		explicit vector(std::vector<mpz_class>);
		vector(std::initializer_list<std::string>);

		mpz_class &operator[](size_t idx) { assert(idx < length()); return vec[idx]; }
		mpz_class const &operator[](size_t idx) const { assert(idx < length()); return vec[idx]; }

		[[nodiscard]] size_t length() const { return vec.size(); }

		[[nodiscard]] bool is_zero() const;

		[[nodiscard]] mpz_class norm2() const;

		vector &reverse_eq();
		[[nodiscard]] vector reverse() const;

		vector &neg_eq();

		vector &operator+=(vector const &v) { return add_eq(v); }
		vector &operator-=(vector const &v) { return sub_eq(v); }
		vector &operator*=(mpz_class const &n) { return mul_scalar_eq(n); }

		[[nodiscard]] vector operator-() const { return neg(); }
		[[nodiscard]] vector operator+(vector const &v) const { return add(v); }
		[[nodiscard]] vector operator-(vector const &v) const { return sub(v); }
		[[nodiscard]] vector operator*(mpz_class const &n) const { return mul_scalar(n); }
		friend vector operator*(mpz_class const &n, vector const &v) { return v.mul_scalar(n); }

		friend mpz_class inner_product(vector const&, vector const&);
		friend std::ostream &operator<<(std::ostream&, vector const&);

	private:
		std::vector<mpz_class> vec;

		vector &add_eq(vector const&);
		vector &sub_eq(vector const&);
		vector &mul_scalar_eq(mpz_class const&);

		[[nodiscard]] vector neg() const;
		[[nodiscard]] vector add(vector const&) const;
		[[nodiscard]] vector sub(vector const&) const;
		[[nodiscard]] vector mul_scalar(mpz_class const&) const;
	};

	[[nodiscard]] bool is_linearly_independent(std::vector<vector>);

}

#endif
