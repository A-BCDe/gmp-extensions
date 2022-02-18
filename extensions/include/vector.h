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

		size_t length() const { return vec.size(); }

		mpz_class norm2() const;

		vector &neg_eq();
		vector &add_eq(vector const&);
		vector &sub_eq(vector const&);

		vector neg() const;
		vector add(vector const&) const;
		vector sub(vector const&) const;

		vector &operator+=(vector const &v) { return add_eq(v); }
		vector &operator-=(vector const &v) { return sub_eq(v); }

		vector operator-() const { return neg(); }
		vector operator+(vector const &v) const { return add(v); }
		vector operator-(vector const &v) const { return sub(v); }

		friend std::ostream &operator<<(std::ostream&, vector const&);

	private:
		std::vector<mpz_class> vec;
	};

}

#endif
