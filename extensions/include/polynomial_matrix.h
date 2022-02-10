#ifndef POLYNOMIAL_PROJECT_POLYNOMIAL_MATRIX_H
#define POLYNOMIAL_PROJECT_POLYNOMIAL_MATRIX_H

#include <polynomial.h>

namespace project {

	class polynomial_matrix {
	public:
		polynomial_matrix(size_t, size_t);

		integer_polynomial &operator()(size_t r, size_t c);
		integer_polynomial const &operator()(size_t r, size_t c) const;

		friend std::ostream &operator<<(std::ostream &os, polynomial_matrix const &M);

	private:
		size_t row, col;
		std::vector<integer_polynomial> mat;
	};

}

#endif
