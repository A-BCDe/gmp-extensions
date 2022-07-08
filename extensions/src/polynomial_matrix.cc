#include "polynomial_matrix.h"

#include <ostream>

namespace gmp_extensions {

	polynomial_matrix::polynomial_matrix(size_t row, size_t col)
	  : row(row), col(col), mat(row * col) {}

	integer_polynomial &polynomial_matrix::operator()(size_t r, size_t c) {
		assert(r < row);
		assert(c < col);
		return mat[r * col + c];
	}

	integer_polynomial const &polynomial_matrix::operator()(size_t r, size_t c) const {
		assert(r < row);
		assert(c < col);
		return mat[r * col + c];
	}

	std::ostream &operator<<(std::ostream &os, polynomial_matrix const &M) {
		std::vector<std::string> out;
		size_t size = 0;
		for(size_t r = 0; r < M.row; r++) {
			for(size_t c = 0; c < M.col; c++) {
				out.push_back(M(r, c).get_str());
				if(out.back().size() > size) size = out.back().size();
			}
		}
		for(size_t r = 0; r < M.row; r++) {
			for(size_t c = 0; c < M.col; c++) {
				os << std::string(size - out[r * M.col + c].size() + 2, ' ') << out[r * M.col + c];
				os << ' ';
			}
			if(r < M.row - 1) os << '\n';
		}
		return os;
	}


}
