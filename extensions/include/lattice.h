#ifndef POLYNOMIAL_PROJECT_LATTICE_H
#define POLYNOMIAL_PROJECT_LATTICE_H

#include <gmpxx.h>

#include <vector.h>
#include <matrix.h>

#include <vector>

namespace project {

	class lattice {
	public:
		lattice() = delete;
		explicit lattice(std::vector<vector>);
		//explicit lattice(integer_matrix const&);

		void LLL(mpz_class const&, mpz_class const&);

		void print(std::ostream &os) const {
			os << bases[0];
			for(size_t i = 1; i < bases.size(); i++) {
				os << '\n' << bases[i];
			}
		}

		vector short_vector() const { return bases[0]; }

	private:
		std::vector<vector> bases;

		static void LLL(std::vector<vector>&, mpz_class const&, mpz_class const&);

		static integer_matrix GSO_computation(std::vector<vector> const&);
	};

}

#endif