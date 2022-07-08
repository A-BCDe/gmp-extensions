#include <vector.h>

namespace gmp_extensions {

	vector::vector(size_t len) : vec(len) {}

	vector::vector(std::vector<mpz_class> vec) : vec(std::move(vec)) {}

	vector::vector(std::initializer_list<std::string> list) {
		vec.reserve(list.size());
		for(auto const &now : list) vec.emplace_back(now);
	}

	bool vector::is_zero() const {
		return std::all_of(vec.begin(), vec.end(), [](mpz_class const &n) {
			return n == 0;
		});
	}

	mpz_class vector::norm2() const {
		mpz_class res = 0;
		for(auto const &now : vec) {
			res += now * now;
		} return res;
	}

	vector &vector::reverse_eq() {
		std::reverse(vec.begin(), vec.end());
		return *this;
	}

	vector vector::reverse() const {
		vector v(*this);
		return v.reverse_eq();
	}

	vector &vector::neg_eq() {
		for(auto &now : vec) {
			now = -now;
		} return *this;
	}

	vector &vector::add_eq(vector const &v) {
		assert(length() == v.length());
		for(size_t i = 0; i < length(); i++) {
			vec[i] += v[i];
		} return *this;
	}

	vector &vector::sub_eq(vector const &v) {
		assert(length() == v.length());
		for(size_t i = 0; i < length(); i++) {
			vec[i] -= v[i];
		} return *this;
	}

	vector &vector::mul_scalar_eq(mpz_class const &n) {
		for(auto &now : vec) {
			now *= n;
		} return *this;
	}

	vector vector::neg() const {
		vector v(*this);
		return v.neg_eq();
	}

	vector vector::add(vector const &v) const {
		vector vec(*this);
		return vec.add_eq(v);
	}

	vector vector::sub(vector const &v) const {
		vector vec(*this);
		return vec.sub_eq(v);
	}

	vector vector::mul_scalar(mpz_class const &n) const {
		vector vec(*this);
		return vec.mul_scalar_eq(n);
	}

	// friend functions

	mpz_class inner_product(vector const &v1, vector const &v2) {
		assert(v1.length() == v2.length());
		mpz_class res(0);
		for(size_t i = 0; i < v1.length(); i++) {
			res += v1[i] * v2[i];
		} return res;
	}

	std::ostream &operator<<(std::ostream &os, vector const& v) {
		if(!v.length()) return os;
		os << v[0].get_str();
		for(size_t i = 1; i < v.length(); i++) {
			os << ' ' << v[i].get_str();
		} return os;
	}

	// outside class

	bool is_linearly_independent(std::vector<vector> vectors) {
		if(vectors.empty()) return true;
		if(vectors[0].length() < vectors.size()) return false;
#ifndef NDEBUG
		for(size_t i = 1; i < vectors.size(); i++) {
			assert(vectors[i].length() == vectors[0].length());
		}
#endif
		auto const row = vectors.size();
		auto const col = vectors[0].length();
		size_t piv = 0;
		mpz_class g, tmp;
		for(size_t i = 0; i < col; i++) {
			for(size_t j = piv; j < row; j++) {
				if(vectors[j][i] != 0) {
					std::swap(vectors[j], vectors[piv]);
					goto NXT;
				}
			} continue;

		NXT:for(size_t j = 0; j < row; j++) {
				if(j == piv) continue;
				g = gcd(vectors[j][i], vectors[piv][i]);
				mpz_divexact(tmp.get_mpz_t(), vectors[piv][i].get_mpz_t(), g.get_mpz_t());
				vectors[j] *= tmp;
				mpz_divexact(tmp.get_mpz_t(), vectors[j][i].get_mpz_t(), vectors[piv][i].get_mpz_t());
				vectors[j][i] -= vectors[piv][i] * tmp;
			}
			piv++;
		}

		for(size_t i = 0; i < row; i++) {
			if(vectors[i].is_zero()) return false;
		} return true;
	}
}