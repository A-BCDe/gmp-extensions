#include "vector.h"

namespace project {

	vector::vector(size_t len) : vec(len) {}

	vector::vector(std::vector<mpz_class> vec) : vec(std::move(vec)) {}

	vector::vector(std::initializer_list<std::string> list) {
		vec.reserve(list.size());
		for(auto const &now : list) vec.emplace_back(now);
	}

	mpz_class vector::norm2() const {
		mpz_class res = 0;
		for(auto const &now : vec) {
			res += now * now;
		} return res;
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
			vec[i] += v[i];
		} return *this;
	}

	vector const &vector::neg() const {
		vector v(*this);
		return v.neg_eq();
	}

	vector const &vector::add(vector const &v) const {
		vector vec(*this);
		return vec.add_eq(v);
	}

	vector const &vector::sub(vector const &v) const {
		vector vec(*this);
		return vec.sub_eq(v);
	}

	// friend functions

	std::ostream &operator<<(std::ostream &os, vector const& v) {
		if(!v.length()) return os;
		os << v[0].get_str();
		for(size_t i = 1; i < v.length(); i++) {
			os << ' ' << v[i].get_str();
		} return os;
	}


}