#include <matrix.h>

#include <gmpxx.h>

#include <number_theoretic.h>

#include <cassert>
#include <initializer_list>
#include <iostream>
#include <ostream>
#include <vector>
#include <unistd.h>

namespace project {

// constructors

    integer_matrix::integer_matrix(size_t row, size_t col)
            : row(row), col(col), mat(row * col) {}

    integer_matrix::integer_matrix(size_t row, size_t col, std::vector<std::string> const &data)
            : row(row), col(col), mat(row * col) {
        assert(row * col == data.size());
        for(size_t i = 0; i < row * col; i++) {
            mat[i] = data[i];
        }
    }

    integer_matrix::integer_matrix(size_t row, size_t col, std::initializer_list<mpz_class> list)
            : row(row), col(col), mat(list) {
        assert(list.size() >= row * col);
    }

	integer_matrix integer_matrix::identity(size_t row, size_t col) {
		integer_matrix I(row, col);
		for(size_t i = 0; i < row && i < col; i++) I.mat[i * col + i] = 1;
		return I;
	}

	integer_matrix integer_matrix::sylvester_matrix(integer_polynomial const &f, integer_polynomial const &g) {
		assert(!f.is_zero());
		assert(!g.is_zero());
		integer_matrix matrix(f.degree() + g.degree(), f.degree() + g.degree());
		size_t r = 0;
		for(; r < g.degree(); r++) {
			for(size_t c = 0; c <= f.degree(); c++) {
				matrix(r, r + c) = f[f.degree() - c];
			}
		}
		for(; r < f.degree() + g.degree(); r++) {
			for(size_t c = 0; c <= g.degree(); c++) {
				matrix(r, r + c - g.degree()) = g[g.degree() - c];
			}
		}
		return matrix;
	}

// operators

    mpz_class &integer_matrix::operator()(size_t r, size_t c) {
        assert(r < row);
        assert(c < col);
        return mat[r * col + c];
    }

    mpz_class const &integer_matrix::operator()(size_t r, size_t c) const {
        assert(r < row);
        assert(c < col);
        return mat[r * col + c];
    }

// member functions

    mpz_class &integer_matrix::get(size_t r, size_t c) {
        assert(r < row);
        assert(c < col);
        return mat[r * col + c];
    }

    mpz_class const &integer_matrix::get(size_t r, size_t c) const {
        assert(r < row);
        assert(c < col);
        return mat[r * col + c];
    }

	std::vector<vector> integer_matrix::to_row_vectors() const {
		std::vector<vector> result;
		result.reserve(row);
		for(size_t r = 0; r < row; r++) {
			vector v(col);
			for(size_t c = 0; c < col; c++) {
				v[c] = get(r, c);
			}
			result.emplace_back(std::move(v));
		}
		return result;
	}

	std::vector<vector> integer_matrix::to_col_vectors() const {
		std::vector<vector> result;
		result.reserve(col);
		for(size_t c = 0; c < col; c++) {
			vector v(row);
			for(size_t r = 0; r < row; r++) {
				v[r] = get(r, c);
			}
			result.emplace_back(std::move(v));
		}
		return result;
	}

    integer_matrix integer_matrix::to_smith_normal_form() {
        auto Y = identity(row, col);
        auto current_row = row;
        auto current_col = col;
        size_t piv = 0;
        while(piv < current_row && piv < current_col) {
            while(piv < current_row && piv < current_col) {
                size_t i;
                // exchange (piv, *) and (i, *) when (i, piv) != 0
                for(i = piv; i < current_col; i++) {
                    if(get(i, piv) != 0) {
                        if(get(i, piv) < 0) {
                            row_negate(i);
                        }
                        row_exchange(piv, i);
                        break;
                    }
                }
                if(i == current_col) {
                    // replace zero column with another column.
                    col_exchange_with_companion_matrix(Y, piv, --current_col);
                    continue;
                }
                // GCD operation on (*, piv)
                for(i = piv + 1; i < current_row; i++) {
                    row_gcd_operation(piv, piv, i);
                }

                // exchange (*, piv) and (*, i) when (piv, i) != 0
                for(i = piv; i < current_row; i++) {
                    if(get(piv, i) != 0) {
                        if(get(i, piv) < 0) {
                            col_negate_with_companion_matrix(Y, i);
                        }
                        col_exchange_with_companion_matrix(Y, piv, i);
                        break;
                    }
                }
                if(i == current_row) {
                    row_exchange(piv, --current_row);
                }
                // GCD operation on (piv, *)
                for(size_t i = piv + 1; i < current_col; i++) {
                    col_gcd_operation_with_companion_matrix(Y, piv, piv, i);
                }

                // check (piv, *)
                for(i = piv + 1; i < current_col; i++) {
                    if(get(piv, i) != 0) goto CON;
                }
                // check (*, piv)
                for(i = piv + 1; i < current_row; i++) {
                    if(get(i, piv) != 0) goto CHK;
                }
                break;

                CHK:if(get(piv, piv) != 0) {
                for(i = piv + 1; i < current_row; i++) {
                    row_multiply_and_add(piv, i, get(i, piv) / get(piv, piv));
                }
            }
                CON:continue;
            }
            piv++;
        }
        while(true) {
            bool flag = false;
            for(size_t i = 0; i < current_row - 1 && i < current_col - 1; i++) {
                if(get(i + 1, i + 1) % get(i, i) == 0) continue;
                flag = true;
                mpz_class const a = get(i, i), b = get(i + 1, i + 1);
                auto const p = extended_euclidean_algorithm(a, b);
                mpz_class const g = a * p.first + b * p.second;
                // ax + by = g
                //
                // | a 0 |
                // | 0 b |
                row_multiply_and_add(i, i + 1, p.first);
                // | a  0 |
                // | ax b |
                col_multiply_and_add_with_companion_matrix(Y, i + 1, i, p.second);
                // | a 0 |
                // | g b |
                row_multiply_and_add(i + 1, i, 1 - a / g);
                // | g (1 - a / g)b |
                // | g b            |
                col_multiply_and_add_with_companion_matrix(Y, i, i + 1, (a / g - 1) * (b / g));
                // | g 0      |
                // | g ab / g |
                row_multiply_and_add(i, i + 1, -1);
                // | g 0      |
                // | 0 ab / g |
            }
            if(!flag) break;
        }
        return Y;
    }

    std::pair<integer_matrix, integer_matrix> integer_matrix::smith_normal_form() const {
        integer_matrix matrix(*this);
        integer_matrix Y(matrix.to_smith_normal_form());
        return { std::move(matrix), std::move(Y) };
    }

    mpz_class integer_matrix::to_congruent_diagonal() {
        assert(is_symmetric());
        size_t piv = 0;
        mpz_class denominator = 1, g;
        while(piv < row) {
            mpz_class const &st = get(piv, piv);
            if(st == 0) {
                for(size_t i = piv + 1; i < row; i++) {
                    if(get(piv, i) != 0) {
                        row_multiply_and_add(i, piv, 1);
                        col_multiply_and_add(i, piv, 1);
                        if(st == 0) {
                            row_multiply_and_add(i, piv, 1);
                            col_multiply_and_add(i, piv, 1);
                        }
                        break;
                    }
                }
            }
            assert(st != 0);

            mpz_class const st2 = st * st;

            denominator *= st2;

            for(size_t i = piv + 1; i < row; i++) {
                for(size_t j = piv + 1; j < col; j++) {
                    get(i, j) *= st2;
                    get(i, j) -= st * get(piv, i) * get(piv, j);
                }
                get(i, piv) = 0;
            }
            for(size_t i = piv + 1; i < row; i++) get(piv, i) = 0;
            for(size_t i = 0; i <= piv; i++) get(i, i) *= st2;

            g = denominator;

            for(size_t i = 0; i <= piv; i++) g = gcd(g, get(i, i));
            for(size_t i = piv + 1; i < row; i++) {
                for(size_t j = i; j < col; j++) {
                    g = gcd(g, get(i, j));
                }
            }

            denominator /= g;
            for(size_t i = 0; i <= piv; i++) get(i, i) /= g;
            for(size_t i = piv + 1; i < row; i++) {
                for(size_t j = piv + 1; j < row; j++) {
                    get(i, j) /= g;
                }
            }

            piv++;
        }
        return denominator;
    }

    std::pair<integer_matrix, mpz_class> integer_matrix::congruent_diagonal() const {
        integer_matrix matrix(*this);
        mpz_class denominator = matrix.to_congruent_diagonal();
        return { std::move(matrix), std::move(denominator) };
    }

    std::pair<mpz_class, mpz_class> integer_matrix::signature() const {
        auto cd = congruent_diagonal();
        auto D = cd.first; auto den = cd.second;
        mpz_class pos(0), neg(0);
        bool sign = (den > 0);
        for(size_t i = 0; i < row; i++) {
            pos += ((!sign) ^ (D(i, i) > 0));
            neg += ((!sign) ^ (D(i, i) < 0));
        }
        return { std::move(pos), std::move(neg) };
    }

    // (l, *) <-> (r, *)
    void integer_matrix::row_exchange(size_t l, size_t r) {
        assert(l < row);
        assert(r < row);
        if(l == r) return;
        auto const L = l * col, R = r * col;
        for(size_t i = 0; i < col; i++) {
            std::swap(mat[L + i], mat[R + i]);
        }
    }

    // (r, *) <- -(r, *)
    void integer_matrix::row_negate(size_t r) {
        assert(r < row);
        size_t const R = r * col;
        for(size_t i = R; i < R + col; i++) {
            mat[i] = -mat[i];
        }
    }

	// (r, *) <- n * (r, *)
	void integer_matrix::row_multiply(size_t r, mpz_class const &n) {
		assert(r < row);
		size_t const R = r * col;
		for(size_t i = R; i < R + col; i++) {
			mat[i] *= n;
		}
	}

    // (to, *) <- mul * (from, *) + (to, *)
    void integer_matrix::row_multiply_and_add(size_t from, size_t to, mpz_class const &mul) {
        assert(from < row);
        assert(to < row);
        for(size_t i = 0; i < col; i++) {
            get(to, i) += mul * get(from, i);
        }
    }

    // calculate gcd of (l, piv) and (r, piv)
    void integer_matrix::row_gcd_operation(size_t piv, size_t l, size_t r) {
        assert(l < row);
        assert(r < row);
        if(get(std::max(l, r), piv) == 0) return;
        if(get(std::min(l, r), piv) == 0) {
            row_exchange(l, r);
            return;
        }
        row_multiply_and_add(r, l, -(get(l, piv) / get(r, piv)));
        row_gcd_operation(piv, r, l);
    }

    // (*, l) <-> (*, r)
    void integer_matrix::col_exchange_with_companion_matrix(integer_matrix &M, size_t l, size_t r) {
        col_exchange(l, r);
        M.row_exchange(l, r);
    }

    // (*, c) <- -(*, c)
    void integer_matrix::col_negate_with_companion_matrix(integer_matrix &M, size_t c) {
        col_negate(c);
        M.row_negate(c);
    }

    // (*, to) <- mul * (*, from) + (*, to)
    void integer_matrix::col_multiply_and_add_with_companion_matrix(integer_matrix &M, size_t from, size_t to, mpz_class const &mul) {
        col_multiply_and_add(from, to, mul);
        M.row_multiply_and_add(to, from, -mul);
    }

    // calculate gcd of (piv, l) and (piv, r)
    void integer_matrix::col_gcd_operation_with_companion_matrix(integer_matrix &M, size_t piv, size_t l, size_t r) {
        assert(l < col);
        assert(r < col);
        if(get(piv, std::max(l, r)) == 0) return;
        if(get(piv, std::min(l, r)) == 0) {
            col_exchange_with_companion_matrix(M, l, r);
            return;
        }
        col_multiply_and_add_with_companion_matrix(M, r, l, -(get(piv, l) / get(piv, r)));
        col_gcd_operation_with_companion_matrix(M, piv, r, l);
    }

    // (*, l) <-> (*, r)
    void integer_matrix::col_exchange(size_t l, size_t r) {
        assert(l < col);
        assert(r < col);
        if(l == r) return;
        for(size_t i = 0; i < row * col; i += col) {
            std::swap(mat[i + l], mat[i + r]);
        }
    }

    // (*, c) <- -(*, c)
    void integer_matrix::col_negate(size_t c) {
        assert(c < col);
        for(size_t i = c; i < row * col; i += col) {
            mat[i] = -mat[i];
        }
    }

	// (*, c) <- n * (*, c)
	void integer_matrix::col_multiply(size_t c, mpz_class const &n) {
		assert(c < col);
		for(size_t i = c; i < row * col; i += col) {
			mat[i] *= n;
		}
	}

    // (*, to) <- mul * (*, from) + (*, to)
    void integer_matrix::col_multiply_and_add(size_t from, size_t to, mpz_class const &mul) {
        assert(from < col);
        assert(to < col);
        for(size_t i = 0; i < row * col; i += col) {
            mat[i + to] += mul * mat[i + from];
        }
    }

    // calculate gcd of (piv, l) and (piv, r)
    void integer_matrix::col_gcd_operation(size_t piv, size_t l, size_t r) {
        assert(l < col);
        assert(r < col);
        if(get(piv, std::max(l, r)) == 0) return;
        if(get(piv, std::min(l, r)) == 0) {
            col_exchange(l, r);
            return;
        }
        col_multiply_and_add(r, l, -(get(piv, l) / get(piv, r)));
        col_gcd_operation(piv, r, l);
    }

// misc

    bool integer_matrix::is_symmetric() const {
        if(row != col) return false;
        for(size_t i = 1; i < row; i++) {
            for(size_t j = i + 1; j < col; j++) {
                if(get(i, j) != get(j, i)) return false;
            }
        }
        return true;
    }

    std::vector<mpz_class> integer_matrix::get_diagonal() const {
        std::vector<mpz_class> v;
        for(size_t i = 0; i < row && i < col; i++) {
            v.push_back(get(i, i));
        }
        return v;
    }

// friends

    std::ostream &operator<<(std::ostream &os, integer_matrix const &M) {
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
                os << std::string(size - out[r * M.col + c].size() + 2, ' ') << out[r * M.col + c] << ' ';
            }
            if(r < M.row - 1) os << '\n';
        }
        return os;
    }

	vector solve(integer_matrix matrix, vector v, mpz_class const &m) {
		assert(m > 0);
		mpz_class inv, tmp;
		for(size_t c = 0; c < matrix.col; c++) {
			for(size_t r = c; r < matrix.row; r++) {
				if(gcd(matrix(r, c), m) == 1) {
					matrix.row_exchange(r, c);
					std::swap(v[r], v[c]);
					goto NXT;
				}
			}
			assert(0);
		NXT:invert(inv, matrix(c, c), m);
			matrix.row_multiply(c, inv);
			v[c] *= inv;
			matrix(c, c) = 1;
			for(size_t r = 0; r < matrix.row; r++) {
				if(r == c) continue;
				tmp = -matrix(r, c);
				matrix.row_multiply_and_add(c, r, tmp);
				v[r] += tmp * v[c];
			}
		}
		for(size_t i = 0; i < v.length(); i++) {
			v[i] = (v[i] % m + m) % m;
		}
		return v;
	}

	std::vector<vector> kernel(integer_matrix mat) {
		size_t piv = 0;
		std::vector<vector> space;
		std::vector<size_t> pivot;
		mpz_class g, mul, tmp, tmp2;
		for(size_t c = 0; c < mat.col; c++) {
			size_t r = piv;
			while(r < mat.row) {
				if(mat(r, c) != 0) {
					pivot.push_back(c);
					mat.row_exchange(piv, r);
					for(r = 0; r < mat.row; r++) {
						if(r == piv) continue;
						if(mat(r, c) == 0) continue;
						g = gcd(mat(piv, c), mat(r, c));
						mpz_divexact(tmp.get_mpz_t(), mat(piv, c).get_mpz_t(), g.get_mpz_t());
						mat.row_multiply(r, tmp);
						mpz_divexact(tmp.get_mpz_t(), mat(r, c).get_mpz_t(), mat(piv, c).get_mpz_t());
						mat.row_multiply_and_add(piv, r, -tmp);
					}
					piv++;
					break;
				}
				r++;
			}
		}

		for(size_t r = 0; r < mat.row; r++) {
			g = abs(mat(r, 0));
			for(size_t c = 1; c < mat.col; c++) g = gcd(g, mat(r, c));
			if(g == 0) continue;
			for(size_t c = 0; c < mat.col; c++) {
				mpz_divexact(tmp.get_mpz_t(), mat(r, c).get_mpz_t(), g.get_mpz_t());
				mat(r, c) = tmp;
			}
		}

		space.reserve(mat.col - pivot.size());
		if(pivot.empty()) {
			for(size_t c = 0; c < mat.col; c++) {
				vector v(mat.col);
				v[c] = 1;
				space.emplace_back(std::move(v));
			}
			return space;
		}

		piv = 0;
		mul = 1;
		for(size_t r = 0; r < pivot.size(); r++) {
			mul *= mat(r, pivot[r]);
		}
		for(size_t c = 0; c < mat.col; c++) {
			bool flag = false;
			if(pivot[piv] == c) {
				piv++;
				continue;
			}
			vector v(mat.col);
			v[c] = mul;
			for(size_t r = 0; r < pivot.size(); r++) {
				v[pivot[r]] = v[c] * mat(r, c);
				mpz_divexact(tmp2.get_mpz_t(), v[pivot[r]].get_mpz_t(), mat(r, pivot[r]).get_mpz_t());
				v[pivot[r]] = tmp2;
			}
			g = v[0];
			for(size_t i = 1; i < mat.col; i++) {
				g = gcd(g, v[i]);
			}
			assert(g != 0);
			for(size_t i = 0; i < mat.col; i++) {
				if(!flag) {
					if(v[i] > 0) flag = true;
					else if(v[i] < 0) flag = true, g = -g;
					else continue;
				}
				mpz_divexact(tmp2.get_mpz_t(), v[i].get_mpz_t(), g.get_mpz_t());
				v[i] = tmp2;
			}
			space.emplace_back(std::move(v));
		}

		return space;
	}

	std::vector<vector> kernel(integer_matrix mat, mpz_class const &p) {
		assert(p != 0);
		size_t piv = 0;
		std::vector<vector> space;
		std::vector<size_t> pivot;
		mpz_class g, mul, tmp, tmp2;
		for(size_t c = 0; c < mat.col; c++) {
			size_t r = piv;
			while(r < mat.row) {
				mat(r, c) %= p;
				if(mat(r, c) != 0) {
					pivot.push_back(c);
					mat.row_exchange(piv, r);
					for(r = 0; r < mat.row; r++) {
						if(r == piv) continue;
						if(mat(r, c) == 0) continue;
						g = gcd(mat(piv, c), mat(r, c));
						mpz_divexact(tmp.get_mpz_t(), mat(piv, c).get_mpz_t(), g.get_mpz_t());
						mat.row_multiply(r, tmp);
						mpz_divexact(tmp.get_mpz_t(), mat(r, c).get_mpz_t(), mat(piv, c).get_mpz_t());
						mat.row_multiply_and_add(piv, r, -tmp);
					}
					piv++;
					break;
				}
				r++;
			}
		}

		for(size_t r = 0; r < mat.row; r++) {
			for(size_t c = 0; c < mat.col; c++) {
				mat(r, c) %= p;
			}
		}

		for(size_t r = 0; r < mat.row; r++) {
			g = abs(mat(r, 0));
			for(size_t c = 1; c < mat.col; c++) g = gcd(g, mat(r, c));
			if(g == 0) continue;
			for(size_t c = 0; c < mat.col; c++) {
				mpz_divexact(tmp.get_mpz_t(), mat(r, c).get_mpz_t(), g.get_mpz_t());
				mat(r, c) = tmp;
			}
		}

		space.reserve(mat.col - pivot.size());
		if(pivot.empty()) {
			for(size_t c = 0; c < mat.col; c++) {
				vector v(mat.col);
				v[c] = 1;
				space.emplace_back(std::move(v));
			}
			return space;
		}

		piv = 0;
		mul = 1;
		for(size_t r = 0; r < pivot.size(); r++) {
			mul *= mat(r, pivot[r]);
		}
		assert(mul % p != 0);
		for(size_t c = 0; c < mat.col; c++) {
			bool flag = false;
			if(pivot[piv] == c) {
				piv++;
				continue;
			}
			vector v(mat.col);
			v[c] = mul;
			for(size_t r = 0; r < pivot.size(); r++) {
				tmp2 = -v[c] * mat(r, c);
				mpz_divexact(v[pivot[r]].get_mpz_t(), tmp2.get_mpz_t(), mat(r, pivot[r]).get_mpz_t());
			}
			g = v[0];
			for(size_t i = 1; i < mat.col; i++) {
				g = gcd(g, v[i]);
			}
			assert(g != 0);
			for(size_t i = 0; i < mat.col; i++) {
				if(!flag) {
					if(v[i] > 0) flag = true;
					else if(v[i] < 0) flag = true, g = -g;
					else continue;
				}
				mpz_divexact(tmp2.get_mpz_t(), v[i].get_mpz_t(), g.get_mpz_t());
				v[i] = (tmp2 % p + p) % p;
			}
			space.emplace_back(std::move(v));
		}

		return space;
	}

	mpz_class determinant(integer_matrix matrix) {
		assert(matrix.row == matrix.col);
		mpz_class det(1); // determinant
		mpz_class denom(1); // denominator
		mpz_class g, tmp;
		for(size_t r = 0; r < matrix.row; r++) {
			for(size_t c = r; c < matrix.col; c++) {
				if(matrix(r, c) != 0) {
					matrix.col_exchange(r, c);
					goto NXT;
				}
			}
			return 0;
			NXT:det *= matrix(r, r);
			for(size_t c = r + 1; c < matrix.col; c++) {
				g = gcd(matrix(r, r), matrix(r, c));
				mpz_divexact(tmp.get_mpz_t(), matrix(r, r).get_mpz_t(), g.get_mpz_t());
				matrix.col_multiply(c, tmp);
				denom *= tmp;
				mpz_divexact(tmp.get_mpz_t(), matrix(r, c).get_mpz_t(), matrix(r, r).get_mpz_t());
				matrix.col_multiply_and_add(r, c, -tmp);
			}
		}
		assert(mpz_divisible_p(det.get_mpz_t(), denom.get_mpz_t()));
		mpz_divexact(tmp.get_mpz_t(), det.get_mpz_t(), denom.get_mpz_t());
		return tmp;
	}
}
