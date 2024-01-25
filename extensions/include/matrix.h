#ifndef INTEGER_MATRIX_INTEGER_MATRIX_H
#define INTEGER_MATRIX_INTEGER_MATRIX_H

#include <gmpxx.h>

#include <polynomial.h>
#include <vector.h>

#include <initializer_list>
#include <vector>

namespace gmp_extensions {
    class integer_matrix {
    public:
        integer_matrix(size_t row, size_t col);
        integer_matrix(size_t row, size_t col, std::vector<std::string> const &data);
        integer_matrix(size_t row, size_t col, std::initializer_list<mpz_class> list);

        [[nodiscard]]
        static integer_matrix identity(size_t row, size_t col);
        [[nodiscard]]
        static integer_matrix sylvester_matrix(integer_polynomial const &, integer_polynomial const &);

        mpz_class &operator()(size_t r, size_t c);
        mpz_class const &operator()(size_t r, size_t c) const;

        [[nodiscard]]
        std::vector<vector> to_row_vectors() const;
        [[nodiscard]]
        std::vector<vector> to_col_vectors() const;

        integer_matrix to_smith_normal_form();
        [[nodiscard]]
        std::pair<integer_matrix, integer_matrix> smith_normal_form() const;

        mpz_class to_congruent_diagonal();
        [[nodiscard]]
        std::pair<integer_matrix, mpz_class> congruent_diagonal() const;

        [[nodiscard]]
        std::pair<size_t, size_t> signature() const;

        [[nodiscard]]
        bool is_symmetric() const;

        [[nodiscard]]
        bool is_square() const {
            return row == col;
        }

        [[nodiscard]]
        std::vector<mpz_class> get_diagonal() const;

        [[nodiscard]]
        vector operate_on(vector const &) const;

        integer_matrix &operator+=(integer_matrix const &matrix) { return add_eq(matrix); }

        integer_matrix &operator-=(integer_matrix const &matrix) { return sub_eq(matrix); }

        integer_matrix &operator*=(integer_matrix const &matrix) { return mul_eq(matrix); }

        integer_matrix &operator*=(mpz_class const &m) { return mul_scalar_eq(m); }

        integer_matrix &operator*=(vector const &vec) { return mul_vector_eq(vec); }

        [[nodiscard]]
        integer_matrix
        operator+(integer_matrix const &matrix) const {
            return add(matrix);
        }

        [[nodiscard]]
        integer_matrix
        operator-(integer_matrix const &matrix) const {
            return sub(matrix);
        }

        [[nodiscard]]
        integer_matrix
        operator*(integer_matrix const &matrix) const {
            return mul(matrix);
        }

        [[nodiscard]]
        integer_matrix
        operator*(mpz_class const &m) const {
            return mul_scalar(m);
        }

        friend integer_matrix operator*(mpz_class const &m, integer_matrix const &matrix) {
            return matrix.mul_scalar(m);
        }

        [[nodiscard]]
        integer_matrix
        operator*(vector const &vec) const {
            return mul_vector(vec);
        }

        friend std::ostream &operator<<(std::ostream &os, integer_matrix const &M);
        friend vector solve(integer_matrix, vector, mpz_class const &);
        friend std::vector<vector> kernel(integer_matrix);
        friend std::vector<vector> kernel(integer_matrix, mpz_class const &);
        friend mpz_class determinant(integer_matrix);

    private:
        size_t row, col;
        std::vector<mpz_class> mat;

        mpz_class &get(size_t r, size_t c);
        [[nodiscard]]
        mpz_class const &get(size_t r, size_t c) const;

        void row_exchange(size_t l, size_t r);
        void row_negate(size_t r);
        void row_multiply(size_t, mpz_class const &);
        void row_multiply_and_add(size_t from, size_t to, mpz_class const &mul);
        void row_gcd_operation(size_t piv, size_t l, size_t r);

        void col_exchange_with_companion_matrix(integer_matrix &M, size_t l, size_t r);
        void col_negate_with_companion_matrix(integer_matrix &M, size_t c);
        void col_multiply_and_add_with_companion_matrix(integer_matrix &M, size_t from, size_t to, mpz_class const &mul);
        void col_gcd_operation_with_companion_matrix(integer_matrix &M, size_t piv, size_t l, size_t r);

        void col_exchange(size_t l, size_t r);
        void col_negate(size_t c);
        void col_multiply(size_t, mpz_class const &);
        void col_multiply_and_add(size_t from, size_t to, mpz_class const &mul);
        void col_gcd_operation(size_t piv, size_t l, size_t r);

        integer_matrix &add_eq(integer_matrix const &);
        integer_matrix &sub_eq(integer_matrix const &);
        integer_matrix &mul_eq(integer_matrix const &);
        integer_matrix &mul_scalar_eq(mpz_class const &);
        integer_matrix &mul_vector_eq(vector const &);

        [[nodiscard]]
        integer_matrix add(integer_matrix const &) const;
        [[nodiscard]]
        integer_matrix sub(integer_matrix const &) const;
        [[nodiscard]]
        integer_matrix mul(integer_matrix const &) const;
        [[nodiscard]]
        integer_matrix mul_scalar(mpz_class const &) const;
        [[nodiscard]]
        integer_matrix mul_vector(vector const &) const;
    };
}  // namespace gmp_extensions

#endif
