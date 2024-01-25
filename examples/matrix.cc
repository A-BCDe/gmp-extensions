#include "matrix.h"
#include "number_theoretic.h"

#include <chrono>
#include <iostream>

int main() {
    std::cout << "Hello, world!" << std::endl;
    auto &pg = gmp_extensions::the_prime_generator();
    auto start = std::chrono::system_clock::now();
    while(true) {
        auto const p = pg.next_prime();
        if(p > 200'000'000) break;
    }
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::cout << "\nElapsed time: " << sec.count() << " s\n";

    /*
    gmp_extensions::integer_matrix mat1(2, 2, {
            "1", "1",
            "1", "0"
    });
    gmp_extensions::vector vec1({ 1, 0 });
    std::cout << mat1 * mat1 * mat1 * vec1 << '\n';
    std::cout << "determinant of\n" << mat1 << "\nis " << determinant(mat1) <<
    "\n\n"; gmp_extensions::integer_matrix mat(3, 3, { "0", "1", "2", "0", "1",
    "3", "1", "2", "3"}); std::cout << "determinant of mat = " <<
    determinant(mat).get_str() << '\n'; gmp_extensions::vector v({ "5", "7", "8"
    }); auto res = solve(mat, v, 991); std::cout << res << '\n'; std::cout <<
    determinant(mat) << '\n';*/
}