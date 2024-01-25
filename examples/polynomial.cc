#include "polynomial.h"
#include <chrono>
#include <iostream>

int main() {
    std::cout << "Hello, world!" << std::endl;
    // x^6 - 3x^5 + x^4 - 3x^3 - x^2 - 3x + 1
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "1",
    // "-3", "-1", "-3", "1", "-3", "1" });
    // (x + 1)(x^2 + 5x + 3)(x^3 + 2x^2 + 3x + 4)
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "12",
    // "41", "54", "41", "23", "8", "1" });
    // (x - 1)(x - 2)
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "2",
    // "-3", "1" });
    // (x - 1)(x - 2)(x + 3)
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "6",
    // "-7", "0", "1" });
    // (x^2 + 1)(x^2 + 2)
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "2", "0",
    // "3", "0", "1" });
    // x^2 + 2
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "2", "0",
    // "1" });
    // (x^2 + 2)(x^3 + 2)
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "4", "0",
    // "2", "2", "0", "1" });
    // (x^2 + 2)(x^3 + 3)
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "6", "0",
    // "3", "2", "0", "1" });
    // (x^2 + 2)(x^2 - 2)
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "-4",
    // "0", "0", "0", "1" });
    // (3x + 1)(3x + 2)(3x + 3)(3x + 4)
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "24",
    // "150", "315", "270", "81" }); (2x^2 + 3x + 4)(7x^2 - x + 53)
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "212",
    // "155", "131", "19", "14" });
    // (x^18 - 1)
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "-1",
    // "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0",
    // "0", "0", "1" });
    // (x^2 + 1)^3
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "1", "0",
    // "3", "0", "3", "0", "1" });
    // (x^2 + 1)^2 (x^2 - x + 1)^3
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "1",
    // "-3", "8", "-13", "19", "-20", "19", "-13", "8", "-3", "1" });
    // (x - 1) (x^2 + 1)^2 (x^3 - 1)^3
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{ "1",
    // "-1", "2", "-5", "4", "-7", "9", "-6", "9", "-7", "4", "-5", "2", "-1", "1"
    // });
    // irreducible polynomial of degree 9
    // gmp_extensions::integer_polynomial poly(std::vector<std::string>{
    // "-28314171263536", "0", "0", "0", "0", "12114901", "0", "-1939601", "0",
    // "21354101" });
    // (x^2 + x - 1)(x^3 - x^2 + 2x - 3)
    gmp_extensions::integer_polynomial poly(std::vector<std::string> {"3", "-5", "0", "0", "0", "1"});

    /* A code for factorizing modulo p
     * Should not work for the test case (3x + 1)(3x + 2)(3x + 3)(3x + 4), as the
     leading coefficient is zero. mpz_class const p = 3; std::cout << poly <<
     '\n'; auto factors = poly.factorize_squarefree(p); std::cout << "Irreducible
     factors of " << poly << " mod " << p << " are:\n"; for(auto const &v :
     factors) { std::cout << v << '\n';
            }
    */
    auto start = std::chrono::system_clock::now();
    auto factors = poly.factorize();
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::cout << "Irreducible factors of " << poly << " are:\n";
    for(auto const &v: factors) { std::cout << v << '\n'; }
    std::cout << "\nElapsed time: " << sec.count() << " s\n";
}