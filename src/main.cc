#include <gmpxx.h>

#include <iostream>
#include <integer_polynomial.h>

int main() {
	std::cout << "Hello, world!" << std::endl;
	mpz_class a("123456789123456789123456789123456789"), b("987654321987654321897645321987654312");
	std::cout << a + b << std::endl;
}