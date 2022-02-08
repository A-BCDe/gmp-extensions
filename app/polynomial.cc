#include <gmpxx.h>

#include <iostream>
#include <prime_generator.h>

int main() {
	std::cout << "Hello, world!" << std::endl;
	project::prime_generator pg;
	while(true) {
		std::cout << pg.next_prime().get_str() << '\n';
	}
}
