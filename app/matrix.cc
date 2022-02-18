#include <matrix.h>

#include <iostream>

int main() {
	std::cout << "Hello, world!" << std::endl;
	project::integer_matrix mat(3, 3, {
		"0", "1", "2",
		"0", "1", "3",
		"1", "2", "3"});
	std::cout << "determinant of mat = " << determinant(mat).get_str() << '\n';
	project::vector v({ "5", "7", "8" });
	auto res = solve(mat, v, 991);
	std::cout << res << '\n';
}