#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

#include <matrix.h>

template<typename INT>
INT skip_whitespaces(std::string const &s, INT idx) {
    while(idx < s.size() && (s[idx] == 127 || s[idx] <= ' ')) idx++;
    return idx;
}

template<typename INT>
INT skip_noninteger(std::string const &s, INT idx) {
    while (idx < s.size() && !(s[idx] == '-' || (s[idx] >= '0' && s[idx] <= '9'))) {
        idx++;
    }
    return idx;
}

template<typename INT>
INT skip_integer(std::string const &s, INT idx) {
    while (idx < s.size() && (s[idx] == '-' || (s[idx] >= '0' && s[idx] <= '9'))) {
        idx++;
    }
    return idx;
}

template<typename INT>
INT skip_tuple(std::string const &s, INT idx) {
    if(idx >= s.size() || s[idx] != '(') return idx;
    auto const start = idx;
    while(idx < s.size() && s[idx] != ')') idx++;
    if(idx == s.size()) return start;
    return idx + 1;
}

int main() {
    while(true) {
        std::vector<std::string> input;
        std::string string_input;
        size_t index = 0;
        size_t start, end;
        size_t matrix_size;

        // Signature modifier
        std::getline(std::cin, string_input);
        index = skip_whitespaces(string_input, index);
        start = index;
        index = skip_integer(string_input, index);
        end = index;
        if(start == end) break;
        mpz_class signature_modified(std::stoi(string_input.substr(start, end)));

        // Matrix
        matrix_size = 0;
        for(size_t i = 0; i < matrix_size || !matrix_size; i++) {
            std::getline(std::cin, string_input);
            index = 0;
            while(true) {
                index = skip_whitespaces(string_input, index);
                start = index;
                index = skip_integer(string_input, index);
                end = index;
                if(start == end) break;
                input.push_back(string_input.substr(start, end - start));
            }
            if(!matrix_size) matrix_size = (int)input.size();
        }

        matrix::integer_matrix M(matrix_size, matrix_size, input);
        std::cout << M << '\n';

        auto cd = M.congruent_diagonal();
        auto const congruent_diagonal = cd.first;
        auto const congruent_gcd = cd.second;
        auto snf = M.smith_normal_form();
        auto SNF = snf.first; auto Y = snf.second;
        mpz_class determinant(1), signature(0);
        assert(congruent_gcd > 0);

        for(size_t i = 0; i < matrix_size; i++) {
            determinant *= SNF(i, i);
            signature += (((congruent_diagonal(i, i) > 0) << 1) - 1);
        }

        // Input last line
        std::getline(std::cin, string_input);
        index = skip_whitespaces(string_input, (size_t)0);
        index = skip_integer(string_input, index);
        index = skip_whitespaces(string_input, index);
        index = skip_integer(string_input, index);
        std::string prefix = string_input.substr(0, index);
        std::cout << prefix << " :";

        for(size_t i = 0; i < matrix_size; i++) {
            if(SNF(i, i) == 1) continue;
            std::cout << '\t';
            std::cout << SNF(i, i) << '(' << Y(i, 0);
            for(size_t j = 1; j < matrix_size; j++) {
                std::cout << ',' << Y(i, j);
            }
            std::cout << ')';
        }
        std::cout << '\n';

        std::cout << string_input << "\t\t" << abs(determinant) << "\t\t" << signature_modified - signature << '\n';
    }
}
