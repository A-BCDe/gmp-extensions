#ifndef INTEGER_EXTENSIONS_INTEGER_RANDOM_H
#define INTEGER_EXTENSIONS_INTEGER_RANDOM_H

#include <gmp.h>
#include <gmpxx.h>

namespace gmp_extensions {

    class integer_random {
    public:
        integer_random() { gmp_randinit_default(randstate); }

        ~integer_random() { gmp_randclear(randstate); }

        integer_random(integer_random const &) = delete;
        integer_random(integer_random &&) = default;
        integer_random &operator=(integer_random const &) = delete;
        integer_random &operator=(integer_random &&) = default;

        mpz_class operator()(mpz_class const &lo, mpz_class const &hi);

    private:
        gmp_randstate_t randstate;
    };

}  // namespace gmp_extensions

#endif
