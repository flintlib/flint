/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <string.h>
#include "ulong_extras.h"
#include "fmpz.h"

char * fmpz_get_str_bsplit_threaded(char * s, const fmpz_t f);

TEST_FUNCTION_START(fmpz_get_str, state)
{
    int i;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t b;
        int base;
        char *str1, *str2;

        fmpz_init(a);
        mpz_init(b);
        fmpz_randtest(a, state, 200);
        base = (int) (n_randint(state, 61) + 2);

        fmpz_get_mpz(b, a);

        str1 = fmpz_get_str(NULL, base, a);
        str2 = mpz_get_str(NULL, base, b);

        if (strcmp(str1, str2))
        {
            flint_printf("FAIL:\n");
            gmp_printf("b = %Zd\n", b);
            flint_printf("base = %d\n", base);
            flint_printf("str1 = %s\n, str2 = %s\n", str1, str2);
            fflush(stdout);
            flint_abort();
        }

        flint_free(str1);
        flint_free(str2);

        fmpz_clear(a);
        mpz_clear(b);
    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t b;
        int base;
        char *str1, *str2;

        flint_set_num_threads(1 + n_randint(state, 5));

        fmpz_init(a);
        mpz_init(b);
        fmpz_randbits(a, state, 1000000 + n_randint(state, 1000000));
        base = 10;

        fmpz_get_mpz(b, a);

        if (n_randint(state, 2))
            str1 = fmpz_get_str(NULL, base, a);
        else
            str1 = fmpz_get_str_bsplit_threaded(NULL, a);

        str2 = mpz_get_str(NULL, base, b);

        if (strcmp(str1, str2))
        {
            flint_printf("FAIL:\n");
            gmp_printf("b = %Zd\n", b);
            flint_printf("base = %d\n", base);
            flint_printf("str1 = %s\n, str2 = %s\n", str1, str2);
            fflush(stdout);
            flint_abort();
        }

        flint_free(str1);
        flint_free(str2);

        fmpz_clear(a);
        mpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
