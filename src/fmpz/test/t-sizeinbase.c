/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <string.h>
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_sizeinbase, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t b;
        char * str;
        int base;
        size_t r1, exact;

        fmpz_init(a);
        mpz_init(b);
        fmpz_randtest(a, state, 200);
        base = (int) (n_randint(state, 61) + 2);

        fmpz_get_mpz(b, a);

        r1 = fmpz_sizeinbase(a, base);

        /* exact size = length of |b| written in base b */
        str = mpz_get_str(NULL, base, b);
        exact = strlen(str) - (str[0] == '-');
        flint_free(str);

        result = (r1 == exact);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("b = %Zd\n", b);
            flint_printf("base = %d\n", base);
            flint_printf("r1 = %wu, exact = %wu\n",
                         (ulong) r1, (ulong) exact);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        mpz_clear(b);
    }

    /* 10^19 - 1 has 19 decimal digits, not 20. */
    {
        fmpz_t a;
        size_t r;

        fmpz_init(a);
        fmpz_set_ui(a, 10);
        fmpz_pow_ui(a, a, 19);
        fmpz_sub_ui(a, a, 1);

        r = fmpz_sizeinbase(a, 10);
        if (r != 19)
        {
            flint_printf("FAIL (10^19 - 1):\n");
            flint_printf("fmpz_sizeinbase(10^19 - 1, 10) = %wu, expected 19\n",
                         (ulong) r);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
    }

    TEST_FUNCTION_END(state);
}
