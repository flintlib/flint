/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "gmpcompat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_sizeinbase, state)
{
    mp_limb_t n;
    int base, size1, size2;
    slong rep;
    mpz_t t;
    char * str;

    mpz_init(t);
    str = flint_malloc((FLINT_BITS + 1) * sizeof(char));

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        n = n_randtest(state);
        base = 2 + n_randint(state, 34);

        size1 = n_sizeinbase(n, base);

        flint_mpz_set_ui(t, n);

        mpz_get_str(str, base, t);
        size2 = strlen(str);

        if (size1 != size2)
        {
            flint_printf("FAIL: n = %wu, base = %d\n", n, base);
            flint_printf("n_sizeinbase: %d, strlen: %d\n", size1, size2);
            fflush(stdout);
            flint_abort();
        }
    }

    flint_free(str);
    mpz_clear(t);

    TEST_FUNCTION_END(state);
}
