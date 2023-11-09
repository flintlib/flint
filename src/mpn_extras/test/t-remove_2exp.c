/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_remove_2exp, state)
{
    int zero, nonzero;
    flint_bitcnt_t check;
    mpz_t a;
    mpz_t b;

    mpz_init(a);
    mpz_init(b);

    for (zero=0; zero<300; zero++)
    {
        for (nonzero=0; nonzero<300; nonzero++)
        {
            flint_mpz_set_ui(a, 1);
            mpz_setbit(a, nonzero);
            mpz_set(b, a);
            mpz_mul_2exp(a, a, zero);
            a->_mp_size = flint_mpn_remove_2exp(a->_mp_d, a->_mp_size, &check);
            if (check != zero || mpz_cmp(a,b))
            {
                gmp_printf("%d %d \n", zero, nonzero);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    mpz_clear(a);
    mpz_clear(b);

    TEST_FUNCTION_END(state);
}
