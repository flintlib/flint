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

void test_exact(int d)
{
    int i, j;
    ulong exp, exact;
    mpz_t a, a2, b, c;
    mpz_init(a);
    mpz_init(a2);
    mpz_init(b);
    mpz_init(c);
    for (i=0; i<100; i++)
    {
        for (j=1; j<100; j++)
        {
            exact = i / j;
            flint_mpz_set_ui(a, d);
            flint_mpz_pow_ui(a, a, i);
            mpz_set(a2, a);
            flint_mpz_set_ui(b, d);
            flint_mpz_pow_ui(b, b, j);
            a->_mp_size = flint_mpn_remove_power_ascending(a->_mp_d, a->_mp_size,
                b->_mp_d, b->_mp_size, &exp);
            flint_mpz_pow_ui(b, b, exact);
            mpz_tdiv_q(c, a2, b);
            if (exp != i/j || mpz_cmp(a, c))
            {
                gmp_printf("%d^%d / %d^%d\n", d, i, d, j);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    mpz_clear(a);
    mpz_clear(a2);
    mpz_clear(b);
    mpz_clear(c);
}

TEST_FUNCTION_START(flint_mpn_remove_power, state)
{
    test_exact(3);
    test_exact(10);
    test_exact(7429);

    TEST_FUNCTION_END(state);
}
