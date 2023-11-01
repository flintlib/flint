/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "ulong_extras.h"

#ifndef check
#define check check
/* Defined in t-abs_fits_ui.c, t-fits_si.c and t-moebius_mu.c */
static void check(fmpz_t input, int output, int expected)
{
    if (output != expected)
    {
        printf("FAIL:\n\n"
               "input = "), fmpz_print(input), printf("\n"
               "output = %d, expected = %d\n", output, expected);
        fflush(stdout);
        flint_abort();
    }
}
#endif

TEST_FUNCTION_START(fmpz_moebius_mu, state)
{
    fmpz_t x;
    ulong p;
    slong i, j, k, l;

    fmpz_init(x);

    for (i = -1000; i < 1000; i++)
    {
        fmpz_set_si(x, i);
        check(x, fmpz_moebius_mu(x), n_moebius_mu(FLINT_ABS(i)));
    }

    for (i = 0; i < 1000; i++)
    {
        fmpz_set_ui(x, 1);
        /* Product of some primes */
        k = n_randtest(state) % 10;
        l = n_randtest(state) % 10;
        for (j = 0; j < k; j++)
        {
            l += (n_randtest(state) % 10) + 1;
            fmpz_mul_ui(x, x, n_nth_prime(l+1));
        }

        check(x, fmpz_moebius_mu(x), (k % 2 ? -1 : 1));
        fmpz_neg(x, x);

        check(x, fmpz_moebius_mu(x), (k % 2 ? -1 : 1));
        fmpz_abs(x, x);

        /* No longer square-free */
        p = n_nth_prime(n_randtest(state) % 100 + 1);
        fmpz_mul_ui(x, x, p*p);
        check(x, fmpz_moebius_mu(x), 0);
    }

    fmpz_clear(x);

    TEST_FUNCTION_END(state);
}
