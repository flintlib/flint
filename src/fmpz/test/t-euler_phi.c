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

TEST_FUNCTION_START(fmpz_euler_phi, state)
{
    slong i;
    ulong n;
    fmpz_t x, y, z;

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(z);

    for (i = 0; i < 100; i++)
    {
        fmpz_set_ui(x, i);
        fmpz_euler_phi(y, x);
        fmpz_euler_phi(x, x);
        fmpz_set_ui(z, n_euler_phi(i));
        if (!fmpz_equal(x, y) || !fmpz_equal(x, z))
        {
            flint_printf("FAIL: %wd\n", i);
            fflush(stdout);
            flint_abort();
        }
    }

    /* Aliasing test */
    for (i = 0; i < 1000; i++)
    {
        fmpz_randtest(x, state, FLINT_BITS);
        fmpz_randtest(y, state, 5);
        fmpz_pow_ui(y, y, n_randtest(state) % 100);
        fmpz_mul(x, x, y);
        fmpz_set(z, x);
        fmpz_euler_phi(y, x);
        fmpz_euler_phi(x, x);
        if (!fmpz_equal(x, y))
        {
            flint_printf("FAIL: ");
            fmpz_print(z);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    /* Power of a single prime, phi(p^n) = (p-1) * p^(n-1) */
    for (i = 0; i < 100; i++)
    {
        n = (n_randtest(state) % 100) + 1;
        fmpz_set_ui(x, n_nth_prime(i+1));
        fmpz_pow_ui(x, x, n);
        fmpz_euler_phi(x, x);
        fmpz_set_ui(y, n_nth_prime(i+1));
        fmpz_pow_ui(y, y, n-1);
        fmpz_mul_ui(y, y, n_nth_prime(i+1)-1);
        if (!fmpz_equal(x, y) || !_fmpz_is_canonical(x))
        {
            flint_printf("FAIL: %wu ^ %wu\n", n_nth_prime(i+1), n);
        }
    }

    /* Something nontrivial */
    fmpz_set_str(x, "10426024348053113487152988625265848110501553295256578345594388516660144", 10);
    fmpz_set_str(y, "2265085829098571747262267425315881590169106756213617459200000000000000", 10);
    fmpz_euler_phi(x, x);
    if (!fmpz_equal(x, y) || !_fmpz_is_canonical(x))
    {
        flint_printf("FAIL: special test value\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(z);

    TEST_FUNCTION_END(state);
}
