/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Opus 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "ulong_extras.h"

#if FLINT_BITS == 64

TEST_FUNCTION_START(n_ll_factor_ecm, state)
{
    slong i;
    ulong count = 0, trials = 0;

    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_t p, q, n, f, t;
        ulong nhi, nlo, fac[2];
        flint_bitcnt_t pbits, qbits;

        fmpz_init(p);
        fmpz_init(q);
        fmpz_init(n);
        fmpz_init(f);
        fmpz_init(t);

        /* a factor ECM has a fair chance of finding, and a large cofactor */
        pbits = 20 + n_randint(state, 18);
        qbits = 128 - pbits - n_randint(state, 20);
        if (pbits + qbits < 65)
            qbits = 65 - pbits;

        fmpz_randprime(p, state, pbits, 0);
        fmpz_randprime(q, state, qbits, 0);
        fmpz_mul(n, p, q);

        if (fmpz_bits(n) > 64 && fmpz_bits(n) <= 128 && fmpz_is_odd(n))
        {
            trials++;
            fmpz_get_uiui(&nhi, &nlo, n);

            if (n_ll_factor_ecm(fac, nhi, nlo, 4, 1000, 100000, state))
            {
                fmpz_set_uiui(f, fac[1], fac[0]);

                if (fmpz_is_one(f) || fmpz_equal(f, n))
                    TEST_FUNCTION_FAIL("trivial factor returned\n"
                                       "n = %{fmpz}\nf = %{fmpz}\n", n, f);

                fmpz_mod(t, n, f);
                if (!fmpz_is_zero(t))
                    TEST_FUNCTION_FAIL("factor does not divide n\n"
                                       "n = %{fmpz}\nf = %{fmpz}\n", n, f);
                count++;
            }
        }

        fmpz_clear(p);
        fmpz_clear(q);
        fmpz_clear(n);
        fmpz_clear(f);
        fmpz_clear(t);
    }

    if (trials > 20 && count < trials / 5)
        TEST_FUNCTION_FAIL("only %wu of %wu numbers factored\n", count, trials);

    TEST_FUNCTION_END(state);
}

#else

TEST_FUNCTION_START(n_ll_factor_ecm, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}

#endif
