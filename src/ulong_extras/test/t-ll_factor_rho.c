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

TEST_FUNCTION_START(n_ll_factor_rho, state)
{
    slong i;
    ulong count = 0, trials = 0;

    for (i = 0; i < 400 * flint_test_multiplier(); i++)
    {
        fmpz_t p, q, n, f, t;
        ulong nhi, nlo, fac[2];
        flint_bitcnt_t pbits, qbits;

        fmpz_init(p);
        fmpz_init(q);
        fmpz_init(n);
        fmpz_init(f);
        fmpz_init(t);

        /* p small enough that rho has a fair chance, q making n two words */
        pbits = 10 + n_randint(state, 22);
        qbits = 65 - pbits + n_randint(state, 128 - 65);
        if (pbits + qbits > 128)
            qbits = 128 - pbits;

        fmpz_randprime(p, state, pbits, 0);
        fmpz_randprime(q, state, qbits, 0);
        fmpz_mul(n, p, q);

        if (fmpz_bits(n) > 64 && fmpz_bits(n) <= 128 && fmpz_is_odd(n))
        {
            trials++;
            fmpz_get_uiui(&nhi, &nlo, n);

            if (n_ll_factor_rho(fac, nhi, nlo, 3, 4096))
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

    /* arguments outside the supported range are rejected */
    {
        ulong fac[2];

        if (n_ll_factor_rho(fac, 0, UWORD(101), 1, 1024))
            TEST_FUNCTION_FAIL("accepted a one limb argument\n");

        if (n_ll_factor_rho(fac, 1, UWORD(0), 1, 1024))
            TEST_FUNCTION_FAIL("accepted an even argument\n");
    }

    /*
       A small modulus with a large budget makes the batched product run into
       every factor at once, so the gcd comes back equal to n and the attempt
       has to be restarted with a different constant.
    */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_t p, q, n, f, t;
        ulong nhi, nlo, fac[2];

        fmpz_init(p);
        fmpz_init(q);
        fmpz_init(n);
        fmpz_init(f);
        fmpz_init(t);

        fmpz_randprime(p, state, 33, 0);
        fmpz_randprime(q, state, 33, 0);
        fmpz_mul(n, p, q);

        if (fmpz_bits(n) > 64 && fmpz_is_odd(n))
        {
            fmpz_get_uiui(&nhi, &nlo, n);

            if (n_ll_factor_rho(fac, nhi, nlo, 4, 65536))
            {
                fmpz_set_uiui(f, fac[1], fac[0]);
                fmpz_mod(t, n, f);

                if (!fmpz_is_zero(t) || fmpz_is_one(f) || fmpz_equal(f, n))
                    TEST_FUNCTION_FAIL("bad factor with a large budget\n"
                                       "n = %{fmpz}\nf = %{fmpz}\n", n, f);
            }
        }

        fmpz_clear(p);
        fmpz_clear(q);
        fmpz_clear(n);
        fmpz_clear(f);
        fmpz_clear(t);
    }

    /*
       Three factors, two of which multiply to more than one limb, so the
       batched gcd can return a factor that needs both limbs.
    */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p, q, r, n, f, t;
        ulong nhi, nlo, fac[2];

        fmpz_init(p);
        fmpz_init(q);
        fmpz_init(r);
        fmpz_init(n);
        fmpz_init(f);
        fmpz_init(t);

        fmpz_randprime(p, state, 34, 0);
        fmpz_randprime(q, state, 34, 0);
        fmpz_randprime(r, state, 20, 0);
        fmpz_mul(n, p, q);
        fmpz_mul(n, n, r);

        if (fmpz_bits(n) > 64 && fmpz_bits(n) <= 128 && fmpz_is_odd(n))
        {
            fmpz_get_uiui(&nhi, &nlo, n);

            if (n_ll_factor_rho(fac, nhi, nlo, 3, 8192))
            {
                fmpz_set_uiui(f, fac[1], fac[0]);
                fmpz_mod(t, n, f);

                if (!fmpz_is_zero(t) || fmpz_is_one(f) || fmpz_equal(f, n))
                    TEST_FUNCTION_FAIL("bad factor for a three factor input\n"
                                       "n = %{fmpz}\nf = %{fmpz}\n", n, f);
            }
        }

        fmpz_clear(p);
        fmpz_clear(q);
        fmpz_clear(r);
        fmpz_clear(n);
        fmpz_clear(f);
        fmpz_clear(t);
    }

    /* the map should succeed on a decent fraction of these */
    if (trials > 100 && count < trials / 4)
        TEST_FUNCTION_FAIL("only %wu of %wu numbers factored\n", count, trials);

    TEST_FUNCTION_END(state);
}

#else

TEST_FUNCTION_START(n_ll_factor_rho, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}

#endif
