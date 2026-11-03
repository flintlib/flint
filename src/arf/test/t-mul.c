/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

/* todo: try some other code here also, e.g. fmpz roundtrip */
int
arf_mul_naive(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
{
    arf_mul(z, x, y, ARF_PREC_EXACT, rnd);
    return arf_set_round(z, z, prec, rnd);
}

/* Random value with few nonzero bits, so that products have long runs of
   zero or one bits around the rounding point (exercising the fallback of
   the truncated multiplication to an exact product). */
static void
arf_randtest_sparse(arf_t x, flint_rand_t state, slong bits, slong mag_bits)
{
    fmpz_t t, u;
    slong i, terms;

    fmpz_init(t);
    fmpz_init(u);

    terms = 1 + n_randint(state, 4);

    for (i = 0; i < terms; i++)
    {
        fmpz_one(u);
        fmpz_mul_2exp(u, u, n_randint(state, bits));
        if (n_randint(state, 2))
            fmpz_add(t, t, u);
        else
            fmpz_sub(t, t, u);
    }

    arf_set_fmpz(x, t);
    if (mag_bits != 0)
        arf_mul_2exp_si(x, x, (slong) n_randint(state, 2 * mag_bits) - mag_bits);

    fmpz_clear(t);
    fmpz_clear(u);
}

/* Checks the relaxed rounding contract: |z - x y| <= tol ulp where
   tol = 100 (ARF_RND_FAST) or 51 (ARF_RND_ACCURATE) hundredths of an ulp
   of the exact product at prec bits; and z must be exact when the
   returned inexact flag is 0. */
static void
check_relaxed(const arf_t z, const arf_t x, const arf_t y, slong prec,
    arf_rnd_t rnd, int inexact, const char * label)
{
    arf_t v, err, tol;
    fmpz_t e;
    int ok;

    arf_init(v);
    arf_init(err);
    arf_init(tol);
    fmpz_init(e);

    arf_mul(v, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);

    if (arf_is_special(v) || prec == ARF_PREC_EXACT)
    {
        ok = arf_equal(z, v);
    }
    else
    {
        arf_sub(err, z, v, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_abs(err, err);
        arf_mul_ui(err, err, 100, ARF_PREC_EXACT, ARF_RND_DOWN);

        fmpz_sub_ui(e, ARF_EXPREF(v), prec);
        arf_set_ui(tol, (rnd == ARF_RND_ACCURATE) ? 51 : 100);
        arf_mul_2exp_fmpz(tol, tol, e);

        ok = (arf_cmp(err, tol) <= 0);
        ok = ok && (inexact || arf_equal(z, v));
    }

    if (!ok)
    {
        flint_printf("FAIL (%s)!\n", label);
        flint_printf("prec = %wd, rnd = %d, inexact = %d\n\n", prec, rnd, inexact);
        flint_printf("x = "); arf_print(x); flint_printf("\n\n");
        flint_printf("y = "); arf_print(y); flint_printf("\n\n");
        flint_printf("z = "); arf_print(z); flint_printf("\n\n");
        flint_printf("v = "); arf_print(v); flint_printf("\n\n");
        flint_abort();
    }

    arf_clear(v);
    arf_clear(err);
    arf_clear(tol);
    fmpz_clear(e);
}

TEST_FUNCTION_START(arf_mul, state)
{
    slong iter, iter2;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z, v;
        slong prec, r1, r2;
        arf_rnd_t rnd;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            slong bits;

            /* mostly small, occasionally very large (FFT range) operands */
            switch (n_randint(state, 20))
            {
                case 0: bits = 50000; break;
                case 1: bits = 6000; break;
                default: bits = 2000; break;
            }

            if (n_randint(state, 4) == 0)
                arf_randtest_sparse(x, state, bits, 100);
            else
                arf_randtest_special(x, state, bits, 100);

            if (n_randint(state, 4) == 0)
                arf_randtest_sparse(y, state, bits, 100);
            else
                arf_randtest_special(y, state, bits, 100);

            prec = 2 + n_randint(state, bits);

            if (n_randint(state, 50) == 0)
                prec = ARF_PREC_EXACT;

            switch (n_randint(state, 5))
            {
                case 0:  rnd = ARF_RND_DOWN; break;
                case 1:  rnd = ARF_RND_UP; break;
                case 2:  rnd = ARF_RND_FLOOR; break;
                case 3:  rnd = ARF_RND_CEIL; break;
                default: rnd = ARF_RND_NEAR; break;
            }

            switch (n_randint(state, 7))
            {
            case 5:
                rnd = n_randint(state, 2) ? ARF_RND_FAST : ARF_RND_ACCURATE;
                if (n_randint(state, 2))
                {
                    r1 = arf_mul(z, x, y, prec, rnd);
                    check_relaxed(z, x, y, prec, rnd, r1, "relaxed");
                }
                else
                {
                    r1 = arf_mul(z, x, x, prec, rnd);
                    check_relaxed(z, x, x, prec, rnd, r1, "relaxed squaring");
                }
                break;

            case 6:
                rnd = n_randint(state, 2) ? ARF_RND_FAST : ARF_RND_ACCURATE;
                arf_set(v, x);
                r1 = arf_mul(x, x, y, prec, rnd);
                check_relaxed(x, v, y, prec, rnd, r1, "relaxed aliasing");
                break;

            case 0:
                r1 = arf_mul(z, x, y, prec, rnd);
                r2 = arf_mul_naive(v, x, y, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            case 1:
                r1 = arf_mul(z, x, x, prec, rnd);
                r2 = arf_mul_naive(v, x, x, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 1)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            case 2:
                r2 = arf_mul_naive(v, x, x, prec, rnd);
                r1 = arf_mul(x, x, x, prec, rnd);
                if (!arf_equal(v, x) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 2)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            case 3:
                r2 = arf_mul_naive(v, x, y, prec, rnd);
                r1 = arf_mul(x, x, y, prec, rnd);
                if (!arf_equal(x, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 3)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            default:
                r2 = arf_mul_naive(v, x, y, prec, rnd);
                r1 = arf_mul(x, y, x, prec, rnd);
                if (!arf_equal(x, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 4)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;
            }
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(v);
    }

    TEST_FUNCTION_END(state);
}
