/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "arb.h"
#include "fixed.h"

/* fixed_exp_diophantine agrees with arb within its budget for x
   anywhere in [0, 1), across the sizes and weight budgets that
   steer the reduction (few rows and small p, q at low weight;
   the full relation table and the negative-residual lifting at
   high weight), and the cached prime logarithms are exact floors
   at every read precision */

static double
_ulps_mp(nn_srcptr y, slong ylen, const arb_t exact, slong n)
{
    arb_t va, d;
    fmpz_t f;
    mag_t mb;
    double u;

    arb_init(va); arb_init(d); fmpz_init(f); mag_init(mb);
    fmpz_set_ui_array(f, y, ylen);
    arb_set_fmpz(va, f);
    arb_mul_2exp_si(va, va, -FLINT_BITS * n);
    arb_sub(d, va, exact, ARF_PREC_EXACT);
    arb_get_mag(mb, d);
    mag_mul_2exp_si(mb, mb, FLINT_BITS * n);
    u = mag_get_d(mb);
    arb_clear(va); arb_clear(d); fmpz_clear(f); mag_clear(mb);
    return u;
}

TEST_FUNCTION_START(fixed_exp_diophantine, state)
{
    slong iter;

    /* the cache: entries at any read precision are exact floors */
    for (iter = 0; iter < 20 + 5 * flint_test_multiplier(); iter++)
    {
        slong nv = 1 + n_randint(state, (iter % 5 == 0) ? 300 : 30);
        slong num = 2 + n_randint(state, (iter % 3 == 0) ? 30 : 13);
        slong j = n_randint(state, num);
        arb_t l;
        fmpz_t f, g;

        if (iter % 7 == 0)
            _fixed_log_primes_clear();

        _fixed_log_primes_ensure(num, nv);

        arb_init(l); fmpz_init(f); fmpz_init(g);
        arb_log_ui(l, n_nth_prime(j + 1), FLINT_BITS * nv + 200);
        arb_mul_2exp_si(l, l, FLINT_BITS * nv);
        arb_floor(l, l, FLINT_BITS * nv + 200);
        if (!arb_get_unique_fmpz(f, l))
            TEST_FUNCTION_FAIL("floor not determined: j = %wd, nv = %wd\n",
                j, nv);
        fmpz_set_ui_array(g, _fixed_log_primes_entry(j, nv), nv + 1);
        if (!fmpz_equal(f, g))
            TEST_FUNCTION_FAIL("cache: j = %wd, nv = %wd\n%{fmpz}\n%{fmpz}\n",
                j, nv, f, g);
        arb_clear(l); fmpz_clear(f); fmpz_clear(g);
    }

    for (iter = 0; iter < 60 + 60 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, (iter % 17 == 0) ? 800 : 60);
        slong zz, k;
        slong prec = FLINT_BITS * n + 128;
        nn_ptr x, y;
        arb_t xa, e;
        fmpz_t f;
        double u, w;

        x = flint_malloc(n * sizeof(ulong));
        y = flint_malloc((n + 1) * sizeof(ulong));
        arb_init(xa); arb_init(e); fmpz_init(f);

        flint_mpn_urandomb(x, state, FLINT_BITS * n);
        if (iter % 11 == 4)
            flint_mpn_store(x, n, ~UWORD(0));
        if (iter % 13 == 6)
        {
            /* few significant bits */
            flint_mpn_zero(x, n);
            x[n - 1] = n_randtest(state);
        }

        /* sweep the leading zeros: a tiny x makes the reduction
           choose nothing beyond the log 2 step and the residual
           very deep */
        switch (iter % 8)
        {
            case 0: zz = 0; break;
            case 1: zz = 1 + n_randint(state, 15); break;
            case 2: zz = 16 + n_randint(state, 16); break;
            case 3: zz = 32 + n_randint(state, 33); break;
            case 4: zz = 65 + n_randint(state, 200); break;
            case 5: zz = n_randint(state, FLINT_BITS * n); break;
            case 6: zz = FLINT_BITS * n - 1
                - n_randint(state, FLINT_MIN(FLINT_BITS * n, 8));
                break;
            default: zz = 0; break;
        }
        {
            slong q = FLINT_MIN(zz / FLINT_BITS, n), i;
            int b = (int) (zz % FLINT_BITS);
            for (i = 0; i < q; i++)
                x[n - 1 - i] = 0;
            if (b && q < n)
                x[n - 1 - q] >>= b;
        }
        if (iter % 23 == 7)
            flint_mpn_zero(x, n);

        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
        arb_exp(e, xa, prec);

        fixed_exp_diophantine(y, x, n);
        u = _ulps_mp(y, n + 1, e, n);
        if (u > (double) FIXED_EXP_DIOPHANTINE_MAX_ERR)
            TEST_FUNCTION_FAIL("n = %wd, zz = %wd, ulp = %f\n", n, zz, u);

        /* other weight budgets, including a tiny one (the descent
           stops early: large residual, the notab fallback) and an
           effectively unbounded one (every row of the table), and
           other prime counts (generated tables) */
        for (k = 0; k < 3; k++)
        {
            slong num;

            switch (k)
            {
                case 0: w = FLINT_BITS * n * 0.1; break;
                case 1: w = FLINT_BITS * n * 4.0; break;
                default: w = 1e30; break;
            }
            if (n_randint(state, 4) == 0)
                w = FLINT_BITS * n * (n_randint(state, 100) / 10.0);

            num = 13;
            if (n_randint(state, 3) == 0)
                num = 2 + n_randint(state, 23);

            _fixed_exp_diophantine_tune(y, x, n, num, w);
            u = _ulps_mp(y, n + 1, e, n);
            if (u > (double) FIXED_EXP_DIOPHANTINE_MAX_ERR)
                TEST_FUNCTION_FAIL("n = %wd, zz = %wd, num = %wd, w = %f, "
                    "ulp = %f\n", n, zz, num, w, u);
        }

        arb_clear(xa); arb_clear(e); fmpz_clear(f);
        flint_free(x); flint_free(y);
    }

    /* one deep case: the bit-burst regime inside fixed_exp_reduced
       and prime powers of thousands of bits */
    {
        slong n = 3000;
        slong prec = FLINT_BITS * n + 128;
        nn_ptr x, y;
        arb_t xa, e;
        fmpz_t f;
        double u;

        x = flint_malloc(n * sizeof(ulong));
        y = flint_malloc((n + 1) * sizeof(ulong));
        arb_init(xa); arb_init(e); fmpz_init(f);

        flint_mpn_urandomb(x, state, FLINT_BITS * n);

        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);

        fixed_exp_diophantine(y, x, n);
        arb_exp(e, xa, prec);
        u = _ulps_mp(y, n + 1, e, n);
        if (u > (double) FIXED_EXP_DIOPHANTINE_MAX_ERR)
            TEST_FUNCTION_FAIL("deep: ulp = %f\n", u);

        arb_clear(xa); arb_clear(e); fmpz_clear(f);
        flint_free(x); flint_free(y);
    }

    TEST_FUNCTION_END(state);
}
