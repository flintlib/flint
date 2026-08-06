/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "fixed.h"

/* fixed_exp_notab and fixed_sin_cos_notab agree with arb within the
   documented budgets for x anywhere in [0, 1): the leading-zero
   count of x steers the halving depth, so it is swept alongside the
   size, and the forced-depth workers run over every depth the tuned
   tables can select (plus the depth-16 minimum) so that no r is
   covered only where the tables happen to choose it */

static double
_ulps(nn_srcptr y, slong ylen, const arb_t exact, slong n)
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

TEST_FUNCTION_START(fixed_exp_notab, state)
{
    slong iter;

    for (iter = 0; iter < 30 + 40 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, (iter % 13 == 0) ? 800 : 40);
        slong zz;
        slong prec = FLINT_BITS * n + 128;
        nn_ptr x, y, ys, yc;
        arb_t xa, e;
        fmpz_t f;
        double u;
        int rr;

        x = flint_malloc(n * sizeof(ulong));
        y = flint_malloc((n + 1) * sizeof(ulong));
        ys = flint_malloc((n + 1) * sizeof(ulong));
        yc = flint_malloc((n + 1) * sizeof(ulong));
        arb_init(xa); arb_init(e); fmpz_init(f);

        flint_mpn_urandomb(x, state, FLINT_BITS * n);
        if (iter % 11 == 4)
            flint_mpn_store(x, n, ~UWORD(0));

        /* sweep the leading zeros across the halving boundaries */
        switch (iter % 8)
        {
            case 0: zz = 0; break;
            case 1: zz = 1 + n_randint(state, 15); break;
            case 2: zz = 16 + n_randint(state, 16); break;
            case 3: zz = 32 + n_randint(state, 33); break;
            case 4: zz = 65 + n_randint(state, 64); break;
            case 5: zz = n_randint(state, FLINT_BITS * n); break;
            case 6: zz = FLINT_BITS * n - 1
                - n_randint(state, FLINT_MIN(FLINT_BITS * n, 8));
                break;
            default: zz = FLINT_BITS * n; break;   /* x = 0 */
        }
        {
            slong q = FLINT_MIN(zz / FLINT_BITS, n), i;
            int b = (int) (zz % FLINT_BITS);
            for (i = 0; i < q; i++)
                x[n - 1 - i] = 0;
            if (b && q < n)
                x[n - 1 - q] >>= b;
        }

        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);

        fixed_exp_notab(y, x, n);
        arb_exp(e, xa, prec);
        u = _ulps(y, n + 1, e, n);
        if (u > (double) FIXED_EXP_NOTAB_MAX_ERR)
            TEST_FUNCTION_FAIL("exp: n = %wd, zz = %wd, "
                "ulp = %f\n", n, zz, u);

        fixed_sin_cos_notab(ys, yc, x, n);
        arb_sin(e, xa, prec);
        u = _ulps(ys, n + 1, e, n);
        if (u > (double) FIXED_SIN_COS_NOTAB_MAX_ERR)
            TEST_FUNCTION_FAIL("sin: n = %wd, zz = %wd, "
                "ulp = %f\n", n, zz, u);
        arb_cos(e, xa, prec);
        u = _ulps(yc, n + 1, e, n);
        if (u > (double) FIXED_SIN_COS_NOTAB_MAX_ERR)
            TEST_FUNCTION_FAIL("cos: n = %wd, zz = %wd, "
                "ulp = %f\n", n, zz, u);

        /* forced depths: every table value plus the minimum */
        for (rr = 0; rr < 3; rr++)
        {
            static const int rs[3] = { 16, 24, 32 };

            _fixed_exp_notab_r(y, x, n, rs[rr]);
            arb_exp(e, xa, prec);
            u = _ulps(y, n + 1, e, n);
            if (u > (double) FIXED_EXP_NOTAB_MAX_ERR)
                TEST_FUNCTION_FAIL("exp r = %d: n = %wd, "
                    "zz = %wd, ulp = %f\n", rs[rr], n, zz, u);

            _fixed_sin_cos_notab_r(ys, yc, x, n, rs[rr]);
            arb_sin(e, xa, prec);
            u = _ulps(ys, n + 1, e, n);
            if (u > (double) FIXED_SIN_COS_NOTAB_MAX_ERR)
                TEST_FUNCTION_FAIL("sin r = %d: n = %wd, "
                    "zz = %wd, ulp = %f\n", rs[rr], n, zz, u);
            arb_cos(e, xa, prec);
            u = _ulps(yc, n + 1, e, n);
            if (u > (double) FIXED_SIN_COS_NOTAB_MAX_ERR)
                TEST_FUNCTION_FAIL("cos r = %d: n = %wd, "
                    "zz = %wd, ulp = %f\n", rs[rr], n, zz, u);
        }

        arb_clear(xa); arb_clear(e); fmpz_clear(f);
        flint_free(x); flint_free(y); flint_free(ys); flint_free(yc);
    }

    /* one deep case through the bit-burst regime and the Newton
       square root */
    {
        slong n = 2500;
        slong prec = FLINT_BITS * n + 128;
        nn_ptr x, y, ys, yc;
        arb_t xa, e;
        fmpz_t f;
        double u;

        x = flint_malloc(n * sizeof(ulong));
        y = flint_malloc((n + 1) * sizeof(ulong));
        ys = flint_malloc((n + 1) * sizeof(ulong));
        yc = flint_malloc((n + 1) * sizeof(ulong));
        arb_init(xa); arb_init(e); fmpz_init(f);

        flint_mpn_urandomb(x, state, FLINT_BITS * n);
        x[n - 1] >>= 1;

        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);

        fixed_exp_notab(y, x, n);
        arb_exp(e, xa, prec);
        u = _ulps(y, n + 1, e, n);
        if (u > (double) FIXED_EXP_NOTAB_MAX_ERR)
            TEST_FUNCTION_FAIL("deep exp: ulp = %f\n", u);

        fixed_sin_cos_notab(ys, yc, x, n);
        arb_sin(e, xa, prec);
        u = _ulps(ys, n + 1, e, n);
        if (u > (double) FIXED_SIN_COS_NOTAB_MAX_ERR)
            TEST_FUNCTION_FAIL("deep sin: ulp = %f\n", u);
        arb_cos(e, xa, prec);
        u = _ulps(yc, n + 1, e, n);
        if (u > (double) FIXED_SIN_COS_NOTAB_MAX_ERR)
            TEST_FUNCTION_FAIL("deep cos: ulp = %f\n", u);

        arb_clear(xa); arb_clear(e); fmpz_clear(f);
        flint_free(x); flint_free(y); flint_free(ys); flint_free(yc);
    }

    TEST_FUNCTION_END(state);
}
