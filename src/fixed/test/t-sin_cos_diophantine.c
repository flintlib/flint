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
#include "fmpq.h"
#include "arb.h"
#include "fixed.h"

/* fixed_sin_cos_diophantine agrees with arb within its budget for x
   anywhere in [0, 1), across sizes, weight budgets and prime counts,
   and the cached angles are exact floors */

static double
_ulps_sc(nn_srcptr y, slong ylen, const arb_t exact, slong n)
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

TEST_FUNCTION_START(fixed_sin_cos_diophantine, state)
{
    slong iter;

    for (iter = 0; iter < 20 + 5 * flint_test_multiplier(); iter++)
    {
        slong nv = 1 + n_randint(state, (iter % 5 == 0) ? 300 : 30);
        slong num = 2 + n_randint(state, (iter % 3 == 0) ? 40 : 13);
        slong j = n_randint(state, num);
        arb_t l;
        fmpz_t f, g;
        fmpq_t q;

        if (iter % 7 == 0)
            _fixed_atan_gauss_clear();

        _fixed_atan_gauss_ensure(num, nv);

        arb_init(l); fmpz_init(f); fmpz_init(g); fmpq_init(q);
        fmpq_set_si(q, _fixed_gaussian_primes[2 * j + 1],
            _fixed_gaussian_primes[2 * j]);
        arb_set_fmpq(l, q, FLINT_BITS * nv + 200);
        arb_atan(l, l, FLINT_BITS * nv + 200);
        arb_mul_2exp_si(l, l, 1);
        arb_mul_2exp_si(l, l, FLINT_BITS * nv);
        arb_floor(l, l, FLINT_BITS * nv + 200);
        if (!arb_get_unique_fmpz(f, l))
            TEST_FUNCTION_FAIL("floor not determined: j = %wd, nv = %wd\n",
                j, nv);
        fmpz_set_ui_array(g, _fixed_atan_gauss_entry(j, nv), nv + 1);
        if (!fmpz_equal(f, g))
            TEST_FUNCTION_FAIL("cache: j = %wd, nv = %wd\n%{fmpz}\n%{fmpz}\n",
                j, nv, f, g);
        arb_clear(l); fmpz_clear(f); fmpz_clear(g); fmpq_clear(q);
    }

    for (iter = 0; iter < 60 + 60 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, (iter % 17 == 0) ? 800 : 60);
        slong zz, k;
        slong prec = FLINT_BITS * n + 128;
        nn_ptr x, ys, yc;
        arb_t xa, s, c;
        fmpz_t f;
        double u, w;

        x = flint_malloc(n * sizeof(ulong));
        ys = flint_malloc((n + 1) * sizeof(ulong));
        yc = flint_malloc((n + 1) * sizeof(ulong));
        arb_init(xa); arb_init(s); arb_init(c); fmpz_init(f);

        flint_mpn_urandomb(x, state, FLINT_BITS * n);
        if (iter % 11 == 4)
            flint_mpn_store(x, n, ~UWORD(0));
        if (iter % 13 == 6)
        {
            flint_mpn_zero(x, n);
            x[n - 1] = n_randtest(state);
        }
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
        arb_sin_cos(s, c, xa, prec);

        fixed_sin_cos_diophantine(ys, yc, x, n);
        u = FLINT_MAX(_ulps_sc(ys, n + 1, s, n), _ulps_sc(yc, n + 1, c, n));
        if (u > (double) FIXED_SIN_COS_DIOPHANTINE_MAX_ERR)
            TEST_FUNCTION_FAIL("n = %wd, zz = %wd, ulp = %f\n", n, zz, u);

        /* one output only, and the tangent */
        fixed_sin_cos_diophantine(ys, NULL, x, n);
        fixed_sin_cos_diophantine(NULL, yc, x, n);
        u = FLINT_MAX(_ulps_sc(ys, n + 1, s, n), _ulps_sc(yc, n + 1, c, n));
        if (u > (double) FIXED_SIN_COS_DIOPHANTINE_MAX_ERR)
            TEST_FUNCTION_FAIL("single: n = %wd, zz = %wd, ulp = %f\n", n, zz, u);
        {
            arb_t ta;
            arb_init(ta);
            arb_tan(ta, xa, prec);
            fixed_tan_diophantine(ys, x, n);
            u = _ulps_sc(ys, n + 1, ta, n);
            if (u > (double) FIXED_SIN_COS_DIOPHANTINE_MAX_ERR)
                TEST_FUNCTION_FAIL("tan: n = %wd, zz = %wd, ulp = %f\n", n, zz, u);
            if (n_randint(state, 2))
            {
                _fixed_tan_diophantine_tune(ys, x, n, 2 + n_randint(state, 19),
                    FLINT_BITS * n * (n_randint(state, 40) / 10.0));
                u = _ulps_sc(ys, n + 1, ta, n);
                if (u > (double) FIXED_SIN_COS_DIOPHANTINE_MAX_ERR)
                    TEST_FUNCTION_FAIL("tan tune: n = %wd, zz = %wd, ulp = %f\n", n, zz, u);
            }
            arb_clear(ta);
        }

        for (k = 0; k < 3; k++)
        {
            slong num;

            switch (k)
            {
                case 0: w = FLINT_BITS * n * 0.05; break;
                case 1: w = FLINT_BITS * n * 2.0; break;
                default: w = 1e30; break;
            }
            if (n_randint(state, 4) == 0)
                w = FLINT_BITS * n * (n_randint(state, 100) / 10.0);

            num = 13;
            if (n_randint(state, 3) == 0)
                num = 2 + n_randint(state, 19);

            _fixed_sin_cos_diophantine_tune(ys, yc, x, n, num, w);
            u = FLINT_MAX(_ulps_sc(ys, n + 1, s, n), _ulps_sc(yc, n + 1, c, n));
            if (u > (double) FIXED_SIN_COS_DIOPHANTINE_MAX_ERR)
                TEST_FUNCTION_FAIL("n = %wd, zz = %wd, num = %wd, w = %f, "
                    "ulp = %f\n", n, zz, num, w, u);
        }

        arb_clear(xa); arb_clear(s); arb_clear(c); fmpz_clear(f);
        flint_free(x); flint_free(ys); flint_free(yc);
    }

    TEST_FUNCTION_END(state);
}
