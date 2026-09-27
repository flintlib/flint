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
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

/* _mp_real_series_rs and _mp_real_series_rs_sin_cos against arb: within
   5 ulps for exp and exp(-x), 2 ulps otherwise, at any r >= 8,
   including aliasing; _mp_real_series_rs_tan within 3 ulps */
TEST_FUNCTION_START(mp_real_series_rs, state)
{
    slong iter;
    arb_t t, y;
    fmpz_t f;

    arb_init(t);
    arb_init(y);
    fmpz_init(f);

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, (iter % 8 == 0) ? 160 : 24);
        slong r = 8 + n_randint(state, (iter % 2) ? 56 : 300);
        slong prec = FLINT_BITS * n + 128, part, nparts;
        int func = 1 + n_randint(state, 9), alias = n_randint(state, 2);
        ulong x[170], y1[340], y2[170];

        if (r > FLINT_BITS * n - 1)
            r = FLINT_BITS * n - 1;
        if (r < 8)
            continue;

        /* x < 2^-r, often near the bound */
        fmpz_randbits(f, state, FLINT_BITS * n);
        fmpz_abs(f, f);
        fmpz_fdiv_q_2exp(f, f, r + n_randint(state, 3));
        if (n_randint(state, 4) == 0)
        {
            fmpz_one(f);
            fmpz_mul_2exp(f, f, FLINT_BITS * n - r);
            fmpz_sub_ui(f, f, 1 + n_randint(state, 5));
        }
        flint_mpn_zero(x, n);
        fmpz_get_ui_array(x, n, f);

        if (func == 9)
        {
            nparts = 2;
            if (alias)
            {
                flint_mpn_copyi(y1, x, n);
                _mp_real_series_rs_sin_cos(y1, y2, y1, n, r, iter & 1);
            }
            else
                _mp_real_series_rs_sin_cos(y1, y2, x, n, r, iter & 1);
        }
        else
        {
            nparts = 1;
            if (alias)
            {
                flint_mpn_copyi(y1, x, n);
                _mp_real_series_rs(y1, y1, n, r, func);
            }
            else
                _mp_real_series_rs(y1, x, n, r, func);
        }

        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(t, f);
        arb_mul_2exp_si(t, t, -FLINT_BITS * n);

        for (part = 0; part < nparts; part++)
        {
            int fu = func;
            nn_srcptr res = y1;
            slong rn = n;
            double err, bound = 2.0;

            if (func == 9)
            {
                fu = part ? ((iter & 1) ? MP_REAL_SERIES_COSH : MP_REAL_SERIES_COS)
                          : ((iter & 1) ? MP_REAL_SERIES_SINH : MP_REAL_SERIES_SIN);
                res = part ? y2 : y1;
            }

            switch (fu)
            {
                case MP_REAL_SERIES_SIN: arb_sin(y, t, prec); break;
                case MP_REAL_SERIES_SINH: arb_sinh(y, t, prec); break;
                case MP_REAL_SERIES_COS:
                    arb_cos(y, t, prec); arb_sub_ui(y, y, 1, prec); arb_neg(y, y); break;
                case MP_REAL_SERIES_COSH:
                    arb_cosh(y, t, prec); arb_sub_ui(y, y, 1, prec); break;
                case MP_REAL_SERIES_ATAN: arb_atan(y, t, prec); break;
                case MP_REAL_SERIES_ATANH: arb_atanh(y, t, prec); break;
                case MP_REAL_SERIES_EXP:
                    arb_exp(y, t, prec); rn = n + 1; bound = 5.0; break;
                case MP_REAL_SERIES_EXP_NEG:
                    arb_neg(y, t); arb_exp(y, y, prec); rn = n + 1; bound = 5.0; break;
                default:    /* TAN: not a series_rs family */
                    continue;
            }

            arb_mul_2exp_si(y, y, FLINT_BITS * n);
            fmpz_set_ui_array(f, res, rn);
            arb_sub_fmpz(y, y, f, prec);
            err = arf_get_d(arb_midref(y), ARF_RND_NEAR);

            if (!(fabs(err) < bound))
                TEST_FUNCTION_FAIL("func = %d, n = %wd, r = %wd, alias = %d, "
                    "err = %f ulps\n", fu, n, r, alias, err);
        }
    }

    /* the chunked tangent: within 3 ulps wherever it applies */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, (iter % 4 == 0) ? 160 : 80);
        slong r = 32 + n_randint(state, (iter % 2) ? 40 : 300);
        slong prec = FLINT_BITS * n + 128;
        ulong x[170], y1[170];
        double err;

        if (r > FLINT_BITS * n - 1 || !_mp_real_series_rs_tan_ok(n, r))
            continue;

        fmpz_randbits(f, state, FLINT_BITS * n);
        fmpz_abs(f, f);
        fmpz_fdiv_q_2exp(f, f, r + n_randint(state, 3));
        if (n_randint(state, 3) == 0)
        {
            fmpz_one(f);
            fmpz_mul_2exp(f, f, FLINT_BITS * n - r);
            fmpz_sub_ui(f, f, 1 + n_randint(state, 5));
        }
        flint_mpn_zero(x, n);
        fmpz_get_ui_array(x, n, f);

        if (iter & 1)
        {
            flint_mpn_copyi(y1, x, n);
            _mp_real_series_rs_tan(y1, y1, n, r);
        }
        else
            _mp_real_series_rs_tan(y1, x, n, r);

        arb_set_fmpz(t, f);
        arb_mul_2exp_si(t, t, -FLINT_BITS * n);
        arb_tan(y, t, prec);
        arb_mul_2exp_si(y, y, FLINT_BITS * n);
        fmpz_set_ui_array(f, y1, n);
        arb_sub_fmpz(y, y, f, prec);
        err = arf_get_d(arb_midref(y), ARF_RND_NEAR);

        if (!(fabs(err) < 3.0))
            TEST_FUNCTION_FAIL("tan: n = %wd, r = %wd, err = %f ulps\n", n, r, err);
    }

    arb_clear(t);
    arb_clear(y);
    fmpz_clear(f);

    TEST_FUNCTION_END(state);
}
