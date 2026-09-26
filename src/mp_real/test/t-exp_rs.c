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
#include "max_err.h"

/* enforce x < 2^-zbits */
static void
mask_top(nn_ptr x, slong n, slong zbits)
{
    slong i = n - 1;

    while (zbits >= FLINT_BITS && i >= 0)
    {
        x[i--] = 0;
        zbits -= FLINT_BITS;
    }
    if (zbits > 0 && i >= 0)
        x[i] >>= zbits;
}

/* reference: the tapered series of series_rs.c at the given size */
static void
ref_exp(nn_ptr res, nn_srcptr x, slong n)
{
    slong zb = _mp_real_elem_lzb(x, n);

    if (zb == WORD_MAX)
    {
        flint_mpn_zero(res, n);
        res[n] = 1;
    }
    else
        _mp_real_series_rs(res, x, n, zb, MP_REAL_SERIES_EXP);
}

TEST_FUNCTION_START(mp_real_exp_rs, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, (iter % 5 == 0) ? 250 : 24);
        int zbits = (iter % 2) ? 64 : 32;
        ulong x[260], xr[261], res[261], ref[262], r2[262], d[262];
        ulong bhi, blo;
        slong i;

        flint_mpn_rrandom(x, state, n);
        mask_top(x, n, zbits);

        _mp_real_exp_rs(res, &test_err_out, x, n);

        /* compare against the tapered series at one extra limb */
        xr[0] = 0;
        flint_mpn_copyi(xr + 1, x, n);
        ref_exp(ref, xr, n + 1);

        r2[0] = 0;
        flint_mpn_copyi(r2 + 1, res, n + 1);

        if (mpn_cmp(r2, ref, n + 2) >= 0)
            mpn_sub_n(d, r2, ref, n + 2);
        else
            mpn_sub_n(d, ref, r2, n + 2);

        if (test_err_out > DOC_EXP_RS_MAX_ERR(n))
            TEST_FUNCTION_FAIL("err = %wu above the documented maximum\n", test_err_out);
        bhi = test_err_out;
        blo = DOC_EXP_RS_MAX_ERR(n + 1);

        for (i = 2; i < n + 2; i++)
            if (d[i] != 0)
                TEST_FUNCTION_FAIL("zbits = %d, n = %wd, limb %wd\n",
                    zbits, n, i);
        if (d[1] > bhi || (d[1] == bhi && d[0] > blo))
            TEST_FUNCTION_FAIL("zbits = %d, n = %wd, d1 = %wu, d0 = %wu\n",
                zbits, n, d[1], d[0]);

        /* absolute check (note: this will become circular once arb_exp
           itself is built on the mp_real module; it can then be
           replaced by an mpfr comparison) */
        if (1)
        {
            arb_t xa, e, va, delta;
            fmpz_t t;
            mag_t mb;
            slong prec = FLINT_BITS * n + 128;
            double u;

            arb_init(xa); arb_init(e); arb_init(va); arb_init(delta);
            fmpz_init(t); mag_init(mb);

            fmpz_set_ui_array(t, x, n);
            arb_set_fmpz(xa, t);
            arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
            arb_exp(e, xa, prec);

            fmpz_set_ui_array(t, res, n + 1);
            arb_set_fmpz(va, t);
            arb_mul_2exp_si(va, va, -FLINT_BITS * n);

            arb_sub(delta, va, e, ARF_PREC_EXACT);
            arb_get_mag(mb, delta);
            mag_mul_2exp_si(mb, mb, FLINT_BITS * n);
            u = mag_get_d(mb);

            arb_clear(xa); arb_clear(e); arb_clear(va); arb_clear(delta);
            fmpz_clear(t); mag_clear(mb);

            if (u > TEST_ERR(DOC_EXP_RS_MAX_ERR(n)))
                TEST_FUNCTION_FAIL("arb: zbits = %d, n = %wd, ulp = %f\n",
                    zbits, n, u);
        }
    }

    TEST_FUNCTION_END(state);
}
