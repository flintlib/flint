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
#include "fixed.h"

static void
mask_top2(nn_ptr x, slong n, slong zbits)
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

/* |(y, ylen) * 2^FLINT_BITS - (ref, ylen + 1)| <= bhi * 2^FLINT_BITS + blo */
static int
close_1ulp(nn_srcptr y, slong ylen, nn_srcptr ref, ulong bhi, ulong blo)
{
    ulong r2[300], d[300];
    slong i;

    r2[0] = 0;
    flint_mpn_copyi(r2 + 1, y, ylen);

    if (mpn_cmp(r2, ref, ylen + 1) >= 0)
        mpn_sub_n(d, r2, ref, ylen + 1);
    else
        mpn_sub_n(d, ref, r2, ylen + 1);

    for (i = 2; i < ylen + 1; i++)
        if (d[i] != 0)
            return 0;
    return (d[1] < bhi || (d[1] == bhi && d[0] <= blo));
}

TEST_FUNCTION_START(fixed_trig_rs, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, (iter % 5 == 0) ? 120 : 24);
        int zbits = (iter % 2) ? 64 : 32;
        int alt = (iter % 4) < 2;
        ulong x[130], xr[131];
        ulong s1[132], c1[132], s2[132], c2[132], refs[133], refc[133];
        ulong a1[131], refa[132];
        ulong bhi, blo;

        flint_mpn_rrandom(x, state, n);
        mask_top2(x, n, zbits);

        xr[0] = 0;
        flint_mpn_copyi(xr + 1, x, n);

        /* sin/cos family: combined and single-output versions
           (zbits shapes the input to exercise both internal paths) */
        if (alt) fixed_sin_cos_rs(s1, c1, x, n);
        else fixed_sinh_cosh_rs(s1, c1, x, n);
        if (alt) { fixed_sin_rs(s2, x, n); fixed_cos_rs(c2, x, n); }
        else { fixed_sinh_rs(s2, x, n); fixed_cosh_rs(c2, x, n); }

        if (mpn_cmp(s1, s2, n + 1) != 0 || mpn_cmp(c1, c2, n + 1) != 0)
            TEST_FUNCTION_FAIL("single/combined mismatch: zbits = %d, "
                "n = %wd, alt = %d\n", zbits, n, alt);

        _fixed_sin_cos_rs_fallback(refs, refc, xr, n + 1, alt);

        bhi = FIXED_SIN_COS_RS_MAX_ERR(n);
        blo = FIXED_SIN_COS_RS_MAX_ERR(n + 1);

        if (!close_1ulp(s1, n + 1, refs, bhi, blo))
            TEST_FUNCTION_FAIL("sin(h): zbits = %d, n = %wd, alt = %d\n",
                zbits, n, alt);
        if (!close_1ulp(c1, n + 1, refc, bhi, blo))
            TEST_FUNCTION_FAIL("cos(h): zbits = %d, n = %wd, alt = %d\n",
                zbits, n, alt);

        /* atan family */
        if (alt) fixed_atan_rs(a1, x, n);
        else fixed_atanh_rs(a1, x, n);

        _fixed_atan_rs_fallback(refa, xr, n + 1, alt);

        bhi = FIXED_ATAN_RS_MAX_ERR(n);
        blo = FIXED_ATAN_RS_MAX_ERR(n + 1);

        if (!close_1ulp(a1, n, refa, bhi, blo))
            TEST_FUNCTION_FAIL("atan(h): zbits = %d, n = %wd, alt = %d\n",
                zbits, n, alt);

        /* absolute check for small n (note: this will become circular
           once arb is built on the fixed module; it can then be
           replaced by an mpfr comparison) */
        if (n <= 16)
        {
            arb_t xa, f, va, delta;
            fmpz_t t;
            mag_t mb;
            slong prec = FLINT_BITS * n + 128;
            slong which;

            arb_init(xa); arb_init(f); arb_init(va); arb_init(delta);
            fmpz_init(t); mag_init(mb);

            fmpz_set_ui_array(t, x, n);
            arb_set_fmpz(xa, t);
            arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);

            for (which = 0; which < 3; which++)
            {
                nn_srcptr y = (which == 0) ? s1 : (which == 1) ? c1 : a1;
                slong ylen = (which == 2) ? n : n + 1;
                double u;

                if (which == 0)
                    (alt ? arb_sin : arb_sinh)(f, xa, prec);
                else if (which == 1)
                    (alt ? arb_cos : arb_cosh)(f, xa, prec);
                else
                    (alt ? arb_atan : arb_atanh)(f, xa, prec);

                fmpz_set_ui_array(t, y, ylen);
                arb_set_fmpz(va, t);
                arb_mul_2exp_si(va, va, -FLINT_BITS * n);

                arb_sub(delta, va, f, ARF_PREC_EXACT);
                arb_get_mag(mb, delta);
                mag_mul_2exp_si(mb, mb, FLINT_BITS * n);
                u = mag_get_d(mb);

                if (u > (double) FIXED_SIN_COS_RS_MAX_ERR(n))
                {
                    arb_clear(xa); arb_clear(f); arb_clear(va);
                    arb_clear(delta); fmpz_clear(t); mag_clear(mb);
                    TEST_FUNCTION_FAIL("arb: which = %wd, zbits = %d, "
                        "n = %wd, alt = %d, ulp = %f\n",
                        which, zbits, n, alt, u);
                }
            }

            arb_clear(xa); arb_clear(f); arb_clear(va); arb_clear(delta);
            fmpz_clear(t); mag_clear(mb);
        }
    }

    TEST_FUNCTION_END(state);
}
