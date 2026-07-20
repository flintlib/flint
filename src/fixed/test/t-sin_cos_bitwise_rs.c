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

TEST_FUNCTION_START(fixed_sin_cos_bitwise_rs, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong det_n = 1 + (iter / 32) % 9;
        slong det_k = iter % 32;
        /* the reduction runs to r + 1, so the boundary form of the
           final step needs r one BELOW a multiple of FLINT_BITS */
        int det_r = (det_k < 8) ? (int) (FLINT_BITS * (det_k + 1))
                  : (det_k < 16) ? (int) (FLINT_BITS * (det_k - 7) + 1)
                  : (det_k < 24) ? (int) (FLINT_BITS * (det_k - 15) - 1)
                  : (det_k == 24) ? 32
                  : (det_k == 25) ? 16
                  : (det_k == 26) ? 0
                  : (int) (FLINT_BITS * det_n - 16 - (det_k - 27));
        int det = (iter < 32 * 9);
        /* n = 1 is unsupported on 32-bit limbs */
        slong nmin = (FLINT_BITS == 64) ? 1 : 2;
        slong n = det ? FLINT_MAX(det_n, nmin)
            : nmin + n_randint(state, (iter % 97 == 1) ? 1800
                : (iter % 5 == 0) ? 80 : 16);
        int r = det ? FLINT_MAX(det_r, ((det_r == 0) ? 0 : 32))
            : (iter % 4 == 0) ? 0 :
            (iter % 9 == 3) ? 32 + (int) n_randint(state, 17) :
            (iter % 7 == 2) ? FLINT_BITS * (1 + (int) n_randint(state, 7))
                            : 32 + (int) n_randint(state, 401);
        nn_ptr x, ysin, ycos;
        arb_t xa, e, va, dd;
        fmpz_t f;
        mag_t mb;
        slong prec = FLINT_BITS * n + 128;
        double u;
        int which;

        x = flint_malloc((n + 2) * sizeof(ulong));
        ysin = flint_malloc((n + 2) * sizeof(ulong));
        ycos = flint_malloc((n + 2) * sizeof(ulong));

        flint_mpn_rrandom(x, state, n);
        if (iter % 13 == 3)
            flint_mpn_zero(x, n);
        if (iter % 7 == 1)
            flint_mpn_store(x, n, ~UWORD(0));

        /* alternate NULL outputs to cover the partial interfaces */
        if (iter % 11 == 5)
            fixed_sin_cos_bitwise_rs(ysin, NULL, x, n, r),
            fixed_sin_cos_bitwise_rs(NULL, ycos, x, n, r);
        else
            fixed_sin_cos_bitwise_rs(ysin, ycos, x, n, r);

        arb_init(xa);
        arb_init(e);
        arb_init(va);
        arb_init(dd);
        fmpz_init(f);
        mag_init(mb);

        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);

        for (which = 0; which < 2; which++)
        {
            if (which == 0)
                arb_sin(e, xa, prec);
            else
                arb_cos(e, xa, prec);

            fmpz_set_ui_array(f, (which == 0) ? ysin : ycos, n + 1);
            arb_set_fmpz(va, f);
            arb_mul_2exp_si(va, va, -FLINT_BITS * n);
            arb_sub(dd, va, e, ARF_PREC_EXACT);
            arb_get_mag(mb, dd);
            mag_mul_2exp_si(mb, mb, FLINT_BITS * n);
            u = mag_get_d(mb);

            if (u > (double) FIXED_SIN_COS_BITWISE_RS_MAX_ERR(n, r))
                TEST_FUNCTION_FAIL(
                    "%s: n = %wd, r = %d, ulp error = %f\n",
                    (which == 0) ? "sin" : "cos", n, r, u);
        }

        arb_clear(xa);
        arb_clear(e);
        arb_clear(va);
        arb_clear(dd);
        fmpz_clear(f);
        mag_clear(mb);
        flint_free(x);
        flint_free(ysin);
        flint_free(ycos);
    }

#if FLINT_BITS == 64
    /* the public per-size wrappers fixed_sin_cos_opt_<n> route to
       the same specialized bodies as the r = 0 dispatch, so the two
       must agree bit for bit; this pins the wrapper wiring and the
       NULL-output entry shapes */
    {
        static void (* const tab[])(nn_ptr, nn_ptr, nn_srcptr) = {
            NULL,
            fixed_sin_cos_opt_1, fixed_sin_cos_opt_2,
            fixed_sin_cos_opt_3, fixed_sin_cos_opt_4,
            fixed_sin_cos_opt_5, fixed_sin_cos_opt_6,
            fixed_sin_cos_opt_7, fixed_sin_cos_opt_8,
            fixed_sin_cos_opt_9, fixed_sin_cos_opt_10,
            fixed_sin_cos_opt_11, fixed_sin_cos_opt_12
        };
        slong n;

        for (n = 1; n <= 12; n++)
        {
            ulong x[12], s1[13], c1[13], s2[13], c2[13];

            flint_mpn_rrandom(x, state, n);
            tab[n](s1, c1, x);
            fixed_sin_cos_bitwise_rs(s2, c2, x, n, 0);
            if (mpn_cmp(s1, s2, n + 1) != 0
                || mpn_cmp(c1, c2, n + 1) != 0)
                TEST_FUNCTION_FAIL("opt wrapper mismatch: "
                    "n = %wd\n", n);
        }
    }
#endif

    TEST_FUNCTION_END(state);
}
