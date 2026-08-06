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

TEST_FUNCTION_START(fixed_exp_bitwise_rs, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        /* the first iterations sweep (n, r) deterministically: all
           window counts, both boundary forms of the extra step
           (r at and off multiples of FLINT_BITS), and the smallest
           and largest legal r for each n */
        slong det_n = 1 + (iter / 24) % 9;
        slong det_k = iter % 24;
        int det_r = (det_k < 8) ? (int) (FLINT_BITS * (det_k + 1))
                  : (det_k < 16) ? (int) (FLINT_BITS * (det_k - 7) + 1)
                  : (det_k == 16) ? 32
                  : (det_k == 17) ? 0
                  : (int) (FLINT_BITS * det_n - 16 - (det_k - 18));
        int det = (iter < 24 * 9);

        /* n = 1 is unsupported on 32-bit limbs (r would clamp below
           the series contract) */
        slong nmin = (FLINT_BITS == 64) ? 1 : 2;
        slong n = det ? FLINT_MAX(det_n, nmin)
            : nmin + n_randint(state, (iter % 97 == 1) ? 2400
                : (iter % 5 == 0) ? 80 : 16);
        int r = det ? FLINT_MAX(det_r, ((det_r == 0) ? 0 : 32))
            : (iter % 4 == 0) ? 0 :
            (iter % 7 == 2) ? FLINT_BITS * (1 + (int) n_randint(state, 7))
                            : 32 + (int) n_randint(state, 417);
        nn_ptr x, res;
        arb_t xa, e, va, delta;
        fmpz_t t;
        mag_t mb;
        slong prec = FLINT_BITS * n + 128;
        double u;

        x = flint_malloc(n * sizeof(ulong));
        res = flint_malloc((n + 2) * sizeof(ulong));
        flint_mpn_rrandom(x, state, n);
        if (iter % 13 == 3)
            flint_mpn_zero(x, n);

        fixed_exp_bitwise_rs(res, x, n, r);

        /* absolute check (note: this will become circular once arb is
           built on the fixed module; it can then be replaced by an
           mpfr comparison) */
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

        if (u > (double) FIXED_EXP_BITWISE_RS_MAX_ERR(n, r))
            TEST_FUNCTION_FAIL("n = %wd, r = %d, ulp = %f\n", n, r, u);
        flint_free(x);
        flint_free(res);
    }

    /* guaranteed large sizes: at the default multiplier the random
       big-n draws land inside the deterministic prefix and never
       fire, leaving the reconstruction's product regime
       (_fixed_exp_recon_prod: >= 192 output limbs with at least 8
       pending factors) unexercised */
    {
        int rs[] = { 32, 0 };
        slong c;

        for (c = 0; c < 2; c++)
        {
            slong n = 400;
            int r = rs[c];
            nn_ptr x, res;
            arb_t xa, e, va, delta;
            fmpz_t t;
            mag_t mb;
            double u;
            slong prec = FLINT_BITS * n + 64;

            x = flint_malloc(n * sizeof(ulong));
            res = flint_malloc((n + 1) * sizeof(ulong));
            flint_mpn_rrandom(x, state, n);

            fixed_exp_bitwise_rs(res, x, n, r);

            arb_init(xa); arb_init(e); arb_init(va);
            arb_init(delta); fmpz_init(t); mag_init(mb);
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

            arb_clear(xa); arb_clear(e); arb_clear(va);
            arb_clear(delta); fmpz_clear(t); mag_clear(mb);

            if (u > (double) FIXED_EXP_BITWISE_RS_MAX_ERR(n,
                (r == 0) ? fixed_exp_bitwise_rs_default_r(n) : r))
                TEST_FUNCTION_FAIL("big n = %wd, r = %d, "
                    "ulp = %f\n", n, r, u);
            flint_free(x);
            flint_free(res);
        }
    }

    TEST_FUNCTION_END(state);
}
