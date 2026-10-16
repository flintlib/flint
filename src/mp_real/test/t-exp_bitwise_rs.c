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

TEST_FUNCTION_START(mp_real_exp_bitwise_rs, state)
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

        _mp_real_exp_bitwise_rs(res, &test_err_out, x, n, r);

        /* absolute check (note: this will become circular once arb is
           built on the mp_real module; it can then be replaced by an
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

        if (u > TEST_ERR(DOC_EXP_BITWISE_RS_MAX_ERR(n, r)))
            TEST_FUNCTION_FAIL("n = %wd, r = %d, ulp = %f\n", n, r, u);
        flint_free(x);
        flint_free(res);
    }

    /* guaranteed large sizes: at the default multiplier the random
       big-n draws land inside the deterministic prefix and never
       fire, leaving the reconstruction's product regime
       (_mp_real_exp_recon_prod: >= 192 output limbs with at least 8
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

            _mp_real_exp_bitwise_rs(res, &test_err_out, x, n, r);

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

            if (u > TEST_ERR(DOC_EXP_BITWISE_RS_MAX_ERR(n,
                (r == 0) ? _mp_real_exp_bitwise_rs_default_r(n) : r)))
                TEST_FUNCTION_FAIL("big n = %wd, r = %d, "
                    "ulp = %f\n", n, r, u);
            flint_free(x);
            flint_free(res);
        }
    }

    /* arguments at the decision boundaries: sums of table entries
       log(1 + 2^-i) (the greedy reduction's boundaries), plus a few
       ulps, or plus 0.75 2^-i (the signed-digit thresholds of
       _mp_real_exp_bitwise_vec) */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong n = 4 + n_randint(state, (iter % 4 == 0) ? 160 : 40);
        int r = (iter % 3 == 0) ? 0 : 32 + (int) n_randint(state, 300);
        slong i, k, prec = FLINT_BITS * n + 128;
        nn_ptr x, res, s;
        arb_t xa, e, va, delta;
        fmpz_t t;
        mag_t mb;
        double u;

        r = FLINT_MIN(r, FLINT_BITS * n - 16);
        x = flint_malloc(n * sizeof(ulong));
        s = flint_malloc(n * sizeof(ulong));
        res = flint_malloc((n + 2) * sizeof(ulong));

        _mp_real_exp_logs_ensure(n, (r == 0) ? 600 : r + 2);
        flint_mpn_zero(x, n);
        k = (r == 0) ? _mp_real_exp_bitwise_rs_default_r(n) : r;
        for (i = 1; i <= k + 1; i += 1 + n_randint(state, 4))
        {
            if (mpn_add_n(s, x, _mp_real_exp_logs_entry(i, n), n) == 0)
                flint_mpn_copyi(x, s, n);
        }
        if (iter % 2)
        {
            /* + 3 2^(-i-2) */
            slong b;
            i = 1 + n_randint(state, FLINT_MIN(k + 1, FLINT_BITS * n - 2));
            flint_mpn_zero(s, n);
            for (b = FLINT_BITS * n - i - 2; b <= FLINT_BITS * n - i - 1; b++)
                mpn_add_1(s + b / FLINT_BITS, s + b / FLINT_BITS,
                    n - b / FLINT_BITS, UWORD(1) << (b % FLINT_BITS));
            if (mpn_add_n(s, x, s, n) == 0)
                flint_mpn_copyi(x, s, n);
        }
        mpn_add_1(x, x, n, n_randint(state, 3));
        if (n_randint(state, 2) && (x[0] != 0 || n > 1))
            mpn_sub_1(x, x, n, n_randint(state, 2));

        _mp_real_exp_bitwise_rs(res, &test_err_out, x, n, r);

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

        if (u > TEST_ERR(DOC_EXP_BITWISE_RS_MAX_ERR(n, k)))
            TEST_FUNCTION_FAIL("boundary: n = %wd, r = %d, ulp = %f\n", n, r, u);

        flint_free(x);
        flint_free(s);
        flint_free(res);
    }

    TEST_FUNCTION_END(state);
}
