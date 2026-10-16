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
#include "bitwise_boundary.h"

TEST_FUNCTION_START(mp_real_sin_cos_bitwise_rs, state)
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
            _mp_real_sin_cos_bitwise_rs(ysin, NULL, &test_err_out, x, n, r),
            _mp_real_sin_cos_bitwise_rs(NULL, ycos, &test_err_out, x, n, r);
        else
            _mp_real_sin_cos_bitwise_rs(ysin, ycos, &test_err_out, x, n, r);

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

            if (u > TEST_ERR(DOC_SIN_COS_BITWISE_RS_MAX_ERR(n, r)))
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

    /* r = 0 at the per-size sizes (the small.h and dense-table kernels
       on 64-bit machines), with the NULL-output entry shapes, which must
       agree bit for bit with the full call */
    {
        slong n, k;

        /* (n = 1 is unsupported on 32-bit limbs) */
        for (n = (FLINT_BITS == 64) ? 1 : 2; n <= 12; n++)
        {
            for (k = 0; k < 16; k++)
            {
                ulong x[12], s1[13], c1[13], s2[13], c2[13];
                arb_t xa, es, ec;
                double u, v;
                ulong e1;

                flint_mpn_rrandom(x, state, n);
                if (k == 0)
                    flint_mpn_zero(x, n);
                if (k == 1)
                    flint_mpn_store(x, n, ~UWORD(0));

                _mp_real_sin_cos_bitwise_rs(s1, c1, &test_err_out, x, n, 0);
                e1 = test_err_out;
                _mp_real_sin_cos_bitwise_rs(s2, NULL, &test_err_out, x, n, 0);
                _mp_real_sin_cos_bitwise_rs(NULL, c2, &test_err_out, x, n, 0);
                if (mpn_cmp(s1, s2, n + 1) != 0 || mpn_cmp(c1, c2, n + 1) != 0
                        || test_err_out != e1)
                    TEST_FUNCTION_FAIL("r = 0, NULL output: n = %wd\n", n);

                arb_init(xa);
                arb_init(es);
                arb_init(ec);
                _test_bitwise_arb(xa, x, n);
                arb_sin_cos(es, ec, xa, FLINT_BITS * n + 128);
                u = _test_bitwise_ulps(s1, n + 1, n, es);
                v = _test_bitwise_ulps(c1, n + 1, n, ec);
                arb_clear(xa);
                arb_clear(es);
                arb_clear(ec);

                if (u > TEST_ERR(DOC_SIN_COS_BITWISE_RS_MAX_ERR(n, 0))
                        || v > TEST_ERR(DOC_SIN_COS_BITWISE_RS_MAX_ERR(n, 0)))
                    TEST_FUNCTION_FAIL("r = 0: n = %wd, ulp = %f, %f\n",
                        n, u, v);
            }
        }
    }

    /* arguments at the decision boundaries (greedy and signed digits) */
    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        slong n;
        int r;
        nn_ptr x, ys, yc;
        arb_t xa, es, ec;
        double u, v;

        _test_bitwise_boundary_nr(&n, &r, iter, state);
        x = flint_malloc(n * sizeof(ulong));
        ys = flint_malloc((n + 1) * sizeof(ulong));
        yc = flint_malloc((n + 1) * sizeof(ulong));
        _test_bitwise_boundary_arg(x, n,
            r ? r : _mp_real_trig_bitwise_rs_default_r(n), BOUNDARY_TRIG,
            state);

        _mp_real_sin_cos_bitwise_rs(ys, yc, &test_err_out, x, n, r);

        arb_init(xa);
        arb_init(es);
        arb_init(ec);
        _test_bitwise_arb(xa, x, n);
        arb_sin_cos(es, ec, xa, FLINT_BITS * n + 128);
        u = _test_bitwise_ulps(ys, n + 1, n, es);
        v = _test_bitwise_ulps(yc, n + 1, n, ec);
        arb_clear(xa);
        arb_clear(es);
        arb_clear(ec);

        if (u > TEST_ERR(DOC_SIN_COS_BITWISE_RS_MAX_ERR(n, r))
                || v > TEST_ERR(DOC_SIN_COS_BITWISE_RS_MAX_ERR(n, r)))
            TEST_FUNCTION_FAIL("boundary: n = %wd, r = %d, ulp = %f, %f\n",
                n, r, u, v);

        flint_free(x);
        flint_free(ys);
        flint_free(yc);
    }

    TEST_FUNCTION_END(state);
}
