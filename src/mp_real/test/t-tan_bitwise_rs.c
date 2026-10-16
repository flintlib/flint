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

TEST_FUNCTION_START(mp_real_tan_bitwise_rs, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        /* the specialized sizes take r = 0; larger ones divide sin by
           cos and so accept any legal r */
        slong nmin = (FLINT_BITS == 64) ? 1 : 2;
        slong n = nmin + n_randint(state, (iter % 97 == 1) ? 1800
            : (iter % 5 == 0) ? 40 : 9);
        int r = (iter % 3 == 0) ? 0
              : (iter % 3 == 1) ? 32 + (int) n_randint(state, 200)
                                : 32 + (int) n_randint(state, 17);
        nn_ptr x, res;
        arb_t xa, e, va, d;
        fmpz_t f;
        mag_t mb;
        slong prec = FLINT_BITS * n + 128;
        double u;

        if (n > 40)
            n = 40;

x = flint_malloc((n + 2) * sizeof(ulong));
        res = flint_malloc((n + 2) * sizeof(ulong));

                flint_mpn_rrandom(x, state, n);
        if (iter % 13 == 3)
            flint_mpn_zero(x, n);
        if (iter % 7 == 1)
            flint_mpn_store(x, n, ~UWORD(0));

        _mp_real_tan_bitwise_rs(res, &test_err_out, x, n, r);

        arb_init(xa);
        arb_init(e);
        arb_init(va);
        arb_init(d);
        fmpz_init(f);
        mag_init(mb);

        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(xa, f);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
        arb_tan(e, xa, prec);

        fmpz_set_ui_array(f, res, n + 1);
        arb_set_fmpz(va, f);
        arb_mul_2exp_si(va, va, -FLINT_BITS * n);
        arb_sub(d, va, e, ARF_PREC_EXACT);
        arb_get_mag(mb, d);
        mag_mul_2exp_si(mb, mb, FLINT_BITS * n);
        u = mag_get_d(mb);

        if (u > TEST_ERR(DOC_TAN_BITWISE_RS_MAX_ERR(n, r)))
            TEST_FUNCTION_FAIL("n = %wd, r = %d, ulp error = %f\n", n, r, u);

        arb_clear(xa);
        arb_clear(e);
        arb_clear(va);
        arb_clear(d);
        fmpz_clear(f);
        mag_clear(mb);
        flint_free(x);
        flint_free(res);
    }

    /* r = 0 at the per-size sizes (on 64-bit machines the sine and
       cosine kernels and a division) */
    {
        slong n, k;

        /* (n = 1 is unsupported on 32-bit limbs) */
        for (n = (FLINT_BITS == 64) ? 1 : 2; n <= 12; n++)
        {
            for (k = 0; k < 16; k++)
            {
                ulong x[12], t1[13];
                arb_t xa, e;
                double u;

                flint_mpn_rrandom(x, state, n);
                if (k == 0)
                    flint_mpn_zero(x, n);
                if (k == 1)
                    flint_mpn_store(x, n, ~UWORD(0));

                _mp_real_tan_bitwise_rs(t1, &test_err_out, x, n, 0);

                arb_init(xa);
                arb_init(e);
                _test_bitwise_arb(xa, x, n);
                arb_tan(e, xa, FLINT_BITS * n + 128);
                u = _test_bitwise_ulps(t1, n + 1, n, e);
                arb_clear(xa);
                arb_clear(e);

                if (u > TEST_ERR(DOC_TAN_BITWISE_RS_MAX_ERR(n, 0)))
                    TEST_FUNCTION_FAIL("r = 0: n = %wd, ulp = %f\n", n, u);
            }
        }
    }

    /* arguments at the decision boundaries (greedy and signed digits) */
    for (iter = 0; iter < 60 * flint_test_multiplier(); iter++)
    {
        slong n;
        int r;
        nn_ptr x, res;
        arb_t xa, e;
        double u;

        _test_bitwise_boundary_nr(&n, &r, iter, state);
        x = flint_malloc(n * sizeof(ulong));
        res = flint_malloc((n + 1) * sizeof(ulong));
        _test_bitwise_boundary_arg(x, n,
            r ? r : _mp_real_trig_bitwise_rs_default_r(n), BOUNDARY_TRIG,
            state);

        _mp_real_tan_bitwise_rs(res, &test_err_out, x, n, r);

        arb_init(xa);
        arb_init(e);
        _test_bitwise_arb(xa, x, n);
        arb_tan(e, xa, FLINT_BITS * n + 128);
        u = _test_bitwise_ulps(res, n + 1, n, e);
        arb_clear(xa);
        arb_clear(e);

        if (u > TEST_ERR(DOC_TAN_BITWISE_RS_MAX_ERR(n, r)))
            TEST_FUNCTION_FAIL("boundary: n = %wd, r = %d, ulp = %f\n", n, r, u);

        flint_free(x);
        flint_free(res);
    }

    TEST_FUNCTION_END(state);
}
