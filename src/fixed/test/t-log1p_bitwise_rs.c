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

TEST_FUNCTION_START(fixed_log1p_bitwise_rs, state)
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
            : nmin + n_randint(state, (iter % 97 == 1) ? 1800
                : (iter % 5 == 0) ? 80 : 16);
        int r = det ? FLINT_MAX(det_r, ((det_r == 0) ? 0 : 32))
            : (iter % 4 == 0) ? 0 :
            (iter % 7 == 2) ? FLINT_BITS * (1 + (int) n_randint(state, 7))
                            : 32 + (int) n_randint(state, 401);
        nn_ptr x, res;
        arb_t xa, e, va, delta;
        fmpz_t t;
        mag_t mb;
        slong prec = FLINT_BITS * n + 128;
        double u;

        x = flint_malloc((n + 1) * sizeof(ulong));
        res = flint_malloc((n + 1) * sizeof(ulong));

        flint_mpn_rrandom(x, state, n);
        if (iter % 13 == 3)
            flint_mpn_zero(x, n);

        fixed_log1p_bitwise_rs(res, x, n, r);

        /* absolute check (note: this will become circular once arb is
           built on the fixed module; it can then be replaced by an
           mpfr comparison) */
        arb_init(xa); arb_init(e); arb_init(va); arb_init(delta);
        fmpz_init(t); mag_init(mb);

        fmpz_set_ui_array(t, x, n);
        arb_set_fmpz(xa, t);
        arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
        arb_add_ui(xa, xa, 1, prec);
        arb_log(e, xa, prec);

        fmpz_set_ui_array(t, res, n);
        arb_set_fmpz(va, t);
        arb_mul_2exp_si(va, va, -FLINT_BITS * n);

        arb_sub(delta, va, e, ARF_PREC_EXACT);
        arb_get_mag(mb, delta);
        mag_mul_2exp_si(mb, mb, FLINT_BITS * n);
        u = mag_get_d(mb);

        arb_clear(xa); arb_clear(e); arb_clear(va); arb_clear(delta);
        fmpz_clear(t); mag_clear(mb);
        flint_free(x);
        flint_free(res);

        if (u > (double) FIXED_LOG1P_BITWISE_RS_MAX_ERR(n, r))
            TEST_FUNCTION_FAIL("n = %wd, r = %d, ulp = %f\n", n, r, u);
    }

#if FLINT_BITS == 64
    /* dense r = 0 sampling at the fully specialized sizes (the
       straight-line code's reduction and borrow branches are
       input-dependent) */
    {
        slong n, k;

        for (n = 1; n <= 7; n++)
        {
            for (k = 0; k < 24; k++)
            {
                ulong x[7], res[7];
                arb_t xa, e, va, delta;
                fmpz_t t;
                mag_t mb;
                double u;
                slong prec = FLINT_BITS * n + 128;

                flint_mpn_rrandom(x, state, n);
                if (k == 0)
                    flint_mpn_zero(x, n);
                if (k == 1)
                    flint_mpn_store(x, n, ~UWORD(0));

                fixed_log1p_bitwise_rs(res, x, n, 0);

                arb_init(xa); arb_init(e); arb_init(va);
                arb_init(delta); fmpz_init(t); mag_init(mb);
                fmpz_set_ui_array(t, x, n);
                arb_set_fmpz(xa, t);
                arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
                arb_log1p(e, xa, prec);
                fmpz_set_ui_array(t, res, n);
                arb_set_fmpz(va, t);
                arb_mul_2exp_si(va, va, -FLINT_BITS * n);
                arb_sub(delta, va, e, ARF_PREC_EXACT);
                arb_get_mag(mb, delta);
                mag_mul_2exp_si(mb, mb, FLINT_BITS * n);
                u = mag_get_d(mb);
                arb_clear(xa); arb_clear(e); arb_clear(va);
                arb_clear(delta); fmpz_clear(t); mag_clear(mb);

                if (u > (double) FIXED_LOG1P_BITWISE_RS_MAX_ERR(n,
                    fixed_log1p_bitwise_rs_default_r(n)))
                    TEST_FUNCTION_FAIL("r = 0: n = %wd, "
                        "ulp = %f\n", n, u);
            }
        }
    }
#endif

    TEST_FUNCTION_END(state);
}
