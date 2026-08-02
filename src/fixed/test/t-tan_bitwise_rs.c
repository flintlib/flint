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

TEST_FUNCTION_START(fixed_tan_bitwise_rs, state)
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

        fixed_tan_bitwise_rs(res, x, n, r);

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

        if (u > (double) FIXED_TAN_BITWISE_RS_MAX_ERR(n, r))
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

#if FLINT_BITS == 64
    /* the public per-size wrappers fixed_tan_opt_<n> route to the
       same specialized bodies as the r = 0 dispatch, so the two
       must agree bit for bit */
    {
        static void (* const tab[])(nn_ptr, nn_srcptr) = {
            NULL,
            fixed_tan_opt_1, fixed_tan_opt_2, fixed_tan_opt_3,
            fixed_tan_opt_4, fixed_tan_opt_5, fixed_tan_opt_6,
            fixed_tan_opt_7, fixed_tan_opt_8, fixed_tan_opt_9,
            fixed_tan_opt_10, fixed_tan_opt_11, fixed_tan_opt_12
        };
        slong n;

        for (n = 1; n <= 12; n++)
        {
            ulong x[12], t1[13], t2[13];

            flint_mpn_rrandom(x, state, n);
            tab[n](t1, x);
            fixed_tan_bitwise_rs(t2, x, n, 0);
            if (mpn_cmp(t1, t2, n + 1) != 0)
                TEST_FUNCTION_FAIL("opt wrapper mismatch: "
                    "n = %wd\n", n);
        }
    }
#endif

    TEST_FUNCTION_END(state);
}
