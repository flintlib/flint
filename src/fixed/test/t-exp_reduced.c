/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "fixed.h"

/* every algorithm of fixed_exp_reduced agrees with arb within the
   documented budget, for random (wn, r) with t < 2^-r */

TEST_FUNCTION_START(fixed_exp_reduced, state)
{
    slong iter;

    for (iter = 0; iter < 20 + 20 * flint_test_multiplier(); iter++)
    {
        slong wn = 1 + n_randint(state, (iter % 17 == 0) ? 700 : 60);
        flint_bitcnt_t r = 16 + n_randint(state,
            FLINT_MIN(FLINT_BITS * (ulong) wn, 3000));
        slong i, alg;
        nn_ptr t, y;
        arb_t ta, ref, got;
        fmpz_t f;

        t = flint_malloc((wn + 2) * sizeof(ulong));
        y = flint_malloc((wn + 2) * sizeof(ulong));
        arb_init(ta); arb_init(ref); arb_init(got); fmpz_init(f);

        flint_mpn_urandomb(t, state, FLINT_BITS * wn);
        /* clear the top r bits so t < 2^-r */
        for (i = 0; i < (slong) (r / FLINT_BITS) && i < wn; i++)
            t[wn - 1 - i] = 0;
        if ((r % FLINT_BITS) && (slong) (r / FLINT_BITS) < wn)
            t[wn - 1 - r / FLINT_BITS] >>= (r % FLINT_BITS);

        fmpz_set_ui_array(f, t, wn);
        arb_set_fmpz(ta, f);
        arb_mul_2exp_si(ta, ta, -FLINT_BITS * wn);
        arb_exp(ref, ta, FLINT_BITS * wn + 128);

        for (alg = 0; alg <= 4; alg++)
        {
            if (alg == 1 || alg == 2)
                if (r < 32)
                    continue;
            fixed_exp_reduced(y, t, wn, r, (int) alg);

            fmpz_set_ui_array(f, y, wn + 1);
            arb_set_fmpz(got, f);
            arb_mul_2exp_si(got, got, -FLINT_BITS * wn);
            arb_sub(got, got, ref, FLINT_BITS * wn + 128);
            arb_mul_2exp_si(got, got, FLINT_BITS * wn);
            arb_abs(got, got);
            {
                arb_t bnd;
                arb_init(bnd);
                arb_set_ui(bnd, FIXED_EXP_REDUCED_MAX_ERR);
                if (!arb_le(got, bnd))
                {
                    arb_clear(bnd);
                    TEST_FUNCTION_FAIL("alg = %wd, wn = %wd, "
                        "r = %wd\n", alg, wn, (slong) r);
                }
                arb_clear(bnd);
            }
        }

        arb_clear(ta); arb_clear(ref); arb_clear(got); fmpz_clear(f);
        flint_free(t); flint_free(y);
    }

    TEST_FUNCTION_END(state);
}
