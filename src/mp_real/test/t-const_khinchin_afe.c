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
#include "mp_real.h"

/* Khinchin's constant by the theta splitting, with split points
   u0 = 4^-a, a = 1 to 7, on 1 to 3 threads with 1 to 3 dual terms at
   once, against the zeta series (_mp_real_const_khinchin_tune) at a
   higher precision: the balls must overlap and be accurate to within a
   few bits of n limbs.  Now and then a few hundred limbs, where the
   long k-series take the limb-array bottom, the transformed products
   and several chunks of the evaluation, and the late dual stages the direct
   summation. */

TEST_FUNCTION_START(mp_real_const_khinchin_afe, state)
{
    slong iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        slong n, a, maxdual, acc;
        mp_real_t x, y;
        arb_t A, B;

        n = 1 + n_randint(state, (iter % 10 == 9) ? 400 : 40);
        a = 1 + n_randint(state, 4);
        /* the extreme split points (many k-terms resp. dual terms) */
        if (iter % 10 == 4)
        {
            n = 1 + n_randint(state, 8);
            a = n_randint(state, 2) ? 7 : 5 + n_randint(state, 2);
        }
        maxdual = n_randint(state, 4);
        flint_set_num_threads(1 + n_randint(state, 3));

        mp_real_init(x);
        mp_real_init(y);
        arb_init(A);
        arb_init(B);

        _mp_real_const_khinchin_afe(x, n, a, maxdual);
        mp_real_get_arb(A, x);
        _mp_real_const_khinchin_tune(y, n + 1 + n_randint(state, 3), 0, 0, 0);
        mp_real_get_arb(B, y);

        acc = arb_rel_accuracy_bits(A);
        if (!arb_overlaps(A, B) || acc < FLINT_BITS * n - 16)
        {
            flint_printf("FAIL: n = %wd, a = %wd, maxdual = %wd, acc = %wd\n", n, a, maxdual, acc);
            flint_printf("A = "); arb_printn(A, 50, 0); flint_printf("\n");
            flint_printf("B = "); arb_printn(B, 50, 0); flint_printf("\n");
            flint_abort();
        }

        mp_real_clear(x);
        mp_real_clear(y);
        arb_clear(A);
        arb_clear(B);
    }

    flint_set_num_threads(1);

    TEST_FUNCTION_END(state);
}
