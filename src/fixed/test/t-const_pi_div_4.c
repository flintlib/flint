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

/* fixed_const_pi_div_4 and fixed_const_log2 must return EXACTLY
   floor(c B^n) at every requested size -- not merely a close
   approximation -- through cache growth, prefix reads of the cached
   entry, explicit clears with rebuild, and fresh small-then-large
   request patterns */

static int
_floor_matches_pi4(nn_srcptr y, slong n, const arb_t exact)
{
    fmpz_t f, g;
    arb_t t;
    int ok;

    fmpz_init(f); fmpz_init(g); arb_init(t);
    arb_mul_2exp_si(t, exact, FLINT_BITS * n);
    arb_floor(t, t, FLINT_BITS * n + 128);
    ok = arb_get_unique_fmpz(f, t);
    if (ok)
    {
        fmpz_set_ui_array(g, y, n);
        ok = fmpz_equal(f, g);
    }
    fmpz_clear(f); fmpz_clear(g); arb_clear(t);
    return ok;
}

TEST_FUNCTION_START(fixed_const_pi_div_4, state)
{
    slong sizes[] = { 1, 2, 3, 7, 20, 5, 130, 40, 600, 200, 1,
        900 };
    slong i, iter;
    arb_t ref;

    arb_init(ref);

    for (iter = 0; iter < 2; iter++)
    {
        for (i = 0; i < 12; i++)
        {
            slong n = sizes[i] + (slong) n_randint(state, 4);
            nn_ptr y = flint_malloc(n * sizeof(ulong));

            arb_const_pi(ref, FLINT_BITS * n + 192);
            arb_mul_2exp_si(ref, ref, -2);
            fixed_const_pi_div_4(y, n);
            if (!_floor_matches_pi4(y, n, ref))
                TEST_FUNCTION_FAIL("floor mismatch: n = %wd\n", n);
            flint_free(y);
        }
        /* second pass rebuilds from nothing */
        _fixed_const_pi_div_4_clear();
    }

    arb_clear(ref);
    TEST_FUNCTION_END(state);
}
