/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fixed.h"
#include "arb.h"

TEST_FUNCTION_START(fixed_tab_bsplit, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        /* every fourth iteration goes big to reach the direct-log
           regime (prec >= 6000, i >= 30); iter = 0 pins one case
           inside it, since iter % 4 == 0 and iter % 2 would
           otherwise anti-correlate the size with the function
           choice and starve that branch */
        slong n = (iter == 0) ? 100
            : 1 + n_randint(state, (iter % 4 == 0) ? 160 : 12);
        slong i = (iter == 0) ? 40
            : 1 + n_randint(state, (iter % 3 == 0) ? 700 : 40);
        int is_log = (iter == 0) ? 1 : (int) n_randint(state, 2);
        nn_ptr got, ref;
        arb_t v, t;
        arf_t lb;
        fmpz_t f;
        slong j;
        int ok;

        got = flint_malloc(n * sizeof(ulong));
        ref = flint_malloc(n * sizeof(ulong));
        arb_init(v); arb_init(t); arf_init(lb); fmpz_init(f);

        if (is_log)
            fixed_log1p_2mexp_ui_bs(got, (ulong) i, n);
        else
            fixed_atan_2mexp_ui_bs(got, (ulong) i, n);

        arb_one(v);
        arb_mul_2exp_si(v, v, -i);
        if (is_log)
            arb_log1p(t, v, FLINT_BITS * n + 64);
        else
            arb_atan(t, v, FLINT_BITS * n + 64);
        arb_get_lbound_arf(lb, t, FLINT_BITS * n + 64);
        arf_mul_2exp_si(lb, lb, FLINT_BITS * n);
        arf_get_fmpz(f, lb, ARF_RND_FLOOR);
        fmpz_get_ui_array(ref, n, f);

        /* one-sided: the exact floor or up to 3 ulps below, never
           above */
        ok = 0;
        for (j = 0; j <= 3 && !ok; j++)
        {
            ok = (mpn_cmp(got, ref, n) == 0);
            mpn_sub_1(ref, ref, n, 1);
        }

        if (!ok)
            TEST_FUNCTION_FAIL("%s: i = %wd, n = %wd\n",
                is_log ? "log1p" : "atan", i, n);

        arb_clear(v); arb_clear(t); arf_clear(lb); fmpz_clear(f);
        flint_free(got); flint_free(ref);
    }

    TEST_FUNCTION_END(state);
}
