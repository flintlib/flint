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

static int
_floor_matches_euler(nn_srcptr y, slong n, const arb_t exact)
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

TEST_FUNCTION_START(fixed_const_euler, state)
{
    slong iter;
    arb_t ref, a;
    fball_t x;

    arb_init(ref);
    arb_init(a);
    fball_init(x);

    /* the balls, through each set of logarithms and the automatic
       choice, including the table below 54 limbs */
    for (iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        flint_set_num_threads(1 + n_randint(state, 4));
        slong n = 1 + n_randint(state, (iter % 8 == 0) ? 1000 : 150);
        int set = (int) n_randint(state, 5) - 1;

        if (set < 0)
            fball_const_euler(x, n);
        else
            _fball_const_euler(x, n, set);

        fball_get_arb(a, x);
        arb_const_euler(ref, FLINT_BITS * n + 64);

        if (!arb_overlaps(a, ref)
            || arb_rel_accuracy_bits(a) < FLINT_BITS * (n - 1) - 32)
        {
            flint_printf("FAIL: n = %wd, set = %d\n", n, set);
            arb_printd(a, 50); flint_printf("\n");
            arb_printd(ref, 50); flint_printf("\n");
            flint_abort();
        }
    }

    /* the verified floors, through cache growth, prefixes and clears */
    for (iter = 0; iter < 2; iter++)
    {
        flint_set_num_threads(1 + n_randint(state, 4));
        slong sizes[] = { 1, 2, 3, 7, 60, 5, 130, 40, 300, 100, 1 };
        slong i;

        for (i = 0; i < 11; i++)
        {
            slong n = sizes[i] + (slong) n_randint(state, 4);
            nn_ptr y = flint_malloc(n * sizeof(ulong));

            arb_const_euler(ref, FLINT_BITS * n + 192);
            fixed_const_euler(y, n);
            if (!_floor_matches_euler(y, n, ref))
                TEST_FUNCTION_FAIL("floor mismatch: n = %wd\n", n);
            flint_free(y);
        }
        _fixed_const_euler_clear();
    }

    arb_clear(ref);
    arb_clear(a);
    fball_clear(x);
    flint_set_num_threads(1);
    TEST_FUNCTION_END(state);
}
