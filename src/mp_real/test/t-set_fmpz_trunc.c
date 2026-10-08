/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "arb.h"
#include "mp_real.h"
#include "ball_helpers.h"

/* mp_real_set_fmpz is exact, with the canonical exact representation
   (top limb nonzero, low zero limbs stripped); mp_real_set_trunc keeps at
   most n limbs, contains its input, has a radius of at most 2 ulps of the
   new bottom limb (at most 1 if the input was exact), is exact when the
   input fits, and may alias its input. */

TEST_FUNCTION_START(mp_real_set_fmpz_trunc, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_t f, g;
        mp_real_t x;
        arb_t a;

        fmpz_init(f);
        fmpz_init(g);
        mp_real_init(x);
        arb_init(a);

        fmpz_randtest(f, state, 2 + n_randint(state, 600));
        if (n_randint(state, 4) == 0)
            fmpz_mul_2exp(f, f, n_randint(state, 500));
        if (n_randint(state, 8) == 0)
            fmpz_set_si(f, n_randint(state, 2) ? COEFF_MAX : COEFF_MIN);
        if (n_randint(state, 4) == 0)
            mp_real_randtest(x, state, 10, 1);      /* overwrite garbage */

        mp_real_set_fmpz(x, f);
        mp_real_get_arb(a, x);
        if (!arb_is_exact(a) || !arb_get_unique_fmpz(g, a) || !fmpz_equal(f, g)
            || (x->size > 0 && (x->d[x->size - 1] == 0 || x->d[0] == 0))
            || (fmpz_is_zero(f) != (x->size == 0)))
        {
            flint_printf("FAIL: set_fmpz\n");
            fmpz_print(f); flint_printf("\n");
            flint_abort();
        }

        fmpz_clear(f);
        fmpz_clear(g);
        mp_real_clear(x);
        arb_clear(a);
    }

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        mp_real_t x, y;
        arb_t a, b;
        slong n = 1 + n_randint(state, 12);
        int alias = n_randint(state, 2);

        mp_real_init(x);
        mp_real_init(y);
        arb_init(a);
        arb_init(b);

        mp_real_randtest(x, state, 16, n_randint(state, 2));
        mp_real_get_arb(a, x);

        if (alias)
        {
            mp_real_set(y, x);
            mp_real_set_trunc(y, y, n);
        }
        else
            mp_real_set_trunc(y, x, n);
        mp_real_get_arb(b, y);

        if (!arb_contains(b, a) || y->size > n
            || (x->size <= n && !arb_equal(a, b))
            || (x->size > n && y->err > (x->err != 0) + 1))
        {
            flint_printf("FAIL: set_trunc, n = %wd, alias = %d\n", n, alias);
            arb_printd(a, 50); flint_printf("\n");
            arb_printd(b, 50); flint_printf("\n");
            flint_abort();
        }

        mp_real_clear(x);
        mp_real_clear(y);
        arb_clear(a);
        arb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
