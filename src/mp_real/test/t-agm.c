/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"
#include "ball_helpers.h"

/* the AGM against arb_agm at the midpoints, for exact and inexact
   inputs, equal inputs, nearly equal inputs, far-apart inputs and
   every order of the finish */
static void
test_agm(flint_rand_t state, slong iters)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        slong n = 1 + n_randint(state, (iter % 10 == 0) ? 300 : 30);
        int m = n_randint(state, 3) ? 0 : 2 + (int) n_randint(state, 15);
        mp_real_t x, y, r;
        arb_t xa, ya, ra, rr;
        int w;

        mp_real_init(x); mp_real_init(y); mp_real_init(r);
        arb_init(xa); arb_init(ya); arb_init(ra); arb_init(rr);

        for (w = 0; w < 2; w++)
        {
            mp_real_struct * v = w ? y : x;
            slong sz = 1 + n_randint(state, n + 2);
            nn_ptr t = flint_malloc(sz * sizeof(ulong));
            flint_mpn_urandomb(t, state, FLINT_BITS * sz);
            t[sz - 1] |= 1;
            if (n_randint(state, 5) == 0)
                mp_real_set_ui(v, 1 + n_randint(state, 10));
            else
                _mp_real_set_mpn_2exp(v, t, sz, -FLINT_BITS * sz
                    + (slong) n_randint(state, 400) - 200);
            if (n_randint(state, 3) == 0)
                _mp_real_add_error_ulps(v, 1 + n_randint(state, 100));
            flint_free(t);
        }
        if (n_randint(state, 6) == 0)
            mp_real_set(y, x);
        if (n_randint(state, 6) == 0)
        {
            mp_real_set(y, x);
            mp_real_add_error_2exp_si(y, mp_real_abs_bound_lt_2exp_si(x)
                - FLINT_BITS * (slong) n_randint(state, n + 3));
        }

        _mp_real_agm_order(r, x, y, n, m);
        mp_real_get_arb(xa, x);
        mp_real_get_arb(ya, y);
        mp_real_get_arb(rr, r);
        arb_get_mid_arb(xa, xa);
        arb_get_mid_arb(ya, ya);
        arb_agm(ra, xa, ya, FLINT_BITS * n + 128);
        if (!arb_overlaps(ra, rr))
        {
            flint_printf("FAIL: agm (iter %wd, n = %wd, m = %d)\n", iter, n, m);
            mp_real_print(x); mp_real_print(y); mp_real_print(r);
            flint_abort();
        }
        if (x->err == 0 && y->err == 0
            && arb_rel_accuracy_bits(rr) < FLINT_BITS * n - 16)
        {
            flint_printf("FAIL: agm accuracy (iter %wd, n = %wd, m = %d): %wd bits\n",
                iter, n, m, arb_rel_accuracy_bits(rr));
            mp_real_print(x); mp_real_print(y); mp_real_print(r);
            flint_abort();
        }

        mp_real_clear(x); mp_real_clear(y); mp_real_clear(r);
        arb_clear(xa); arb_clear(ya); arb_clear(ra); arb_clear(rr);
    }
}

TEST_FUNCTION_START(mp_real_agm, state)
{
    test_agm(state, 300 + 300 * flint_test_multiplier());

    TEST_FUNCTION_END(state);
}
