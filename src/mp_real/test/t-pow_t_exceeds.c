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
#include "mp_real/impl.h"

/* _mp_real_pow_t_exceeds: when it claims |y log |x|| > 2^lim for |y| >=
   2^(ye - 1), the logarithm confirms it; mantissas with long leading
   runs of ones (near 1 from below) and zeros (from above) included */
TEST_FUNCTION_START(mp_real_pow_t_exceeds, state)
{
    slong iter, hits = 0;

    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        slong xn = 1 + n_randint(state, 3), xe, ye, lim, i, run;
        nn_ptr d = flint_malloc(xn * sizeof(ulong));
        arf_t x;
        arb_t t;
        mag_t lb, bound;

        flint_mpn_rrandom(d, state, xn);
        xe = (slong) n_randint(state, 9) - 4;
        if (n_randint(state, 2))
            xe = n_randint(state, 2);       /* 0 or 1: near 1 */

        /* a leading run below the top bit: ones (xe = 0) or zeros (xe = 1) */
        run = n_randint(state, FLINT_BITS * xn);
        for (i = 1; i <= run; i++)
        {
            slong b = FLINT_BITS * xn - 1 - i;
            if (xe == 0)
                d[b / FLINT_BITS] |= UWORD(1) << (b % FLINT_BITS);
            else
                d[b / FLINT_BITS] &= ~(UWORD(1) << (b % FLINT_BITS));
        }
        d[xn - 1] |= UWORD(1) << (FLINT_BITS - 1);

        ye = n_randint(state, 400);
        lim = ye - (slong) n_randint(state, 400) + 5;

        arf_init(x);
        arb_init(t);
        mag_init(lb);
        mag_init(bound);
        arf_set_mpn(x, d, xn, 0);
        arf_mul_2exp_si(x, x, xe - FLINT_BITS * xn);

        if (!arf_is_one(x) && _mp_real_pow_t_exceeds(d, xn, xe, ye, lim))
        {
            hits++;
            arb_set_arf(t, x);
            arb_log(t, t, FLINT_BITS * xn + 128);
            arb_mul_2exp_si(t, t, ye - 1);
            arb_get_mag_lower(lb, t);
            mag_one(bound);
            mag_mul_2exp_si(bound, bound, lim);
            if (mag_cmp(lb, bound) <= 0)
            {
                flint_printf("FAIL: xe = %wd, ye = %wd, lim = %wd\n", xe, ye, lim);
                flint_printf("x = "); arf_printd(x, 30); flint_printf("\n");
                flint_printf("t = "); arb_printd(t, 30); flint_printf("\n");
                flint_abort();
            }
        }

        arf_clear(x);
        arb_clear(t);
        mag_clear(lb);
        mag_clear(bound);
        flint_free(d);
    }

    if (hits == 0)
    {
        flint_printf("FAIL: no argument decided\n");
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
