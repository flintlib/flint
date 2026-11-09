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

/* arb_sin_cos_pi and arb_tan_pi at dyadic x = m / 2^k with a small
   denominator (the algebraic methods of mp_real behind the fixed-point
   kernels: closed forms, half angles, the Chebyshev iteration), at
   precisions where they apply: overlap with sin and cos of pi x by the
   general functions at a higher precision, and full accuracy */
TEST_FUNCTION_START(arb_sin_cos_pi_dyadic, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        arb_t x, s, c, t, rs, rc, rt;
        slong prec, wp, k;

        arb_init(x);
        arb_init(s);
        arb_init(c);
        arb_init(t);
        arb_init(rs);
        arb_init(rc);
        arb_init(rt);

        k = 1 + n_randint(state, n_randint(state, 4) ? 8 : 24);
        arb_set_si(x, 2 * (slong) n_randint(state, UWORD(1) << FLINT_MIN(k + 3, 30)) + 1);
        if (n_randint(state, 2))
            arb_neg(x, x);
        arb_mul_2exp_si(x, x, -k);

        prec = 2 + n_randint(state, (iter % 4 == 0) ? 40000 : 6000);
        wp = 2 * prec + 100;

        arb_sin_cos_pi(s, c, x, prec);
        arb_tan_pi(t, x, prec);

        arb_const_pi(rt, wp);
        arb_mul(rt, rt, x, wp);
        arb_sin_cos(rs, rc, rt, wp);
        arb_div(rt, rs, rc, wp);

        if (!arb_overlaps(s, rs) || !arb_overlaps(c, rc) || !arb_overlaps(t, rt))
        {
            flint_printf("FAIL: overlap, k = %wd, prec = %wd\n", k, prec);
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n");
            flint_printf("s = "); arb_printd(s, 30); flint_printf("  "); arb_printd(rs, 30); flint_printf("\n");
            flint_printf("c = "); arb_printd(c, 30); flint_printf("  "); arb_printd(rc, 30); flint_printf("\n");
            flint_printf("t = "); arb_printd(t, 30); flint_printf("  "); arb_printd(rt, 30); flint_printf("\n");
            flint_abort();
        }

        if ((!arb_is_exact(s) && arb_rel_accuracy_bits(s) < prec - 4)
            || (!arb_is_exact(c) && arb_rel_accuracy_bits(c) < prec - 4)
            || (arb_is_finite(t) && !arb_is_exact(t) && arb_rel_accuracy_bits(t) < prec - 4))
        {
            flint_printf("FAIL: accuracy, k = %wd, prec = %wd\n", k, prec);
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n");
            flint_printf("s = "); arb_printd(s, 30); flint_printf("\n");
            flint_printf("c = "); arb_printd(c, 30); flint_printf("\n");
            flint_printf("t = "); arb_printd(t, 30); flint_printf("\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(s);
        arb_clear(c);
        arb_clear(t);
        arb_clear(rs);
        arb_clear(rc);
        arb_clear(rt);
    }

    TEST_FUNCTION_END(state);
}
