/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "elem_check.h"

/* atan2(y, x) against arb_atan2 on the exact arguments, with the
   magnitude ratio of the arguments spread over the tiny, normal and huge
   ranges */
TEST_FUNCTION_START(arf_atan2, state)
{
    slong iter;

    for (iter = 0; iter < 4000 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, r, s;
        arb_t b, bx, by;
        slong prec;
        arf_rnd_t rnd;
        int inexact, inexact2;

        arf_init(x);
        arf_init(y);
        arf_init(r);
        arf_init(s);
        arb_init(b);
        arb_init(bx);
        arb_init(by);

        prec = 2 + n_randint(state, 300);
        if (n_randint(state, 8) == 0)
            prec = 2 + n_randint(state, 3000);
        if (n_randint(state, 100) == 0)
            prec = 2 + n_randint(state, 12000);

        rnd = (arf_rnd_t) n_randint(state, 7);

        arf_elem_randtest_arg(x, state, prec);
        arf_elem_randtest_arg(y, state, prec);

        switch (n_randint(state, 6))
        {
            case 0:
                /* |y| = |x| */
                arf_set(y, x);
                if (n_randint(state, 2))
                    arf_neg(y, y);
                break;
            case 1:
                /* the same magnitude */
                arf_mul_2exp_fmpz(y, y, ARF_EXPREF(x));
                if (arf_is_normal(y))
                    arf_mul_2exp_fmpz(y, y, ARF_EXPREF(x)),
                    arf_mul_2exp_si(y, y, -1);
                break;
            case 2:
                if (n_randint(state, 2))
                    arf_set_si(x, n_randint(state, 3) - 1);
                else
                    arf_set_si(y, n_randint(state, 3) - 1);
                break;
            default:
                break;
        }

        inexact = arf_atan2(r, y, x, prec, rnd);

        /* aliasing */
        arf_set(s, y);
        inexact2 = arf_atan2(s, s, x, prec, rnd);
        if (!arf_equal(r, s) || (inexact != inexact2 && !arf_is_nan(r)))
        {
            flint_printf("FAIL (aliasing y)\n");
            flint_abort();
        }
        arf_set(s, x);
        inexact2 = arf_atan2(s, y, s, prec, rnd);
        if (!arf_equal(r, s) || (inexact != inexact2 && !arf_is_nan(r)))
        {
            flint_printf("FAIL (aliasing x)\n");
            flint_abort();
        }

        arb_set_arf(bx, x);
        arb_set_arf(by, y);
        arb_atan2(b, by, bx, prec + 64);

        if (arf_is_finite(x) && arf_is_finite(y))
            arf_elem_verify2("arf_atan2", x, y, r, inexact, b, prec, rnd, WORD_MAX);
        else if (arf_is_nan(x) || arf_is_nan(y))
        {
            if (!arf_is_nan(r))
            {
                flint_printf("FAIL (nan)\n");
                flint_abort();
            }
        }
        else
        {
            /* the conventions for infinities (arb's result is
               indeterminate): c pi/4 with the sign of y */
            int c, xs = arf_sgn(x), ys = arf_sgn(y);

            if (ys == 0)
                c = (xs >= 0) ? 0 : 4;
            else if (xs == 0 || (arf_is_inf(y) && arf_is_finite(x)))
                c = 2;
            else if (arf_is_finite(y))
                c = (xs > 0) ? 0 : 4;
            else
                c = (xs > 0) ? 1 : 3;

            arb_const_pi(b, prec + 64);
            arb_mul_2exp_si(b, b, -2);
            arb_mul_si(b, b, (ys < 0) ? -c : c, prec + 64);

            if (arf_is_nan(r))
            {
                flint_printf("FAIL (nan for an infinity)\n");
                flint_printf("x = "); arf_debug(x); flint_printf("\n");
                flint_printf("y = "); arf_debug(y); flint_printf("\n");
                flint_abort();
            }

            arf_elem_verify2("arf_atan2 (special)", x, y, r, inexact, b, prec, rnd, WORD_MAX);
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(r);
        arf_clear(s);
        arb_clear(b);
        arb_clear(bx);
        arb_clear(by);
    }

    TEST_FUNCTION_END(state);
}
