/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "elem_check.h"

/* x^y against arb_pow on the exact arguments: random arguments, small
   integer and dyadic exponents, perfect powers and powers of two */
TEST_FUNCTION_START(arf_pow, state)
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

        switch (n_randint(state, 8))
        {
            case 0:
                /* a small integer or half-integer exponent */
                arf_set_si(y, (slong) n_randint(state, 41) - 20);
                if (n_randint(state, 2))
                    arf_mul_2exp_si(y, y, -(slong) n_randint(state, 3));
                break;
            case 1:
                /* a perfect power of a short base */
                {
                    arf_t t;
                    ulong k = 2 + n_randint(state, 7);
                    arf_init(t);
                    arf_urandom(t, state, 1 + n_randint(state, 40), ARF_RND_DOWN);
                    arf_mul_2exp_si(t, t, (slong) n_randint(state, 20) - 10);
                    if (arf_is_zero(t))
                        arf_one(t);
                    arf_one(x);
                    while (k-- > 0)
                        arf_mul(x, x, t, ARF_PREC_EXACT, ARF_RND_DOWN);
                    arf_clear(t);
                    arf_set_si(y, 1);
                    arf_mul_2exp_si(y, y, -(slong) (1 + n_randint(state, 3)));
                    if (n_randint(state, 2))
                        arf_mul_si(y, y, (slong) n_randint(state, 21) - 10, ARF_PREC_EXACT, ARF_RND_DOWN);
                }
                break;
            case 2:
                /* a power of two */
                arf_one(x);
                arf_mul_2exp_si(x, x, (slong) n_randint(state, 2000) - 1000);
                if (n_randint(state, 2))
                    arf_neg(x, x);
                break;
            case 3:
                /* x close to 1 */
                arf_set_ui_2exp_si(x, 1, -(slong) n_randint(state, 2 * prec + 10));
                if (n_randint(state, 2))
                    arf_neg(x, x);
                arf_add_ui(x, x, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
                break;
            default:
                break;
        }

        inexact = arf_pow(r, x, y, prec, rnd);

        /* aliasing */
        arf_set(s, x);
        inexact2 = arf_pow(s, s, y, prec, rnd);
        if (!arf_equal(r, s) || (inexact != inexact2 && !arf_is_nan(r)))
        {
            flint_printf("FAIL (aliasing x): prec = %wd, rnd = %d\n", prec, rnd);
            flint_printf("x = "); arf_debug(x); flint_printf("\n");
            flint_printf("y = "); arf_debug(y); flint_printf("\n");
            flint_printf("r = "); arf_debug(r); flint_printf("\n");
            flint_printf("s = "); arf_debug(s); flint_printf("\n");
            flint_abort();
        }
        arf_set(s, y);
        inexact2 = arf_pow(s, x, s, prec, rnd);
        if (!arf_equal(r, s) || (inexact != inexact2 && !arf_is_nan(r)))
        {
            flint_printf("FAIL (aliasing y)\n");
            flint_abort();
        }

        if (arf_is_finite(x) && arf_is_finite(y))
        {
            arb_set_arf(bx, x);
            arb_set_arf(by, y);
            arb_pow(b, bx, by, prec + 64);

            /* (arb_pow of a negative base rounds the midpoint of a power
               before checking the sign: a nan from arb is skipped) */
            if (arf_is_nan(r) && !arf_is_nan(arb_midref(b)) && !(arf_sgn(x) < 0 && arf_is_int(y)))
            {
                /* x^y beyond the exponent range: |y log x| must be huge */
                arb_t l;
                arb_init(l);
                arb_abs(l, bx);
                arb_log(l, l, 64);
                arb_mul(l, l, by, 64);
                if (arb_rel_accuracy_bits(l) > 10 && arf_cmpabs_2exp_si(arb_midref(l), FLINT_MAX(128, 2 * prec) - 4) < 0)
                {
                    flint_printf("FAIL (nan): prec = %wd\n", prec);
                    flint_printf("x = "); arf_debug(x); flint_printf("\n");
                    flint_printf("y = "); arf_debug(y); flint_printf("\n");
                    flint_abort();
                }
                arb_clear(l);
            }
            else if (!arf_is_nan(r))
                arf_elem_verify2("arf_pow", x, y, r, inexact, b, prec, rnd, WORD_MAX);
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
