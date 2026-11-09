/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "test_helpers.h"
#include "arb.h"

/* arb_tan against sin / cos at a higher precision and against MPFR for
   exact arguments (containment, and the relative accuracy of about
   prec bits for exact arguments away from the poles), with aliasing */
TEST_FUNCTION_START(arb_tan, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, a, b, s, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 300);
        if (n_randint(state, 10) == 0)
            prec1 = 2 + n_randint(state, 2000);
        prec2 = prec1 + 30;

        arb_init(x);
        arb_init(a);
        arb_init(b);
        arb_init(s);
        arb_init(c);

        arb_randtest(x, state, 1 + n_randint(state, 300), 2 + n_randint(state, 10));
        if (n_randint(state, 3) == 0)
            mag_zero(arb_radref(x));
        if (n_randint(state, 4) == 0)
            arb_mul_2exp_si(x, x, -(slong) n_randint(state, 2 * prec1));

        arb_tan(a, x, prec1);
        arb_tan(b, x, prec2);

        if (!arb_overlaps(a, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_sin_cos(s, c, x, prec2);
        arb_div(c, s, c, prec2);

        if (!arb_overlaps(a, c))
        {
            flint_printf("FAIL: sin / cos\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_is_exact(x) && !arf_is_zero(arb_midref(x))
            && arf_cmpabs_2exp_si(arb_midref(x), 20) < 0
            && arf_cmpabs_2exp_si(arb_midref(x), -prec1 - 100) > 0)
        {
            /* MPFR at prec1 + 64 bits: contained (its ulp as the error),
               and the accuracy of arb's result unless near a pole */
            mpfr_t m;
            arf_t t;
            mpfr_init2(m, FLINT_MAX(prec1 + 64, arf_bits(arb_midref(x))));
            arf_init(t);
            arf_get_mpfr(m, arb_midref(x), MPFR_RNDN);
            mpfr_tan(m, m, MPFR_RNDN);
            arf_set_mpfr(t, m);
            arb_set_arf(c, t);
            arf_mag_set_ulp(arb_radref(c), t, prec1 + 60);

            if (!arb_overlaps(a, c) || (arb_rel_accuracy_bits(a) < prec1 - 4
                && !(arb_is_finite(a) && arf_cmpabs_2exp_si(t, prec1 / 2) > 0)))
            {
                flint_printf("FAIL: mpfr, prec = %wd, acc = %wd\n\n", prec1, arb_rel_accuracy_bits(a));
                flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                flint_printf("a = "); arb_printd(a, 30); flint_printf("\n\n");
                flint_printf("c = "); arb_printd(c, 30); flint_printf("\n\n");
                flint_abort();
            }
            mpfr_clear(m);
            arf_clear(t);
        }

        arb_tan(x, x, prec1);

        if (!arb_overlaps(a, x))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(a);
        arb_clear(b);
        arb_clear(s);
        arb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
