/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sin_cos_generic, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        fmpq_t q;
        mpfr_t t, u;
        slong prec0, prec;

        prec0 = 400;
        if (iter % 100 == 0)
            prec0 = 8000;

        prec = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        fmpq_init(q);
        mpfr_init2(t, prec + 200);
        mpfr_init2(u, prec + 200);

        arb_randtest(a, state, 1 + n_randint(state, prec0), 6);
        arb_randtest(b, state, 1 + n_randint(state, prec0), 6);
        arb_randtest(c, state, 1 + n_randint(state, prec0), 6);
        arb_get_rand_fmpq(q, state, a, 1 + n_randint(state, prec0));

        fmpq_get_mpfr(t, q, MPFR_RNDN);
        mpfr_sin_cos(t, u, t, MPFR_RNDN);

        arb_sin_cos_generic(b, c, a, prec);

        if (!arb_contains_mpfr(b, t))
        {
            flint_printf("FAIL: containment (sin)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        if (!arb_contains_mpfr(c, u))
        {
            flint_printf("FAIL: containment (cos)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        fmpq_clear(q);
        mpfr_clear(t);
        mpfr_clear(u);
    }

    /* check large arguments */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d, e;
        slong prec0, prec1, prec2, prec3;

        if (iter % 10 == 0)
            prec0 = 10000;
        else
            prec0 = 1000;

        prec1 = 2 + n_randint(state, prec0);
        prec2 = 2 + n_randint(state, prec0);

        if (iter % 10 == 0)
            prec3 = 50000;
        else
            prec3 = 100;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(e);

        arb_randtest_special(a, state, 1 + n_randint(state, prec0), prec3);
        arb_randtest_special(b, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(d, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(e, state, 1 + n_randint(state, prec0), 100);

        arb_sin_cos_generic(b, c, a, prec1);
        arb_sin_cos_generic(d, e, a, prec2);

        if (!arb_overlaps(b, d) || !arb_overlaps(c, e))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_printf("e = "); arb_print(e); flint_printf("\n\n");
            flint_abort();
        }

        /* check sin(a)^2 + cos(a)^2 = 1 */
        arb_mul(d, b, b, prec1);
        arb_mul(e, c, c, prec1);
        arb_add(d, d, e, prec1);
        arb_sub_ui(d, d, 1, prec1);

        if (!arb_contains_zero(d))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        arb_clear(e);
    }

    /* check accuracy */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, s, c;
        slong prec0, prec;
        mag_t allow;

        if (iter % 10 == 0)
            prec0 = 10000;
        else
            prec0 = 400;

        prec = 2 + n_randint(state, prec0);

        arb_init(a);
        arb_init(s);
        arb_init(c);
        mag_init(allow);

        arb_randtest_special(a, state, 1 + n_randint(state, prec0), 12);
        arb_randtest_special(s, state, 1 + n_randint(state, prec0), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, prec0), 100);

        if (n_randint(state, 2))
        {
            arb_sin_cos_generic(s, c, a, prec);
        }
        else
        {
            arb_sin_cos_generic(s, NULL, a, prec);
            arb_sin_cos_generic(NULL, c, a, prec);
        }

        if (!arb_is_finite(a))
        {
            mag_inf(allow);
        }
        else
        {
            mag_set_ui_2exp_si(allow, 1, -prec + 1);
            mag_max(allow, allow, arb_radref(a));
            mag_mul_2exp_si(allow, allow, 1);
        }

        if (mag_cmp(arb_radref(s), allow) > 0 ||
            mag_cmp(arb_radref(c), allow) > 0)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("a = "); arb_printn(a, 500, 0); flint_printf("\n\n");
            flint_printf("s = "); arb_printn(s, 500, 0); flint_printf("\n\n");
            flint_printf("c = "); arb_printn(c, 500, 0); flint_printf("\n\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("allow = "); mag_printd(allow, 5); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(s);
        arb_clear(c);
        mag_clear(allow);
    }

    TEST_FUNCTION_END(state);
}
