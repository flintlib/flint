/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

static void
acb_inv_naive(acb_t z, const acb_t x, slong prec)
{
#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(z)
#define d acb_imagref(z)

    if (arb_is_zero(b))
    {
        arb_inv(c, a, prec);
        arb_zero(d);
    }
    else if (arb_is_zero(a))
    {
        arb_inv(d, b, prec);
        arb_neg(d, d);
        arb_zero(c);
    }
    else
    {
        arb_t t;
        arb_init(t);

        arb_mul(t, a, a, prec);
        arb_addmul(t, b, b, prec);

        arb_div(c, a, t, prec);
        arb_div(d, b, t, prec);

        arb_neg(d, d);

        arb_clear(t);
    }

#undef a
#undef b
#undef c
#undef d
}

TEST_FUNCTION_START(acb_inv, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, d, e, f;
        arf_t t;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        acb_init(e);
        acb_init(f);
        arf_init(t);

        prec = 2 + n_randint(state, 1000);

        acb_randtest_special(a, state, 1 + n_randint(state, 1000), 100);
        acb_randtest_special(b, state, 1 + n_randint(state, 1000), 100);

        acb_inv(b, a, prec);
        acb_inv_naive(c, a, prec);

        if (!acb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_set(c, a);
        acb_inv(c, c, prec);

        if (!acb_equal(b, c))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_randtest(a, state, 1 + n_randint(state, 1000), 10);
        acb_randtest(b, state, 1 + n_randint(state, 1000), 10);

        acb_zero(d);
        arf_set_mag(t, arb_radref(acb_realref(a)));
        if (n_randint(state, 2))
            arf_neg(t, t);
        arf_add(arb_midref(acb_realref(d)),
            arb_midref(acb_realref(a)), t, ARF_PREC_EXACT, ARF_RND_DOWN);

        arf_set_mag(t, arb_radref(acb_imagref(a)));
        if (n_randint(state, 2))
            arf_neg(t, t);
        arf_add(arb_midref(acb_imagref(d)),
            arb_midref(acb_imagref(a)), t, ARF_PREC_EXACT, ARF_RND_DOWN);

        acb_inv(b, a, 2 + n_randint(state, 1000));
        acb_inv(d, d, 2 + n_randint(state, 1000));

        if (!acb_overlaps(b, d))
        {
            flint_printf("FAIL: corner test\n\n");
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("d = "); acb_printd(d, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
        acb_clear(e);
        acb_clear(f);
        arf_clear(t);
    }

    TEST_FUNCTION_END(state);
}
