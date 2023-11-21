/*
    Copyright (C) 2012-2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_add_error, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        arf_t m, r;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arf_init(m);
        arf_init(r);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 10);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 10);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 10);
        arf_randtest_special(m, state, 1 + n_randint(state, 2000), 10);
        arf_randtest_special(r, state, 1 + n_randint(state, 2000), 10);

        /* c = a plus error bounds */
        arb_set(c, a);
        arf_set(arb_midref(b), m);
        arf_get_mag(arb_radref(b), r);
        arb_add_error(c, b);

        /* b = a + random point */
        arb_set(b, a);

        if (n_randint(state, 2))
            arf_add(arb_midref(b), arb_midref(b), m, ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_sub(arb_midref(b), arb_midref(b), m, ARF_PREC_EXACT, ARF_RND_DOWN);

        if (n_randint(state, 2))
            arf_add(arb_midref(b), arb_midref(b), r, ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_sub(arb_midref(b), arb_midref(b), r, ARF_PREC_EXACT, ARF_RND_DOWN);

        /* should this be done differently? */
        if (arf_is_nan(arb_midref(b)))
            arf_zero(arb_midref(b));

        if (!arb_contains(c, b))
        {
            flint_printf("FAIL (arb_add_error)\n\n");
            flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); arb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("c = "); arb_printn(c, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arf_clear(m);
        arf_clear(r);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        arf_t m;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arf_init(m);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 10);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 10);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 10);
        arf_randtest_special(m, state, 1 + n_randint(state, 2000), 10);

        /* c = a plus error bounds */
        arb_set(c, a);
        arb_add_error_arf(c, m);

        /* b = a + random point */
        arb_set(b, a);

        if (n_randint(state, 2))
            arf_add(arb_midref(b), arb_midref(b), m, ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_sub(arb_midref(b), arb_midref(b), m, ARF_PREC_EXACT, ARF_RND_DOWN);

        /* should this be done differently? */
        if (arf_is_nan(arb_midref(b)))
            arf_zero(arb_midref(b));

        if (!arb_contains(c, b))
        {
            flint_printf("FAIL (arb_add_error_arf)\n\n");
            flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); arb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("c = "); arb_printn(c, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arf_clear(m);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        arf_t t;
        mag_t r;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        mag_init(r);
        arf_init(t);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 10);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 10);
        mag_randtest(r, state, 10);

        /* c = a plus error bounds */
        arb_set(c, a);
        arb_add_error_mag(c, r);

        /* b = a + random point */
        arb_set(b, a);
        arf_set_mag(t, r);
        if (n_randint(state, 2))
            arf_add(arb_midref(b), arb_midref(b), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_sub(arb_midref(b), arb_midref(b), t, ARF_PREC_EXACT, ARF_RND_DOWN);

        /* should this be done differently? */
        if (arf_is_nan(arb_midref(b)))
            arf_zero(arb_midref(b));

        if (!arb_contains(c, b))
        {
            flint_printf("FAIL (arb_add_error_mag)\n\n");
            flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); arb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("c = "); arb_printn(c, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        mag_clear(r);
        arf_clear(t);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        arf_t t;
        slong e;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arf_init(t);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 10);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 10);
        e = n_randint(state, 10) - 10;

        /* c = a plus error bounds */
        arb_set(c, a);
        arb_add_error_2exp_si(c, e);

        /* b = a + random point */
        arb_set(b, a);
        arf_one(t);
        arf_mul_2exp_si(t, t, e);
        if (n_randint(state, 2))
            arf_add(arb_midref(b), arb_midref(b), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_sub(arb_midref(b), arb_midref(b), t, ARF_PREC_EXACT, ARF_RND_DOWN);

        /* should this be done differently? */
        if (arf_is_nan(arb_midref(b)))
            arf_zero(arb_midref(b));

        if (!arb_contains(c, b))
        {
            flint_printf("FAIL (arb_add_error_2exp_si)\n\n");
            flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); arb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("c = "); arb_printn(c, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arf_clear(t);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        arf_t t;
        fmpz_t e;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arf_init(t);
        fmpz_init(e);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 10);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 10);
        fmpz_randtest(e, state, 10);

        /* c = a plus error bounds */
        arb_set(c, a);
        arb_add_error_2exp_fmpz(c, e);

        /* b = a + random point */
        arb_set(b, a);
        arf_one(t);
        arf_mul_2exp_fmpz(t, t, e);
        if (n_randint(state, 2))
            arf_add(arb_midref(b), arb_midref(b), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_sub(arb_midref(b), arb_midref(b), t, ARF_PREC_EXACT, ARF_RND_DOWN);

        /* should this be done differently? */
        if (arf_is_nan(arb_midref(b)))
            arf_zero(arb_midref(b));

        if (!arb_contains(c, b))
        {
            flint_printf("FAIL (arb_add_error_2exp_fmpz)\n\n");
            flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); arb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("c = "); arb_printn(c, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arf_clear(t);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
