/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

static void
arb_supremum(arf_t res, const arb_t x)
{
    if (arf_is_nan(arb_midref(x)))
    {
        arf_nan(res);
    }
    else if (mag_is_inf(arb_radref(x)))
    {
        arf_pos_inf(res);
    }
    else
    {
        arf_set_mag(res, arb_radref(x));
        arf_add(res, res, arb_midref(x), ARF_PREC_EXACT, ARF_RND_CEIL);
    }
}

static void
arb_infimum(arf_t res, const arb_t x)
{
    if (arf_is_nan(arb_midref(x)))
    {
        arf_nan(res);
    }
    else if (mag_is_inf(arb_radref(x)))
    {
        arf_neg_inf(res);
    }
    else
    {
        arf_set_mag(res, arb_radref(x));
        arf_sub(res, arb_midref(x), res, ARF_PREC_EXACT, ARF_RND_FLOOR);
    }
}

int
arb_richcmp(const arb_t x, const arb_t y, int op)
{
    switch (op)
    {
        case 0: return arb_eq(x, y);
        case 1: return arb_ne(x, y);
        case 2: return arb_le(x, y);
        case 3: return arb_lt(x, y);
        case 4: return arb_ge(x, y);
        default: return arb_gt(x, y);
    }
}

int
arb_richcmp_fallback(const arb_t x, const arb_t y, int op)
{
    arf_t xa, xb, ya, yb;
    int res;

    arf_init(xa);
    arf_init(xb);
    arf_init(ya);
    arf_init(yb);

    arb_infimum(xa, x);
    arb_supremum(xb, x);
    arb_infimum(ya, y);
    arb_supremum(yb, y);

    if (arf_is_nan(xa) || arf_is_nan(ya))
    {
        res = 0;
    }
    else
    {
        if (op == 0) /* eq */
        {
            res = arf_equal(xa, xb) && arf_equal(ya, yb) && arf_equal(xa, ya);
        }
        else if (op == 1) /* ne */
        {
            res = (arf_cmp(yb, xa) < 0) || (arf_cmp(xb, ya) < 0);
        }
        else if (op == 2) /* le */
        {
            res = (arf_cmp(xb, ya) <= 0);
        }
        else if (op == 3) /* lt */
        {
            res = (arf_cmp(xb, ya) < 0);
        }
        else if (op == 4) /* ge */
        {
            res = (arf_cmp(xa, yb) >= 0);
        }
        else              /* gt */
        {
            res = (arf_cmp(xa, yb) > 0);
        }
    }

    arf_clear(xa);
    arf_clear(xb);
    arf_clear(ya);
    arf_clear(yb);

    return res;
}

TEST_FUNCTION_START(arb_richcmp, state)
{
    slong iter;

    for (iter = 0; iter < 1000000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        int op, res1, res2;

        arb_init(a);
        arb_init(b);

        arb_randtest_special(a, state, 1 + n_randint(state, 100), 5);
        arb_randtest_special(b, state, 1 + n_randint(state, 100), 5);

        op = n_randint(state, 6);

        res1 = arb_richcmp(a, b, op);
        res2 = arb_richcmp_fallback(a, b, op);

        if (res1 != res2)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("a = "); arb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 30); flint_printf("\n\n");
            flint_printf("op = %d\n\n", op);
            flint_printf("res1 (cmp)      = %d\n\n", res1);
            flint_printf("res2 (fallback) = %d\n\n", res2);
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
