/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"
#include "double_interval.h"

/* Defined in t-fast_div.c and t-fast_mul.c */
#ifndef arf_min2
#define arf_min2 arf_min2
void
arf_min2(arf_t res, const arf_t x, const arf_t y)
{
    if (arf_is_nan(x))
        arf_set(res, y);
    else if (arf_is_nan(y))
        arf_set(res, x);
    else
        arf_min(res, x, y);
}
#endif

/* Defined in t-fast_div.c and t-fast_mul.c */
#ifndef arf_max2
#define arf_max2 arf_max2
void
arf_max2(arf_t res, const arf_t x, const arf_t y)
{
    if (arf_is_nan(x))
        arf_set(res, y);
    else if (arf_is_nan(y))
        arf_set(res, x);
    else
        arf_max(res, x, y);
}
#endif

TEST_FUNCTION_START(di_fast_div, state)
{
    slong iter;

    for (iter = 0; iter < 1000000 * 0.1 * flint_test_multiplier(); iter++)
    {
        di_t x, y, z;
        arf_t a, b, c, d, ra, rb, rc, rd, va, vb, za, zb;

        arf_init(a); arf_init(b); arf_init(c); arf_init(d);
        arf_init(ra); arf_init(rb); arf_init(rc); arf_init(rd);
        arf_init(va); arf_init(vb);
        arf_init(za);
        arf_init(zb);

        x = di_randtest(state);
        y = di_randtest(state);
        z = di_fast_div(x, y);

        DI_CHECK(z)

        arf_set_d(a, x.a);
        arf_set_d(b, x.b);
        arf_set_d(c, y.a);
        arf_set_d(d, y.b);

        arf_div(ra, a, c, 106, ARF_RND_DOWN);
        arf_div(rb, a, d, 106, ARF_RND_DOWN);
        arf_div(rc, b, c, 106, ARF_RND_DOWN);
        arf_div(rd, b, d, 106, ARF_RND_DOWN);

        arf_min2(va, ra, rb);
        arf_min2(va, va, rc);
        arf_min2(va, va, rd);
        arf_max2(vb, ra, rb);
        arf_max2(vb, vb, rc);
        arf_max2(vb, vb, rd);

        arf_set_d(za, z.a);
        arf_set_d(zb, z.b);

        if (arf_cmp(va, za) < 0 || arf_cmp(vb, zb) > 0)
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); di_print(x); printf("\n");
            flint_printf("y = "); di_print(y); printf("\n");
            flint_printf("z = "); di_print(z); printf("\n");
            flint_printf("a = "); arf_printd(a, 20); printf("\n");
            flint_printf("b = "); arf_printd(b, 20); printf("\n");
            flint_printf("c = "); arf_printd(c, 20); printf("\n");
            flint_printf("d = "); arf_printd(d, 20); printf("\n");
            flint_printf("ra = "); arf_printd(ra, 20); printf("\n");
            flint_printf("rb = "); arf_printd(rb, 20); printf("\n");
            flint_printf("rc = "); arf_printd(rc, 20); printf("\n");
            flint_printf("rd = "); arf_printd(rd, 20); printf("\n");
            flint_printf("va = "); arf_printd(va, 20); printf("\n");
            flint_printf("vb = "); arf_printd(vb, 20); printf("\n");
            flint_printf("za = "); arf_printd(za, 20); printf("\n");
            flint_printf("zb = "); arf_printd(zb, 20); printf("\n");
            flint_abort();
        }

        arf_clear(a); arf_clear(b); arf_clear(c); arf_clear(d);
        arf_clear(ra); arf_clear(rb); arf_clear(rc); arf_clear(rd);
        arf_clear(va); arf_clear(vb);
        arf_clear(za);
        arf_clear(zb);
    }

    TEST_FUNCTION_END(state);
}
