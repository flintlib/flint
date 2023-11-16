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

TEST_FUNCTION_START(di_fast_add, state)
{
    slong iter;

    for (iter = 0; iter < 1000000 * 0.1 * flint_test_multiplier(); iter++)
    {
        di_t x, y, z;
        arf_t a, b, t, za, zb;

        arf_init(a);
        arf_init(b);
        arf_init(t);
        arf_init(za);
        arf_init(zb);

        x = di_randtest(state);
        y = di_randtest(state);

        z = di_fast_add(x, y);

        DI_CHECK(z)

        arf_set_d(a, x.a);
        arf_set_d(t, y.a);
        arf_add(a, a, t, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_d(b, x.b);
        arf_set_d(t, y.b);
        arf_add(b, b, t, ARF_PREC_EXACT, ARF_RND_DOWN);

        arf_set_d(za, z.a);
        arf_set_d(zb, z.b);

        if (arf_cmp(a, za) < 0 || arf_cmp(b, zb) > 0)
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); di_print(x); printf("\n");
            flint_printf("y = "); di_print(y); printf("\n");
            flint_printf("z = "); di_print(z); printf("\n");
            flint_printf("a = "); arf_printd(a, 20); printf("\n");
            flint_printf("b = "); arf_printd(b, 20); printf("\n");
            flint_printf("za = "); arf_printd(za, 20); printf("\n");
            flint_printf("zb = "); arf_printd(zb, 20); printf("\n");
            flint_abort();
        }

        arf_clear(a);
        arf_clear(b);
        arf_clear(t);
        arf_clear(za);
        arf_clear(zb);
    }

    TEST_FUNCTION_END(state);
}
