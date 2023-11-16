/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb_modular.h"

TEST_FUNCTION_START(acb_modular_psl2z_mul, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        psl2z_t f, g, h, u, v;

        psl2z_init(f);
        psl2z_init(g);
        psl2z_init(h);
        psl2z_init(u);
        psl2z_init(v);

        psl2z_randtest(f, state, n_randint(state, 100));
        psl2z_randtest(g, state, n_randint(state, 100));
        psl2z_randtest(h, state, n_randint(state, 100));
        psl2z_randtest(u, state, n_randint(state, 100));
        psl2z_randtest(v, state, n_randint(state, 100));

        /* test (f*g)*h = f*(g*h) */

        psl2z_mul(u, f, g);
        psl2z_mul(u, u, h);

        psl2z_mul(v, g, h);
        psl2z_mul(v, f, v);

        if (!psl2z_equal(u, v) || !psl2z_is_correct(u) || !psl2z_is_correct(v))
        {
            flint_printf("FAIL\n");
            flint_printf("f = "); psl2z_print(f); flint_printf("\n");
            flint_printf("g = "); psl2z_print(g); flint_printf("\n");
            flint_printf("h = "); psl2z_print(h); flint_printf("\n");
            flint_printf("u = "); psl2z_print(u); flint_printf("\n");
            flint_printf("v = "); psl2z_print(v); flint_printf("\n");
            flint_abort();
        }

        psl2z_clear(f);
        psl2z_clear(g);
        psl2z_clear(h);
        psl2z_clear(u);
        psl2z_clear(v);
    }

    TEST_FUNCTION_END(state);
}
