/*
    Copyright (C) 2010, 2012 William Hart
    Copyright (C) 2011, 2012 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE2_FUNCTION_START(T, poly_evaluate, vec_fast, state)
{
    int i, result = 1;

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) P;
        TEMPLATE(T, struct) * x, * y, * z;
        slong j, n, npoints;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        npoints = n_randint(state, 10);
        n = n_randint(state, 10);

        TEMPLATE(T, poly_init)(P, ctx);
        x = _TEMPLATE(T, vec_init)(npoints, ctx);
        y = _TEMPLATE(T, vec_init)(npoints, ctx);
        z = _TEMPLATE(T, vec_init)(npoints, ctx);

        TEMPLATE(T, poly_randtest)(P, state, n, ctx);

        for (j = 0; j < npoints; j++)
            TEMPLATE(T, randtest)(x + j, state, ctx);

        TEMPLATE4(T, poly_evaluate, T, vec_iter)(y, P, x, npoints, ctx);
        TEMPLATE4(T, poly_evaluate, T, vec_fast)(z, P, x, npoints, ctx);

        result = _TEMPLATE(T, vec_equal)(y, z, npoints, ctx);

        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, ctx_print)(ctx);
            flint_printf("\nn=%wd, npoints=%wd\n\n", n, npoints);
            flint_printf("P: "); TEMPLATE(T, poly_print_pretty)(P, "x", ctx);
            flint_printf("\n\n");
            for (j = 0; j < npoints; j++)
                TEMPLATE(T, print_pretty)(y + j, ctx), flint_printf(" ");
            flint_printf("\n");
            for (j = 0; j < npoints; j++)
                TEMPLATE(T, print_pretty)(z + j, ctx), flint_printf(" ");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear)(P, ctx);
        _TEMPLATE(T, vec_clear)(x, npoints, ctx);
        _TEMPLATE(T, vec_clear)(y, npoints, ctx);
        _TEMPLATE(T, vec_clear)(z, npoints, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
