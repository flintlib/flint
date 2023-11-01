/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

void test_poly(
    TEMPLATE(T, poly_factor_t) roots,
    const TEMPLATE(T, poly_t) f,
    int want_mult,
    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, multiplicity;
    TEMPLATE(T, poly_t) q, qt, r;

    TEMPLATE(T, poly_init)(q, ctx);
    TEMPLATE(T, poly_init)(qt, ctx);
    TEMPLATE(T, poly_init)(r, ctx);
    TEMPLATE(T, poly_set)(q, f, ctx);

    TEMPLATE(T, poly_roots)(roots, f, want_mult, ctx);

    for (i = 0; i < roots->num; i++)
    {
        if (TEMPLATE(T, poly_degree)(roots->poly + i, ctx) != 1)
        {
            flint_printf("FAILED:\ncheck root is linear\n");
            fflush(stdout);
            flint_abort();
        }

        if (TEMPLATE(T, is_one)(roots->poly[i].coeffs + 1, ctx) == 0)
        {
            flint_printf("FAILED:\ncheck root is monic\n");
            fflush(stdout);
            flint_abort();
        }

        multiplicity = 0;
        while (TEMPLATE(T, poly_divrem)(qt, r, q, roots->poly + i, ctx),
               TEMPLATE(T, poly_is_zero)(r, ctx))
        {
            TEMPLATE(T, poly_swap)(q, qt, ctx);
            multiplicity++;
        }

        if (multiplicity <= 0)
        {
            flint_printf("FAILED:\ncheck root is a root\n");
            fflush(stdout);
            flint_abort();
        }

        if (roots->exp[i] != (want_mult ? multiplicity : 1))
        {
            flint_printf("FAILED:\ncheck root multiplicity\n");
            fflush(stdout);
            flint_abort();
        }
    }

    TEMPLATE(T, poly_roots)(roots, q, want_mult, ctx);
    if (roots->num > 0)
    {
        flint_printf("FAILED:\ncheck missing roots\n");
        fflush(stdout);
        flint_abort();
    }

    TEMPLATE(T, poly_clear)(q, ctx);
    TEMPLATE(T, poly_clear)(qt, ctx);
    TEMPLATE(T, poly_clear)(r, ctx);
}

TEST_TEMPLATE_FUNCTION_START(T, poly_factor_roots, state)
{
    slong i, j, k, l;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong n, m;

        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) f;
        TEMPLATE(T, poly_factor_t) r;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        TEMPLATE(T, poly_init)(f, ctx);
        TEMPLATE(T, poly_factor_init)(r, ctx);

        n = 1 + n_randint(state, 10);
        m = n_randint(state, 5);

        for (j = 0; j < 4; j++)
        {
            TEMPLATE(T, poly_randtest_not_zero)(f, state, n, ctx);

            for (k = 0; k < m; k++)
            {
                TEMPLATE(T, poly_t) ff;
                TEMPLATE(T, poly_init)(ff, ctx);
                TEMPLATE(T, poly_randtest_not_zero)(ff, state, 2, ctx);
                for (l = 1 + n_randint(state, 5); l > 0; l--)
                    TEMPLATE(T, poly_mul)(f, f, ff, ctx);
                TEMPLATE(T, poly_clear)(ff, ctx);
            }

            if (n_randint(state, 2))
            {
                test_poly(r, f, 1, ctx);
                test_poly(r, f, 0, ctx);
            }
            else
            {
                test_poly(r, f, 0, ctx);
                test_poly(r, f, 1, ctx);
            }
        }

        TEMPLATE(T, poly_factor_clear)(r, ctx);
        TEMPLATE(T, poly_clear)(f, ctx);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
