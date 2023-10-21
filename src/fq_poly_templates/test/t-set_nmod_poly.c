/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"
#include "nmod_poly.h"
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_set_nmod_poly, state)
{
    int i, result;

    /* Check litfed polynomials by evaluating at random points */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) r, s;
        TEMPLATE(T, poly_t) a;
        nmod_poly_t b;
        ulong p;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest)(ctx, state);
        TEMPLATE(T, init)(r, ctx); TEMPLATE(T, init)(s, ctx);
        TEMPLATE(T, poly_init)(a, ctx);
        nmod_poly_init(b, fmpz_get_ui(TEMPLATE(T, ctx_prime)(ctx)));

        nmod_poly_randtest(b, state, len);
        p = n_randint(state, 10);

        TEMPLATE(T, poly_set_nmod_poly)(a, b, ctx);
        TEMPLATE(T, set_ui)(r, p, ctx);
        TEMPLATE3(T, poly_evaluate, T)(r, a, r, ctx);
        p = nmod_poly_evaluate_nmod(b, p);
        TEMPLATE(T, set_ui)(s, p, ctx);

        result = TEMPLATE(T, equal)(r, s, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), TEMPLATE(T, ctx_print(ctx)),
                flint_printf("\n");
            flint_printf("b = "), nmod_poly_print_pretty(b, "X"),
                flint_printf("\n");
            flint_printf("p = %u\n", p);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(r, ctx); TEMPLATE(T, clear)(s, ctx);
        nmod_poly_clear(b);
        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
