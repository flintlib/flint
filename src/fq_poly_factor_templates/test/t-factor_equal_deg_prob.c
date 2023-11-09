/*
    Copyright (C) 2012 Lina Kulakova
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

TEST_TEMPLATE_FUNCTION_START(T, poly_factor_equal_deg_prob, state)
{
    int iter;

    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) poly1, poly2, q, r;
        slong length;
        int i, num;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, poly_init) (poly1, ctx);
        TEMPLATE(T, poly_init) (poly2, ctx);

        length = n_randint(state, 10) + 2;
        do
        {
            TEMPLATE(T, poly_randtest) (poly1, state, length, ctx);
            if (poly1->length)
                TEMPLATE(T, poly_make_monic) (poly1, poly1, ctx);
        }
        while ((poly1->length != length)
               || (!TEMPLATE(T, poly_is_irreducible) (poly1, ctx)));

        num = n_randint(state, 5) + 1;

        for (i = 0; i < num; i++)
        {
            do
            {
                TEMPLATE(T, poly_randtest) (poly2, state, length, ctx);
                if (poly2->length)
                    TEMPLATE(T, poly_make_monic) (poly2, poly2, ctx);
            }
            while ((poly2->length != length)
                   || (!TEMPLATE(T, poly_is_irreducible) (poly2, ctx)));

            TEMPLATE(T, poly_mul) (poly1, poly1, poly2, ctx);
        }

        while (!TEMPLATE(T, poly_factor_equal_deg_prob)
               (poly2, state, poly1, length - 1, ctx))
        {
        };

        TEMPLATE(T, poly_divrem) (q, r, poly1, poly2, ctx);
        if (!TEMPLATE(T, poly_is_zero) (r, ctx))
        {
            flint_printf("FAIL:\n");
            flint_printf
                ("Error: factor does not divide original polynomial\n");
            flint_printf("factor:\n");
            TEMPLATE(T, poly_print) (poly2, ctx);
            flint_printf("\n\n");
            flint_printf("polynomial:\n");
            TEMPLATE(T, poly_print) (poly1, ctx);
            flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);
        TEMPLATE(T, poly_clear) (poly1, ctx);
        TEMPLATE(T, poly_clear) (poly2, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
