/*
    Copyright (C) 2011 Fredrik Johansson
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

TEST_TEMPLATE_FUNCTION_START(T, poly_factor_is_irreducible_ddf, state)
{
    int iter;

    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) poly1, poly2;
        slong length;
        int i, num;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (poly1, ctx);
        TEMPLATE(T, poly_init) (poly2, ctx);

        length = n_randint(state, 7) + 2;
        do
        {
            TEMPLATE(T, poly_randtest) (poly1, state, length, ctx);
            if (!TEMPLATE(T, poly_is_zero) (poly1, ctx))
                TEMPLATE(T, poly_make_monic) (poly1, poly1, ctx);
        }
        while ((!TEMPLATE(T, poly_is_irreducible_ddf) (poly1, ctx))
               || (poly1->length < 2));

        num = n_randint(state, 5) + 1;

        for (i = 0; i < num; i++)
        {
            do
            {
                TEMPLATE(T, poly_randtest) (poly2, state, length, ctx);
                if (!TEMPLATE(T, poly_is_zero) (poly2, ctx))
                    TEMPLATE(T, poly_make_monic) (poly2, poly2, ctx);
            }
            while ((!TEMPLATE(T, poly_is_irreducible_ddf) (poly2, ctx))
                   || (poly2->length < 2));

            TEMPLATE(T, poly_mul) (poly1, poly1, poly2, ctx);
        }

        if (TEMPLATE(T, poly_is_irreducible_ddf) (poly1, ctx))
        {
            flint_printf
                ("Error: reducible polynomial declared irreducible!\n");
            flint_printf("poly:\n");
            TEMPLATE(T, poly_print) (poly1, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (poly1, ctx);
        TEMPLATE(T, poly_clear) (poly2, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
