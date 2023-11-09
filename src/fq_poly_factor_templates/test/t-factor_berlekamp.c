/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
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

TEST_TEMPLATE_FUNCTION_START(T, poly_factor_berlekamp, state)
{
    int iter;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) poly1, poly, q, r, product;
        TEMPLATE(T, poly_factor_t) res;
        slong i, j, length, num;
        slong exp[5];

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (poly1, ctx);
        TEMPLATE(T, poly_init) (poly, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);

        TEMPLATE(T, poly_zero) (poly1, ctx);
        TEMPLATE(T, poly_one) (poly1, ctx);

        length = n_randint(state, 4) + 2;
        do
        {
            TEMPLATE(T, poly_randtest) (poly, state, length, ctx);
            if (poly->length)
                TEMPLATE(T, poly_make_monic) (poly, poly, ctx);
        }
        while ((poly->length < 2)
               || (!TEMPLATE(T, poly_is_irreducible) (poly, ctx)));

        exp[0] = n_randint(state, 5) + 1;
        for (i = 0; i < exp[0]; i++)
            TEMPLATE(T, poly_mul) (poly1, poly1, poly, ctx);

        num = n_randint(state, 5) + 1;

        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 5) + 2;
                TEMPLATE(T, poly_randtest) (poly, state, length, ctx);
                if (poly->length)
                {
                    TEMPLATE(T, poly_make_monic) (poly, poly, ctx);
                    TEMPLATE(T, poly_divrem) (q, r, poly1, poly, ctx);
                }
            }
            while ((poly->length < 2)
                   || (!TEMPLATE(T, poly_is_irreducible) (poly, ctx))
                   || (r->length == 0));

            exp[i] = n_randint(state, 5) + 1;
            for (j = 0; j < exp[i]; j++)
                TEMPLATE(T, poly_mul) (poly1, poly1, poly, ctx);
        }

        TEMPLATE(T, poly_factor_init) (res, ctx);
        TEMPLATE(T, poly_factor_berlekamp) (res, poly1, ctx);

        if (res->num != num)
        {
            flint_printf("Error: number of factors incorrect: %wd != %wd\n",
                         res->num, num);
            TEMPLATE(T, ctx_print) (ctx);
            flint_printf("\n");
            TEMPLATE(T, poly_print_pretty) (poly1, "x", ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_init) (product, ctx);
        TEMPLATE(T, poly_one) (product, ctx);
        for (i = 0; i < res->num; i++)
            for (j = 0; j < res->exp[i]; j++)
                TEMPLATE(T, poly_mul) (product, product, res->poly + i, ctx);

        TEMPLATE3(T, poly_scalar_mul, T) (product, product,
                                          poly1->coeffs + poly1->length - 1,
                                          ctx);

        if (!TEMPLATE(T, poly_equal) (poly1, product, ctx))
        {
            flint_printf
                ("Error: product of factors does not equal to the original polynomial\n");
            flint_printf("poly:\n");
            TEMPLATE(T, poly_print) (poly1, ctx);
            flint_printf("\n");
            flint_printf("product:\n");
            TEMPLATE(T, poly_print) (product, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (product, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);
        TEMPLATE(T, poly_clear) (poly1, ctx);
        TEMPLATE(T, poly_clear) (poly, ctx);
        TEMPLATE(T, poly_factor_clear) (res, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
