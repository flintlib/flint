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
#include "ulong_extras.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_factor_squarefree, state)
{
    int iter;

    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        int result = 1;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) pol1, poly, quot, rem;
        TEMPLATE(T, poly_factor_t) res;
        slong exp[5], prod1;
        slong length, i, j, num;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (pol1, ctx);
        TEMPLATE(T, poly_init) (poly, ctx);
        TEMPLATE(T, poly_init) (quot, ctx);
        TEMPLATE(T, poly_init) (rem, ctx);

        TEMPLATE(T, poly_one) (pol1, ctx);

        length = n_randint(state, 5) + 2;

        do
        {
            TEMPLATE(T, poly_randtest) (poly, state, length, ctx);
            TEMPLATE(T, poly_make_monic) (poly, poly, ctx);
        }
        while ((poly->length != length)
               || (!TEMPLATE(T, poly_is_irreducible) (poly, ctx)));
        exp[0] = n_randprime(state, 5, 0);

        prod1 = exp[0];
        for (i = 0; i < exp[0]; i++)
            TEMPLATE(T, poly_mul) (pol1, pol1, poly, ctx);

        num = n_randint(state, 5) + 1;
        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 7) + 2;
                TEMPLATE(T, poly_randtest) (poly, state, length, ctx);
                if (poly->length)
                {
                    TEMPLATE(T, poly_make_monic) (poly, poly, ctx);
                    TEMPLATE(T, poly_divrem) (quot, rem, pol1, poly, ctx);
                }
            }
            while ((!TEMPLATE(T, poly_is_irreducible) (poly, ctx)) ||
                   (poly->length != length) || (rem->length == 0));

            do
                exp[i] = n_randprime(state, 5, 0);
            while (prod1 % exp[i] == 0);

            prod1 *= exp[i];
            for (j = 0; j < exp[i]; j++)
                TEMPLATE(T, poly_mul) (pol1, pol1, poly, ctx);
        }

        TEMPLATE(T, poly_factor_init) (res, ctx);
        TEMPLATE(T, poly_factor_squarefree) (res, pol1, ctx);

        result &= (res->num == num);
        if (result)
        {
            ulong prod2 = 1;
            for (i = 0; i < num; i++)
                prod2 *= res->exp[i];
            result &= (prod1 == prod2);
        }

        if (!result)
        {
            flint_printf("Error: exp don't match. Ctx = ");
            TEMPLATE(T, ctx_print) (ctx);
            flint_printf("\n");
            for (i = 0; i < res->num; i++)
                flint_printf("%wd ", res->exp[i]);
            flint_printf("\n");
            for (i = 0; i < num; i++)
                flint_printf("%wd ", exp[i]);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (quot, ctx);
        TEMPLATE(T, poly_clear) (rem, ctx);
        TEMPLATE(T, poly_clear) (pol1, ctx);
        TEMPLATE(T, poly_clear) (poly, ctx);
        TEMPLATE(T, poly_factor_clear) (res, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
