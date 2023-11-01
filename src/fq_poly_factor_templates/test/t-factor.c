/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
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
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_factor, state)
{
    int iter;

    /* Default algorithm */
    for (iter = 0; iter < flint_test_multiplier(); iter++)
    {
        int result = 1;
        TEMPLATE(T, poly_t) pol1, poly, quot, rem, product;
        TEMPLATE(T, poly_factor_t) res;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) randlead;
        TEMPLATE(T, t) lead;
        slong length, num, i, j;
        ulong exp[5];

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, init) (randlead, ctx);

        TEMPLATE(T, randtest_not_zero) (randlead, state, ctx);

        TEMPLATE(T, poly_init) (pol1, ctx);
        TEMPLATE(T, poly_init) (poly, ctx);
        TEMPLATE(T, poly_init) (quot, ctx);
        TEMPLATE(T, poly_init) (rem, ctx);

        TEMPLATE(T, poly_zero) (pol1, ctx);
        TEMPLATE3(T, poly_set, T) (pol1, randlead, ctx);

        length = n_randint(state, 4) + 2;
        TEMPLATE(T, poly_randtest_irreducible) (poly, state, length, ctx);

        exp[0] = n_randint(state, 3) + 1;
        for (i = 0; i < exp[0]; i++)
            TEMPLATE(T, poly_mul) (pol1, pol1, poly, ctx);

        num = n_randint(state, 3) + 1;
        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 3) + 2;
                TEMPLATE(T, poly_randtest_irreducible) (poly, state, length,
                                                        ctx);
                TEMPLATE(T, poly_divrem) (quot, rem, pol1, poly, ctx);
            }
            while ((poly->length < 2) || (rem->length == 0));

            exp[i] = n_randint(state, 3) + 1;

            for (j = 0; j < exp[i]; j++)
                TEMPLATE(T, poly_mul) (pol1, pol1, poly, ctx);
        }

        TEMPLATE(T, poly_factor_init) (res, ctx);

        TEMPLATE(T, init) (lead, ctx);

        switch (n_randint(state, 4))
        {
            case 0:
                TEMPLATE(T, poly_factor) (res, lead, pol1, ctx);
                break;
            case 1:
                TEMPLATE(T, poly_factor_with_berlekamp) (res, lead, pol1, ctx);
                break;
            case 2:
                if (fmpz_is_even(TEMPLATE(T, ctx_prime) (ctx)))
                    TEMPLATE(T, poly_factor) (res, lead, pol1, ctx);
                else
                    TEMPLATE(T, poly_factor_with_cantor_zassenhaus) (res, lead,
                                                                     pol1,
                                                                     ctx);
                break;
            case 3:
                TEMPLATE(T, poly_factor_with_kaltofen_shoup) (res, lead, pol1,
                                                              ctx);
                break;

        }
        fflush(stdout);

        result &= (res->num == num);
        if (!result)
        {
            flint_printf("Error: number of factors incorrect, %wd, %wd\n",
                         res->num, num);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_init) (product, ctx);
        TEMPLATE(T, poly_one) (product, ctx);
        for (i = 0; i < res->num; i++)
            for (j = 0; j < res->exp[i]; j++)
                TEMPLATE(T, poly_mul) (product, product, res->poly + i, ctx);
        TEMPLATE(T, TEMPLATE(poly_scalar_mul, T)) (product, product, lead,
                                                   ctx);
        result &= TEMPLATE(T, poly_equal) (pol1, product, ctx);
        if (!result)
        {
            flint_printf
                ("Error: product of factors does not equal original polynomial\n");
            TEMPLATE(T, poly_print_pretty) (pol1, "x", ctx);
            flint_printf("\n");
            TEMPLATE(T, poly_print_pretty) (product, "x", ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
        TEMPLATE(T, poly_clear) (product, ctx);

        TEMPLATE(T, poly_clear) (quot, ctx);
        TEMPLATE(T, poly_clear) (rem, ctx);
        TEMPLATE(T, poly_clear) (pol1, ctx);
        TEMPLATE(T, poly_clear) (poly, ctx);
        TEMPLATE(T, poly_factor_clear) (res, ctx);
        TEMPLATE(T, clear) (lead, ctx);
        TEMPLATE(T, clear) (randlead, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Test deflation trick */
    for (iter = 0; iter < flint_test_multiplier(); iter++)
    {
        TEMPLATE(T, poly_t) pol1, poly, quot, rem;
        TEMPLATE(T, poly_factor_t) res, res2;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) lead;
        slong length, num, i, j;
        slong exp[5];
        ulong inflation;
        int found;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (pol1, ctx);
        TEMPLATE(T, poly_init) (poly, ctx);
        TEMPLATE(T, poly_init) (quot, ctx);
        TEMPLATE(T, poly_init) (rem, ctx);

        TEMPLATE(T, poly_zero) (pol1, ctx);
        TEMPLATE(T, poly_one) (pol1, ctx);

        inflation = n_randint(state, 7) + 1;

        length = n_randint(state, 4) + 2;
        TEMPLATE(T, poly_randtest_irreducible) (poly, state, length, ctx);
        TEMPLATE(T, poly_inflate) (poly, poly, inflation, ctx);

        exp[0] = n_randint(state, 6) + 1;
        for (i = 0; i < exp[0]; i++)
            TEMPLATE(T, poly_mul) (pol1, pol1, poly, ctx);

        num = n_randint(state, 5) + 1;
        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 6) + 2;
                TEMPLATE(T, poly_randtest_irreducible) (poly, state, length,
                                                        ctx);
                TEMPLATE(T, poly_divrem) (quot, rem, pol1, poly, ctx);
            }
            while ((poly->length < 2) || (rem->length == 0));
            exp[i] = n_randint(state, 6) + 1;
            TEMPLATE(T, poly_inflate) (poly, poly, inflation, ctx);

            for (j = 0; j < exp[i]; j++)
                TEMPLATE(T, poly_mul) (pol1, pol1, poly, ctx);
        }

        TEMPLATE(T, poly_factor_init) (res, ctx);
        TEMPLATE(T, poly_factor_init) (res2, ctx);

        TEMPLATE(T, init) (lead, ctx);

        switch (n_randint(state, 4))
        {
            case 0:
                TEMPLATE(T, poly_factor) (res, lead, pol1, ctx);
                break;
            case 1:
                TEMPLATE(T, poly_factor_with_berlekamp) (res, lead, pol1, ctx);
                break;
            case 2:
                TEMPLATE(T, poly_factor_with_cantor_zassenhaus) (res, lead,
                                                                 pol1, ctx);
                break;
            case 3:
                TEMPLATE(T, poly_factor_with_kaltofen_shoup) (res, lead, pol1,
                                                              ctx);
                break;
        }

        TEMPLATE(T, poly_factor_cantor_zassenhaus) (res2, pol1, ctx);

        if (res->num != res2->num)
        {
            flint_printf("FAIL: different number of factors found\n");
            fflush(stdout);
            flint_abort();
        }

        for (i = 0; i < res->num; i++)
        {
            found = 0;
            for (j = 0; j < res2->num; j++)
            {
                if (TEMPLATE(T, poly_equal)
                    (res->poly + i, res2->poly + j, ctx)
                    && res->exp[i] == res2->exp[j])
                {
                    found = 1;
                    break;
                }
            }

            if (!found)
            {
                flint_printf("FAIL: factor not found\n");
                fflush(stdout);
                flint_abort();
            }
        }

        TEMPLATE(T, poly_clear) (quot, ctx);
        TEMPLATE(T, poly_clear) (rem, ctx);
        TEMPLATE(T, poly_clear) (pol1, ctx);
        TEMPLATE(T, poly_clear) (poly, ctx);
        TEMPLATE(T, poly_factor_clear) (res, ctx);
        TEMPLATE(T, poly_factor_clear) (res2, ctx);

        TEMPLATE(T, clear) (lead, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
