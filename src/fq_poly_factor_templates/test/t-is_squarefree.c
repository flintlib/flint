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
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_factor_is_squarefree, state)
{
    int iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) poly, Q, R, t;
        fmpz_t x;
        slong i, num_factors, exp, max_exp;
        int v, result;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (poly, ctx);
        TEMPLATE(T, poly_init) (t, ctx);
        TEMPLATE(T, poly_init) (Q, ctx);
        TEMPLATE(T, poly_init) (R, ctx);

        fmpz_init(x);
        fmpz_randtest_mod(x, state, TEMPLATE(T, ctx_prime) (ctx));

        TEMPLATE(T, poly_set_coeff_fmpz) (poly, 0, x, ctx);
        num_factors = n_randint(state, 5);

        max_exp = 0;
        for (i = 0; i < num_factors; i++)
        {
            do
            {
                TEMPLATE(T, poly_randtest) (t, state, n_randint(state, 10),
                                            ctx);
            } while (!TEMPLATE(T, poly_is_irreducible) (t, ctx)
                     || (TEMPLATE(T, poly_length) (t, ctx) < 2));

            exp = n_randint(state, 4) + 1;
            if (n_randint(state, 2) == 0)
                exp = 1;

            TEMPLATE(T, poly_divrem) (Q, R, poly, t, ctx);
            if (!TEMPLATE(T, poly_is_zero) (R, ctx))
            {
                TEMPLATE(T, poly_pow) (t, t, exp, ctx);
                TEMPLATE(T, poly_mul) (poly, poly, t, ctx);
                max_exp = FLINT_MAX(exp, max_exp);
            }
        }

        v = TEMPLATE(T, poly_is_squarefree) (poly, ctx);

        if (v == 1)
            result = (max_exp <= 1 && !TEMPLATE(T, poly_is_zero) (poly, ctx));
        else
            result = (max_exp > 1 || TEMPLATE(T, poly_is_zero) (poly, ctx));

        if (!result)
        {
            flint_printf("FAIL: ");
            TEMPLATE(T, ctx_print) (ctx);
            flint_printf(" %wd, %d\n", max_exp, v);
            TEMPLATE(T, poly_print) (poly, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (poly, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, poly_clear) (Q, ctx);
        TEMPLATE(T, poly_clear) (R, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
