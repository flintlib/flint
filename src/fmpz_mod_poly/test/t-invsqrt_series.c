/*
    Copyright (C) 2011, 2021 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_invsqrt_series, state)
{
    int i, result;

    /* Check 1/g^2 = h mod x^m */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t h, g, r;
        slong m;
        fmpz_mod_ctx_t ctx;

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, FLINT_BITS);

        if (!fmpz_equal_ui(fmpz_mod_ctx_modulus(ctx), 2))
        {
            fmpz_mod_poly_init(h, ctx);
            fmpz_mod_poly_init(g, ctx);
            fmpz_mod_poly_init(r, ctx);

            do fmpz_mod_poly_randtest(h, state, n_randint(state, 1000), ctx);
            while (h->length == 0);
            fmpz_mod_poly_set_coeff_ui(h, 0, 1, ctx);

            m = n_randint(state, h->length) + 1;

            fmpz_mod_poly_invsqrt_series(g, h, m, ctx);

            fmpz_mod_poly_mullow(r, g, g, m, ctx);
            fmpz_mod_poly_inv_series(r, r, m, ctx);
            fmpz_mod_poly_truncate(h, m, ctx);

            result = (fmpz_cmp_ui(ctx->n, 2) == 0) || (fmpz_mod_poly_equal(r, h, ctx));
            if (!result)
            {
                flint_printf("FAIL:\n");
                fmpz_mod_poly_print(h, ctx), flint_printf("\n\n");
                fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
                fmpz_mod_poly_print(r, ctx), flint_printf("\n\n");
                flint_printf("n = ");
                fmpz_print(ctx->n);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_poly_clear(h, ctx);
            fmpz_mod_poly_clear(g, ctx);
            fmpz_mod_poly_clear(r, ctx);
        }

        fmpz_mod_ctx_clear(ctx);
    }

    /* Check aliasing of h and g */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t g, h;
        slong m;
        fmpz_mod_ctx_t ctx;

        fmpz_mod_ctx_init_rand_bits_prime(ctx, state, FLINT_BITS);

        if (!fmpz_equal_ui(fmpz_mod_ctx_modulus(ctx), 2))
        {
            fmpz_mod_poly_init(h, ctx);
            fmpz_mod_poly_init(g, ctx);

            do fmpz_mod_poly_randtest(h, state, n_randint(state, 500), ctx);
            while (h->length == 0);

            fmpz_mod_poly_set_coeff_ui(h, 0, 1, ctx);

            m = n_randint(state, h->length) + 1;

            fmpz_mod_poly_invsqrt_series(g, h, m, ctx);
            fmpz_mod_poly_invsqrt_series(h, h, m, ctx);

            result = (fmpz_cmp_ui(ctx->n, 2) == 0) || fmpz_mod_poly_equal(g, h, ctx);
            if (!result)
            {
                flint_printf("FAIL:\n");
                fmpz_mod_poly_print(h, ctx), flint_printf("\n\n");
                fmpz_mod_poly_print(g, ctx), flint_printf("\n\n");
                flint_printf("n = ");
                fmpz_print(ctx->n);
                flint_printf("\n");
                flint_printf("m = %wd\n", m);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_poly_clear(g, ctx);
            fmpz_mod_poly_clear(h, ctx);
        }

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
