/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"

TEST_FUNCTION_START(fmpz_mod_poly_factor_is_irreducible_ddf, state)
{
    int iter;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        fmpz_mod_poly_t poly1, poly2;
        fmpz_t modulus;
        slong length;
        int i, num;

        fmpz_init_set_ui(modulus, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, modulus);

        fmpz_mod_poly_init(poly1, ctx);
        fmpz_mod_poly_init(poly2, ctx);

        length = n_randint(state, 10) + 2;
        do
        {
            fmpz_mod_poly_randtest(poly1, state, length, ctx);
            if (!fmpz_mod_poly_is_zero(poly1, ctx))
                fmpz_mod_poly_make_monic(poly1, poly1, ctx);
        }
        while ((!fmpz_mod_poly_is_irreducible_ddf(poly1, ctx)) || (poly1->length < 2));

        num = n_randint(state, 5) + 1;

        for (i = 0; i < num; i++)
        {
            do
            {
                fmpz_mod_poly_randtest(poly2, state, length, ctx);
                if (!fmpz_mod_poly_is_zero(poly1, ctx))
                    fmpz_mod_poly_make_monic(poly2, poly2, ctx);
            }
            while ((!fmpz_mod_poly_is_irreducible_ddf(poly2, ctx)) || (poly2->length < 2));

            fmpz_mod_poly_mul(poly1, poly1, poly2, ctx);
        }

        if (fmpz_mod_poly_is_irreducible_ddf(poly1, ctx))
        {
            flint_printf("Error: reducible polynomial declared irreducible!\n");
            flint_printf("poly:\n");
            fmpz_mod_poly_print(poly1, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(modulus);
        fmpz_mod_poly_clear(poly1, ctx);
        fmpz_mod_poly_clear(poly2, ctx);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
