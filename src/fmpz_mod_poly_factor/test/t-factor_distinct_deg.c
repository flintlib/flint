/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
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

/* Defined in t-factor_distinct_deg.c and t-factor_distinct_deg_threaded.c */
#define MAX_DEG 7

TEST_FUNCTION_START(fmpz_mod_poly_factor_distinct_deg, state)
{
    int iter;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        fmpz_mod_poly_t poly1, poly, q, r, product;
        fmpz_mod_poly_factor_t res;
        fmpz_t modulus;
        slong i, length, num;
        slong *degs;
        slong num_of_deg[MAX_DEG + 1];

        for (i = 0; i < MAX_DEG + 1; i++)
            num_of_deg[i] = 0;

        fmpz_init(modulus);
        fmpz_set_ui(modulus, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, modulus);

        fmpz_mod_poly_init(poly1, ctx);
        fmpz_mod_poly_init(poly, ctx);
        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);

        fmpz_mod_poly_zero(poly1, ctx);
        fmpz_mod_poly_set_coeff_ui(poly1, 0, 1, ctx);

        length = n_randint(state, MAX_DEG) + 2;
        do
        {
            fmpz_mod_poly_randtest(poly, state, length, ctx);
            if (poly->length)
                fmpz_mod_poly_make_monic(poly, poly, ctx);
        }
        while ((poly->length < 2) || (!fmpz_mod_poly_is_irreducible(poly, ctx)));

        fmpz_mod_poly_mul(poly1, poly1, poly, ctx);

        num_of_deg[fmpz_mod_poly_degree(poly, ctx)]++;

        num = n_randint(state, 6) + 1;

        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, MAX_DEG) + 2;
                fmpz_mod_poly_randtest(poly, state, length, ctx);
                if (poly->length)
                {
                    fmpz_mod_poly_make_monic(poly, poly, ctx);
                    fmpz_mod_poly_divrem(q, r, poly1, poly, ctx);
                }
            }
            while ((poly->length < 2) || (!fmpz_mod_poly_is_irreducible(poly, ctx))
                   || (r->length == 0));

            fmpz_mod_poly_mul(poly1, poly1, poly, ctx);
            num_of_deg[fmpz_mod_poly_degree(poly, ctx)]++;
        }

        if (!(degs = flint_malloc((poly1->length - 1) * sizeof(slong))))
        {
            flint_printf("Fatal error: not enough memory.");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mod_poly_factor_init(res, ctx);
        fmpz_mod_poly_factor_distinct_deg(res, poly1, &degs, ctx);

        fmpz_mod_poly_init(product, ctx);
        fmpz_mod_poly_set_coeff_ui(product, 0, 1, ctx);
        for (i = 0; i < res->num; i++)
        {
            fmpz_mod_poly_mul(product, product, res->poly + i, ctx);

            if (fmpz_mod_poly_degree(res->poly + i, ctx) != degs[i]*num_of_deg[degs[i]])
            {
               flint_printf("Error: product of factors of degree %w incorrect\n", degs[i]);
               flint_printf("Degree %w != %w * %w\n", fmpz_mod_poly_degree(res->poly + i, ctx), degs[i], num_of_deg[degs[i]]);
               fflush(stdout);
               flint_abort();
            }
        }

        fmpz_mod_poly_scalar_mul_fmpz(product, product,
                                      &(poly1->coeffs[poly1->length - 1]), ctx);

        if (!fmpz_mod_poly_equal(poly1, product, ctx))
        {
            flint_printf
                ("Error: product of factors does not equal to the original polynomial\n");
            flint_printf("poly:\n");
            fmpz_mod_poly_print(poly1, ctx);
            flint_printf("\n");
            flint_printf("product:\n");
            fmpz_mod_poly_print(product, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(degs);
        fmpz_clear(modulus);
        fmpz_mod_poly_clear(product, ctx);
        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_mod_poly_clear(poly1, ctx);
        fmpz_mod_poly_clear(poly, ctx);
        fmpz_mod_poly_factor_clear(res, ctx);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
}
#undef MAX_DEG
