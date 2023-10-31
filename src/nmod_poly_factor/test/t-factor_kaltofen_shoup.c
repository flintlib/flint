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
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

TEST_FUNCTION_START(nmod_poly_factor_kaltofen_shoup, state)
{
    int iter;

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        nmod_poly_t poly1, poly, q, r, product;
        nmod_poly_factor_t res;
        mp_limb_t modulus, lead;
        slong i, j, length, num;
        slong exp[5];

        modulus = n_randtest_prime(state, 0);

        nmod_poly_init(poly1, modulus);
        nmod_poly_init(poly, modulus);
        nmod_poly_init(q, modulus);
        nmod_poly_init(r, modulus);

        nmod_poly_zero(poly1);
        nmod_poly_set_coeff_ui(poly1, 0, 1);

        length = n_randint(state, 7) + 2;
        do
        {
            nmod_poly_randtest(poly, state, length);
            if (poly->length)
                nmod_poly_make_monic(poly, poly);
        }
        while ((poly->length < 2) || (!nmod_poly_is_irreducible(poly)));

        exp[0] = n_randint(state, 30) + 1;
        for (i = 0; i < exp[0]; i++)
            nmod_poly_mul(poly1, poly1, poly);

        num = n_randint(state, 5) + 1;

        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 7) + 2;
                nmod_poly_randtest(poly, state, length);
                if (poly->length)
                {
                    nmod_poly_make_monic(poly, poly);
                    nmod_poly_divrem(q, r, poly1, poly);
                }
            }
            while ((poly->length < 2) || (!nmod_poly_is_irreducible(poly)) ||
                (r->length == 0));

            exp[i] = n_randint(state, 30) + 1;
            for (j = 0; j < exp[i]; j++)
                nmod_poly_mul(poly1, poly1, poly);
        }

        nmod_poly_factor_init(res);
        nmod_poly_factor_kaltofen_shoup(res, poly1);

        if (res->num != num)
        {
            flint_printf("Error: number of factors incorrect: %wd != %wd\n", res->num, num);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_init_preinv(product, poly1->mod.n, poly1->mod.ninv);
        nmod_poly_set_coeff_ui(product, 0, 1);
        for (i = 0; i < res->num; i++)
            for (j = 0; j < res->exp[i]; j++)
                nmod_poly_mul(product, product, res->p + i);

        lead = poly1->coeffs[poly1->length - 1];
        nmod_poly_scalar_mul_nmod(product, product, lead);

        if (!nmod_poly_equal(poly1, product))
        {
            flint_printf("Error: product of factors does not equal to the original polynomial\n");
            flint_printf("poly:\n"); nmod_poly_print(poly1); flint_printf("\n");
            flint_printf("product:\n"); nmod_poly_print(product); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(product);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
        nmod_poly_clear(poly1);
        nmod_poly_clear(poly);
        nmod_poly_factor_clear(res);
    }

    TEST_FUNCTION_END(state);
}
