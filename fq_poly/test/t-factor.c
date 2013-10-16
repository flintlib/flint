/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

int
main(void)
{
    int iter;
    flint_rand_t state;
    flint_randinit(state);

    flint_printf("factor....");
    fflush(stdout);

    /* Default algorithm */
    for (iter = 0; iter < 2 * flint_test_multiplier(); iter++)
    {
        int result = 1;
        fq_poly_t pol1, poly, quot, rem, product;
        fq_poly_factor_t res;
        fq_ctx_t ctx;
        fq_t lead;
        slong length, num, i, j;
        ulong exp[5], prod1;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(pol1);
        fq_poly_init(poly);
        fq_poly_init(quot);
        fq_poly_init(rem);

        fq_poly_zero(pol1, ctx);
        fq_poly_one(pol1, ctx);

        length = n_randint(state, 7) + 2;
        fq_poly_randtest_irreducible(poly, state, length, ctx);

        exp[0] = n_randint(state, 30) + 1;
        prod1 = exp[0];
        for (i = 0; i < exp[0]; i++)
            fq_poly_mul(pol1, pol1, poly, ctx);

        num = n_randint(state, 5) + 1;
        for (i = 1; i < num; i++)
        {
            do 
            {
                length = n_randint(state, 7) + 2;
                fq_poly_randtest_irreducible(poly, state, length, ctx);
                fq_poly_divrem(quot, rem, pol1, poly, ctx);
            }
            while ((poly->length < 2) || (rem->length == 0));
            
            exp[i] = n_randint(state, 30) + 1;
            prod1 *= exp[i];

            for (j = 0; j < exp[i]; j++)
                fq_poly_mul(pol1, pol1, poly, ctx);
        }

        fq_poly_factor_init(res, ctx);

        fq_init(lead);

        switch (n_randint(state, 4))
        {
            case 0:
                fq_poly_factor(res, lead, pol1, ctx);
                break;
            case 1:
                fq_poly_factor_with_berlekamp(res, lead, pol1, ctx);
                break;
            case 2:
                if (fmpz_is_even(fq_ctx_prime(ctx)))
                    fq_poly_factor(res, lead, pol1, ctx);
                else
                    fq_poly_factor_with_cantor_zassenhaus(res, lead, pol1, ctx);
                break;
            case 3:
                fq_poly_factor_with_kaltofen_shoup(res, lead, pol1, ctx);
                break;

        }
        fflush(stdout);

        result &= (res->num == num);
        if (!result)
        {
            flint_printf("Error: number of factors incorrect, %wd, %wd\n",
                         res->num, num);
            abort();
        }

        fq_poly_init(product);
        fq_poly_one(product, ctx);
        for (i = 0; i < res->num; i++)
            for (j = 0; j < res->exp[i]; j++)
                fq_poly_mul(product, product, res->poly + i, ctx);
        fq_poly_scalar_mul_fq(product, product, lead, ctx);
        result &= fq_poly_equal(pol1, product);
        if (!result)
        {
            flint_printf("Error: product of factors does not equal original polynomial\n");
            fq_poly_print_pretty(pol1, "x", ctx); flint_printf("\n");
            fq_poly_print_pretty(product, "x", ctx); flint_printf("\n");
            abort();
        }
        fq_poly_clear(product);

        fq_poly_clear(quot);
        fq_poly_clear(rem);
        fq_poly_clear(pol1);
        fq_poly_clear(poly);
        fq_poly_factor_clear(res);
        fq_clear(lead);
        fq_ctx_clear(ctx);
    }

    /* Test deflation trick */
    for (iter = 0; iter < 2 * flint_test_multiplier(); iter++)
    {
        fq_poly_t pol1, poly, quot, rem;
        fq_poly_factor_t res, res2;
        fq_ctx_t ctx;
        fq_t lead;
        slong length, num, i, j;
        slong exp[5], prod1;
        ulong inflation;
        int found;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(pol1);
        fq_poly_init(poly);
        fq_poly_init(quot);
        fq_poly_init(rem);

        fq_poly_zero(pol1, ctx);
        fq_poly_one(pol1, ctx);

        inflation = n_randint(state, 7) + 1;

        length = n_randint(state, 7) + 2;
        fq_poly_randtest_irreducible(poly, state, length, ctx);
        fq_poly_inflate(poly, poly, inflation, ctx);

        exp[0] = n_randint(state, 6) + 1;
        prod1 = exp[0];
        for (i = 0; i < exp[0]; i++)
            fq_poly_mul(pol1, pol1, poly, ctx);

        num = n_randint(state, 5) + 1;
        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 6) + 2;
                fq_poly_randtest_irreducible(poly, state, length, ctx); 
                fq_poly_divrem(quot, rem, pol1, poly, ctx);
            }
            while ((poly->length < 2) || (rem->length == 0));
            exp[i] = n_randint(state, 6) + 1;
            prod1 *= exp[i];
            fq_poly_inflate(poly, poly, inflation, ctx);

            for (j = 0; j < exp[i]; j++)
                fq_poly_mul(pol1, pol1, poly, ctx);
        }

        fq_poly_factor_init(res, ctx);
        fq_poly_factor_init(res2, ctx);

        fq_init(lead);

        switch (n_randint(state, 4))
        {
            case 0:
                fq_poly_factor(res, lead, pol1, ctx);
                break;
            case 1:
                fq_poly_factor_with_berlekamp(res, lead, pol1, ctx);
                break;
            case 2:
                fq_poly_factor_with_cantor_zassenhaus(res, lead, pol1, ctx);
                break;
            case 3:
                fq_poly_factor_with_kaltofen_shoup(res, lead, pol1, ctx);
                break;
        }

        fq_poly_factor_cantor_zassenhaus(res2, pol1, ctx);

        if (res->num != res2->num)
        {
            flint_printf("FAIL: different number of factors found\n");
            abort();
        }

        for (i = 0; i < res->num; i++)
        {
            found = 0;
            for (j = 0; j < res2->num; j++)
            {
                if (fq_poly_equal(res->poly + i, res2->poly + j) &&
                        res->exp[i] == res2->exp[j])
                {
                    found = 1;
                    break;
                }
            }

            if (!found)
            {
                flint_printf("FAIL: factor not found\n");
                abort();
            }
        }

        fq_poly_clear(quot);
        fq_poly_clear(rem);
        fq_poly_clear(pol1);
        fq_poly_clear(poly);
        fq_poly_factor_clear(res);
        fq_poly_factor_clear(res2);

        fq_clear(lead);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_printf("PASS\n");
    return 0;
}
