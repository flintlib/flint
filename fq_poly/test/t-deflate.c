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

    flint_printf("deflate....");
    fflush(stdout);

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fq_poly_t poly1, poly2, poly3;
        fq_ctx_t ctx;
        ulong infl1, infl, deflation;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(poly1, ctx);
        fq_poly_init(poly2, ctx);
        fq_poly_init(poly3, ctx);

        fq_poly_randtest(poly1, state, n_randint(state, 15), ctx);

        if (fq_poly_length(poly1, ctx) <= 1)
        {
            if (fq_poly_deflation(poly1, ctx) != fq_poly_length(poly1, ctx))
            {
                flint_printf("FAIL: wrong deflation for constant polynomial\n");
                abort();
            }

            fq_poly_deflate(poly2, poly1, n_randint(state, 5) + 1, ctx);
            if (!fq_poly_equal(poly2, poly1, ctx))
            {
                flint_printf("FAIL: constant polynomial changed on deflation\n");
                abort();
            }
        }
        else
        {

            infl = n_randint(state, 13) + 1;
            infl1 = fq_poly_deflation(poly1, ctx);

            fq_poly_inflate(poly2, poly1, infl, ctx);

            deflation = fq_poly_deflation(poly2, ctx);

            if (deflation != infl * infl1)
            {
                flint_printf("FAIL: deflation = %wu, inflation: %wu, %wu\n",
                    deflation, infl, infl1);
                flint_printf("poly1:\n"); fq_poly_print(poly1, ctx); flint_printf("\n\n");
                flint_printf("poly2:\n"); fq_poly_print(poly2, ctx); flint_printf("\n\n");
                abort();
            }

            fq_poly_deflate(poly3, poly2, infl, ctx);
            if (!fq_poly_equal(poly3, poly1, ctx))
            {
                flint_printf("FAIL: deflation = %wu, inflation: %wu, %wu\n",
                    deflation, infl, infl1);
                flint_printf("Deflated polynomial not equal to input:\n");
                flint_printf("poly1:\n"); fq_poly_print(poly1, ctx); flint_printf("\n\n");
                flint_printf("poly2:\n"); fq_poly_print(poly2, ctx); flint_printf("\n\n");
                flint_printf("poly3:\n"); fq_poly_print(poly3, ctx); flint_printf("\n\n");
                abort();
            }

            fq_poly_deflate(poly2, poly2, infl, ctx);
            if (!fq_poly_equal(poly3, poly2, ctx))
            {
                flint_printf("FAIL: aliasing\n");
                abort();
            }
        }

        fq_poly_clear(poly1, ctx);
        fq_poly_clear(poly2, ctx);
        fq_poly_clear(poly3, ctx);
        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_printf("PASS\n");
    return 0;
}
