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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include <stdlib.h>
#include "ulong_extras.h"
#include "fq_poly.h"

int
main(void)
{
    int iter;
    flint_rand_t state;
    flint_randinit(state);

    flint_printf("factor_distinct_deg....");
    fflush(stdout);

    for (iter = 0; iter < 200; iter++)
    {
        fq_ctx_t ctx;
        fq_poly_t poly1, poly, q, r, product;
        fq_poly_factor_t res;
        slong i, length, num;
        slong *degs;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(poly1);
        fq_poly_init(poly);
        fq_poly_init(q);
        fq_poly_init(r);

        fq_poly_zero(poly1, ctx);
        fq_poly_one(poly1, ctx);

        length = n_randint(state, 7) + 2;
        do
        {
            fq_poly_randtest(poly, state, length, ctx);
            if (poly->length)
                fq_poly_make_monic(poly, poly, ctx);
        }
        while ((poly->length < 2) || (!fq_poly_is_irreducible(poly, ctx)));

        fq_poly_mul(poly1, poly1, poly, ctx);

        num = n_randint(state, 5) + 1;

        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 7) + 2;
                fq_poly_randtest(poly, state, length, ctx);
                if (poly->length)
                {
                    fq_poly_make_monic(poly, poly, ctx);
                    fq_poly_divrem(q, r, poly1, poly, ctx);
                }
            }
            while ((poly->length < 2) || (!fq_poly_is_irreducible(poly, ctx))
                   || (r->length == 0));

            fq_poly_mul(poly1, poly1, poly, ctx);
        }

        if (!(degs = flint_malloc((poly1->length - 1) * sizeof(slong))))
        {
            flint_printf("Fatal error: not enough memory.");
            abort();
        }
        fq_poly_factor_init(res, ctx);
        fq_poly_factor_distinct_deg(res, poly1, &degs, ctx);

        fq_poly_init(product);
        fq_poly_one(product, ctx);
        for (i = 0; i < res->num; i++)
            fq_poly_mul(product, product, res->poly + i, ctx);

        fq_poly_scalar_mul_fq(product, product,
                              poly1->coeffs + (poly1->length - 1), ctx);

        if (!fq_poly_equal(poly1, product))
        {
            flint_printf
                ("Error: product of factors does not equal to the original polynomial\n");
            flint_printf("poly:\n");
            fq_poly_print(poly1, ctx);
            flint_printf("\n");
            flint_printf("product:\n");
            fq_poly_print(product, ctx);
            flint_printf("\n");
            abort();
        }

        flint_free(degs);
        fq_poly_clear(product);
        fq_poly_clear(q);
        fq_poly_clear(r);
        fq_poly_clear(poly1);
        fq_poly_clear(poly);
        fq_poly_factor_clear(res);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return 0;
}
