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

    flint_printf("factor_equal_deg_prob....");
    fflush(stdout);

    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        fq_ctx_t ctx;
        fq_poly_t poly1, poly2, q, r;
        slong length;
        int i, num;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(q, ctx);
        fq_poly_init(r, ctx);
        fq_poly_init(poly1, ctx);
        fq_poly_init(poly2, ctx);

        length = n_randint(state, 10) + 2;
        do
        {
            fq_poly_randtest(poly1, state, length, ctx);
            if (poly1->length)
                fq_poly_make_monic(poly1, poly1, ctx);
        }
        while ((poly1->length != length)
               || (!fq_poly_is_irreducible(poly1, ctx)));

        num = n_randint(state, 5) + 1;

        for (i = 0; i < num; i++)
        {
            do
            {
                fq_poly_randtest(poly2, state, length, ctx);
                if (poly2->length)
                    fq_poly_make_monic(poly2, poly2, ctx);
            }
            while ((poly2->length != length)
                   || (!fq_poly_is_irreducible(poly2, ctx)));

            fq_poly_mul(poly1, poly1, poly2, ctx);
        }

        while (!fq_poly_factor_equal_deg_prob
               (poly2, state, poly1, length - 1, ctx))
        {
        };

        fq_poly_divrem(q, r, poly1, poly2, ctx);
        if (!fq_poly_is_zero(r, ctx))
        {
            flint_printf("FAIL:\n");
            flint_printf
                ("Error: factor does not divide original polynomial\n");
            flint_printf("factor:\n");
            fq_poly_print(poly2, ctx);
            flint_printf("\n\n");
            flint_printf("polynomial:\n");
            fq_poly_print(poly1, ctx);
            flint_printf("\n\n");
            abort();
        }

        fq_poly_clear(q, ctx);
        fq_poly_clear(r, ctx);
        fq_poly_clear(poly1, ctx);
        fq_poly_clear(poly2, ctx);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return 0;
}
