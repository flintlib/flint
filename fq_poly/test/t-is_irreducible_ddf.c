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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "fq_poly.h"

int
main(void)
{
    int iter;
    flint_rand_t state;
    flint_randinit(state);

    flint_printf("is_irreducible_ddf....");
    fflush(stdout);

    for (iter = 0; iter < 50; iter++)
    {
        fq_ctx_t ctx;
        fq_poly_t poly1, poly2;
        slong length;
        int i, num;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(poly1, ctx);
        fq_poly_init(poly2, ctx);

        length = n_randint(state, 10) + 2;
        do
        {
            fq_poly_randtest(poly1, state, length, ctx);
            if (!fq_poly_is_zero(poly1, ctx))
                fq_poly_make_monic(poly1, poly1, ctx);
        }
        while ((!fq_poly_is_irreducible_ddf(poly1, ctx)) || (poly1->length < 2));

        num = n_randint(state, 5) + 1;

        for (i = 0; i < num; i++)
        {
            do
            {
                fq_poly_randtest(poly2, state, length, ctx);
                if (!fq_poly_is_zero(poly2, ctx))
                    fq_poly_make_monic(poly2, poly2, ctx);
            }
            while ((!fq_poly_is_irreducible_ddf(poly2, ctx)) || (poly2->length < 2));

            fq_poly_mul(poly1, poly1, poly2, ctx);
        }

        if (fq_poly_is_irreducible_ddf(poly1, ctx))
        {
            flint_printf("Error: reducible polynomial declared irreducible!\n");
            flint_printf("poly:\n");
            fq_poly_print(poly1, ctx);
            flint_printf("\n");
            abort();
        }

        fq_poly_clear(poly1, ctx);
        fq_poly_clear(poly2, ctx);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");
    return 0;
}
