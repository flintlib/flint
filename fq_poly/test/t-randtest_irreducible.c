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


#include "fq_poly.h"

int
main(void)
{
    int iter;
    flint_rand_t state;
    flint_randinit(state);

    printf("randtest_irreducible....");
    fflush(stdout);

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        fq_ctx_t ctx;
        fq_poly_t poly;
        slong length;

        fq_ctx_randtest(ctx, state);

        fq_poly_init(poly);

        length = n_randint(state, 20) + 2;
        fq_poly_randtest_irreducible(poly, state, length, ctx);

        if (!fq_poly_is_irreducible(poly, ctx))
        {
            printf("Error: reducible polynomial created!\n");
            printf("poly:\n");
            fq_poly_print_pretty(poly, "x", ctx);
            printf("\n");
            abort();
        }

        fq_poly_clear(poly);

        fq_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
