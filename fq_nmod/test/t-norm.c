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

#include "fq_nmod.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("norm... ");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with product of Galois conjugates */
    for (i = 0; i < 2000; i++)
    {
        fq_nmod_ctx_t ctx;

        fq_nmod_t a, b, c;
        mp_limb_t x, y;
        long j;

        fq_nmod_ctx_randtest(ctx, state);

        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx);
        fq_nmod_init(c, ctx);

        fq_nmod_randtest(a, state, ctx);
        fq_nmod_reduce(a, ctx);

        x = fq_nmod_norm(a, ctx);

        fq_nmod_one(b, ctx);
        for (j = 0; j < fq_nmod_ctx_degree(ctx); j++)
        {
            fq_nmod_frobenius(c, a, j, ctx);
            fq_nmod_mul(b, b, c, ctx);
        }
        y = nmod_poly_get_coeff_ui(b, 0);

        result = (x == y);
        if (!result)
        {
            printf("FAIL:\n\n");
            printf("a = "), fq_nmod_print_pretty(a, ctx), printf("\n");
            printf("b = "), fq_nmod_print_pretty(b, ctx), printf("\n");
            printf("x = %lu\n", x);
            printf("y = %lu\n", y);
            for (j = 0; j < fq_nmod_ctx_degree(ctx); j++)
            {
                fq_nmod_frobenius(c, a, j, ctx);
                printf("sigma^%ld = ", j), fq_nmod_print_pretty(c, ctx), printf("\n");
            }
            abort();
        }

        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);
        fq_nmod_clear(c, ctx);
        fq_nmod_ctx_clear(ctx);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
