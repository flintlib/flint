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

    Copyright (C) 2012 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "fq_zech.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int j, i, result;
    flint_rand_t state;
    fq_zech_ctx_t ctx;

    printf("norm... ");
    fflush(stdout);

    flint_randinit(state);

    for (j = 0; j < 50; j++)
    {
        fq_zech_ctx_randtest(ctx, state);

        for (i = 0; i < 200; i++)
        {
            fmpz_t t1, t2;
            fq_nmod_t aa;
            fq_zech_t a;

            fmpz_init(t1);
            fmpz_init(t2);
            
            fq_nmod_init(aa, ctx->fq_nmod_ctx);
            fq_zech_init(a, ctx);

            fq_nmod_randtest(aa, state, ctx->fq_nmod_ctx);
            fq_zech_set_fq_nmod(a, aa, ctx);

            fq_nmod_norm(t1, aa, ctx->fq_nmod_ctx);
            fq_zech_norm(t2, a, ctx);

            result = fmpz_equal(t1, t2);
            if (!result)
            {
                printf("FAIL:\n\n");
                fq_zech_ctx_print(ctx);
                printf("\n");
                printf("aa = ");
                fq_nmod_print_pretty(aa, ctx->fq_nmod_ctx);
                printf("\n");
                printf("N(aa) = "); fmpz_print(t1);
                printf("\na = ");
                fq_zech_print_pretty(a, ctx);
                printf("\n");
                printf("N(a) = "); fmpz_print(t2);
                printf("\n");
                abort();
            }

            fq_zech_clear(a, ctx);
            fq_nmod_clear(aa, ctx->fq_nmod_ctx);
            fmpz_clear(t1);
            fmpz_clear(t2);
        }

        fq_zech_ctx_clear(ctx);
    }


    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
