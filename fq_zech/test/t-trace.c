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

#include "fq_zech.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int j, i, result;
    flint_rand_t state;
    fq_zech_ctx_t ctx;

    printf("trace... ");
    fflush(stdout);

    flint_randinit(state);

    for (j = 0; j < 50; j++)
    {
        fq_zech_ctx_randtest(ctx, state);

        for (i = 0; i < 200; i++)
        {
            mp_limb_t t1, t2;
            fq_nmod_t aa;
            fq_zech_t a;

            fq_nmod_init(aa, ctx->fq_nmod_ctx);
            fq_zech_init(a, ctx);

            fq_nmod_randtest(aa, state, ctx->fq_nmod_ctx);
            fq_zech_set_fq_nmod(a, aa, ctx);

            t1 = fq_nmod_trace(aa, ctx->fq_nmod_ctx);
            t2 = fq_zech_trace(a, ctx);

            result = (t1 == t2);
            if (!result)
            {
                printf("FAIL:\n\n");
                fq_zech_ctx_print(ctx);
                printf("\n");
                printf("aa = ");
                fq_nmod_print_pretty(aa, ctx->fq_nmod_ctx);
                printf("\n");
                printf("Tr(aa) = %lu\n", t1);
                printf("a = ");
                fq_zech_print_pretty(a, ctx);
                printf("\n");
                printf("Tr(a) = %lu\n", t2);
                abort();
            }

            fq_zech_clear(a, ctx);
            fq_nmod_clear(aa, ctx->fq_nmod_ctx);
        }

        fq_zech_ctx_clear(ctx);
    }


    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
