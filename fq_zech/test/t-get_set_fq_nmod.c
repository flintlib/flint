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
    int i, j, result;
    fq_zech_ctx_t ctx;
    FLINT_TEST_INIT(state);
    
    flint_printf("get_fq_nmod/set_fq_nmod... ");
    fflush(stdout);

    for (j = 0; j < 10; j++)
    {
        fq_zech_ctx_randtest(ctx, state);

        /* Check aliasing: a = -a */
        for (i = 0; i < 200; i++)
        {
            fq_zech_t a, b;
            fq_nmod_t c;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fq_nmod_init(c, ctx->fq_nmod_ctx);

            
            fq_zech_randtest(a, state, ctx);
            fq_zech_get_fq_nmod(c, a, ctx);
            fq_zech_set_fq_nmod(b, c, ctx);

            result = (fq_zech_equal(a, b, ctx));
            if (!result)
            {
                flint_printf("FAIL:n\n");
                fq_zech_ctx_print(ctx);
                flint_printf("\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), fq_nmod_print_pretty(c, ctx->fq_nmod_ctx), flint_printf("\n");
                flint_printf("table = %wd\n", ctx->eval_table[a->value]);
                abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(b, ctx);
            fq_nmod_clear(c, ctx->fq_nmod_ctx);
        }

        fq_zech_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
