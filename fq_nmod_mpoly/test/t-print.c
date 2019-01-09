/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fq_nmod_mpoly.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("print_parse....");
    fflush(stdout);

    {
        char * fs;
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        const char * vars[] = {"x", "y", "z"};

        fq_nmod_mpoly_ctx_init(ctx, 3, 7, 4);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);

        fq_nmod_mpoly_gen(f, 0, ctx);
        fq_nmod_mpoly_gen(g, 1, ctx);
        fq_nmod_mpoly_gen(h, 2, ctx);

printf("f: "); fq_nmod_mpoly_print_pretty(f, NULL, ctx); printf("\n");
printf("g: "); fq_nmod_mpoly_print_pretty(g, NULL, ctx); printf("\n");
printf("h: "); fq_nmod_mpoly_print_pretty(h, NULL, ctx); printf("\n");

        fq_nmod_mpoly_set_str_pretty(f, "(# + 1)*x*y+#^2*z+x+1", vars, ctx);
printf("f: "); fq_nmod_mpoly_print_pretty(f, vars, ctx); printf("\n");

        fs = fq_nmod_mpoly_get_str_pretty(f, vars, ctx);

flint_printf("f: %s\n", fs);

        flint_free(fs);

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
