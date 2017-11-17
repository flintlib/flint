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
#include <gmp.h>
#include "flint.h"
#include "nmod_mpoly.h"

int
main(void)
{
    FLINT_TEST_INIT(state);

    flint_printf("print_parse....\n");
    fflush(stdout);

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        const char * vars[] = {"x","y","z","w","u","v"};

        nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, 6);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        nmod_mpoly_gen(f, 0, ctx);
        nmod_mpoly_scalar_mul_ui(f, f, 2, ctx);
        nmod_mpoly_add_ui(f, f, 1, ctx);
        printf("f: "); nmod_mpoly_print_pretty(f, vars, ctx); printf("\n");

        nmod_mpoly_gen(g, 1, ctx);
        nmod_mpoly_scalar_mul_ui(g, g, 3, ctx);
        nmod_mpoly_add_ui(g, g, 1, ctx);
        printf("g: "); nmod_mpoly_print_pretty(g, vars, ctx); printf("\n");

        nmod_mpoly_add(h, f, g, ctx);
        printf("h: "); nmod_mpoly_print_pretty(h, vars, ctx); printf("\n");

        nmod_mpoly_mul_johnson(h, f, g, ctx);
        printf("h: "); nmod_mpoly_print_pretty(h, vars, ctx); printf("\n");

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

