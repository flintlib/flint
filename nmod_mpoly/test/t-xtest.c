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
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"
#include "profiler.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("print_parse....");
    fflush(stdout);


printf("\n");
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, fp, gp, h;
        timeit_t time;
        const char * vars[] = {"x","y","z","t","u","v"};

        fmpz_mpoly_ctx_init(ctx, 5, ORD_LEX);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);

        fmpz_mpoly_set_str_pretty(f, "1 + x + y^2 + z + t + u", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "1 + u + t^2 + z^3 + y + x", vars, ctx);
        fmpz_mpoly_pow_fps(fp, f, 3, ctx);
        fmpz_mpoly_pow_fps(gp, g, 3, ctx);

for (i = 0; i < 4; i++)
{
timeit_start(time);
        fmpz_mpoly_mul_johnson(h, fp, gp, ctx);
timeit_stop(time);
flint_printf("time: %wd\n",time->wall);
}

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(fp, ctx);
        fmpz_mpoly_clear(gp, ctx);

        fmpz_mpoly_ctx_clear(ctx);
    }

    {
        slong exp_bound, length;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, fp, gp, h;
        timeit_t time;
        const char * vars[] = {"x","y","z","t","u","v"};

        nmod_mpoly_ctx_init(ctx, 5, ORD_LEX, 17);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(fp, ctx);
        nmod_mpoly_init(gp, ctx);

        nmod_mpoly_set_str_pretty(f, "2*x", vars, ctx);
        printf("f: "); nmod_mpoly_print_pretty(f, vars, ctx); printf("\n");


        nmod_mpoly_set_str_pretty(f, "1 + x + 2*y^2 + 3*z + 4*t + 5*u", vars, ctx);
        nmod_mpoly_set_str_pretty(g, "1 + u + 2*t^2 + 3*z^3 + 4*y + 5*x", vars, ctx);

        printf("f: "); nmod_mpoly_print_pretty(f, vars, ctx); printf("\n");
        printf("g: "); nmod_mpoly_print_pretty(g, vars, ctx); printf("\n");


        nmod_mpoly_pow_rmul(fp, f, 3, ctx);
        nmod_mpoly_pow_rmul(gp, g, 3, ctx);

for (i = 0; i < 4; i++)
{
timeit_start(time);
        nmod_mpoly_mul_johnson(h, fp, gp, ctx);
timeit_stop(time);
flint_printf("time: %wd\n",time->wall);
}

        length = 6;
        exp_bound = 5;
        nmod_mpoly_randtest(f, state, length, exp_bound, ctx);
        printf("f: "); nmod_mpoly_print_pretty(f, vars, ctx); printf("\n");

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(fp, ctx);
        nmod_mpoly_clear(gp, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

