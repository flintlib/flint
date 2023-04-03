/*
    Copyright 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "profiler.h"

int iter_list[21] = {0, 300000, 50000, 5000, 1000, 300, 100, 40, 10, 5, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1};

int main(void)
{
    nmod_mpoly_t f, g, h, k;
    timeit_t timer;
    slong iters, j, n;
    nmod_mpoly_ctx_t ctx;
    mp_limb_t p = n_nextprime(UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX), 1);
    const char * vars[] = {"x", "y", "z", "t", "u"};

    nmod_mpoly_ctx_init(ctx, 5, ORD_LEX, p);

    nmod_mpoly_init(f, ctx);
    nmod_mpoly_init(g, ctx);
    nmod_mpoly_init(h, ctx);
    nmod_mpoly_init(k, ctx);

    nmod_mpoly_set_str_pretty(f, "(1+x+y+2*z^2+3*t^3+5*u^5)", vars, ctx);

    printf("Timing sqrt(p^2) where p = (1+x+y+2z^2+3t^3+5u^5)^n\n\n");
    printf("LEX ordering\n\n");

    for (n = 1; n <= 20; n++)
    {
        nmod_mpoly_pow_ui(g, f, 2*n, ctx);

        iters = iter_list[n];

        timeit_start(timer);

        for (j = 0; j < iters; j++)
            nmod_mpoly_sqrt(h, g, ctx);

        timeit_stop(timer);

        flint_printf("n = %wd: %.10lf s\n", n, (((double)timer->wall)/iters)/1000);
    }

    nmod_mpoly_clear(k, ctx);
    nmod_mpoly_clear(h, ctx);
    nmod_mpoly_clear(g, ctx);
    nmod_mpoly_clear(f, ctx);

    nmod_mpoly_ctx_clear(ctx);
    nmod_mpoly_ctx_init(ctx, 5, ORD_DEGREVLEX, p);

    nmod_mpoly_init(f, ctx);
    nmod_mpoly_init(g, ctx);
    nmod_mpoly_init(h, ctx);
    nmod_mpoly_init(k, ctx);

    nmod_mpoly_set_str_pretty(f, "(1+x+y+2*z^2+3*t^3+5*u^5)", vars, ctx);

    printf("Timing sqrt(p^2) where p = (1+x+y+2z^2+3t^3+5u^5)^n\n\n");
    printf("DEGREVLEX ordering:\n\n");

    for (n = 1; n <= 20; n++)
    {
        nmod_mpoly_pow_ui(g, f, 2*n, ctx);

        iters = iter_list[n];

        timeit_start(timer);

        for (j = 0; j < iters; j++)
            nmod_mpoly_sqrt(h, g, ctx);

        timeit_stop(timer);

        flint_printf("n = %wd: %.10lf s\n", n, (((double)timer->wall)/iters)/1000);
    }

    nmod_mpoly_clear(k, ctx);
    nmod_mpoly_clear(h, ctx);
    nmod_mpoly_clear(g, ctx);
    nmod_mpoly_clear(f, ctx);
    nmod_mpoly_ctx_clear(ctx);

    return 0;
}

