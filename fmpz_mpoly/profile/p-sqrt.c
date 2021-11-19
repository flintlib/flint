/*
    Copyright 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "profiler.h"

int iter_list[21] = {0, 3000000, 500000, 50000, 10000, 3000, 1000, 400, 100, 50, 30, 20, 10, 5, 2, 1, 1, 1, 1, 1, 1};

int main(void)
{
    fmpz_mpoly_t f, g, h, k;
    timeit_t timer;
    slong iters, j, n;
    fmpz_mpoly_ctx_t ctx;
    const char * vars[] = {"x", "y", "z", "t", "u"};

    fmpz_mpoly_ctx_init(ctx, 5, ORD_LEX);

    fmpz_mpoly_init(f, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(h, ctx);
    fmpz_mpoly_init(k, ctx);

    fmpz_mpoly_set_str_pretty(f, "(1+x+y+2*z^2+3*t^3+5*u^5)", vars, ctx);

    printf("Timing sqrt(p^2) where p = (1+x+y+2z^2+3t^3+5u^5)^n\n\n");
    printf("LEX ordering\n\n");

    for (n = 1; n <= 20; n++)
    {
        fmpz_mpoly_pow_ui(k, f, n, ctx);

        fmpz_mpoly_pow_ui(g, f, 2*n, ctx);

        iters = iter_list[n];

        timeit_start(timer);

        for (j = 0; j < iters; j++)
            fmpz_mpoly_sqrt(h, g, ctx);

        timeit_stop(timer);

        if (!fmpz_mpoly_equal(h, k, ctx))
        {
            printf("ERROR\n");
            flint_abort();
        }

        flint_printf("n = %wd: %.10lf s", n, (((double)timer->wall)/iters)/1000);
        fflush(stdout);

        timeit_start(timer);

        for (j = 0; j < iters; j++)
            fmpz_mpoly_sqrt_heap(h, g, ctx, 0);

        timeit_stop(timer);

        if (!fmpz_mpoly_equal(h, k, ctx))
        {
            printf("ERROR\n");
            flint_abort();
        }

        flint_printf("  unchecked %.10lf s\n", n, (((double)timer->wall)/iters)/1000);
    }

    fmpz_mpoly_clear(k, ctx);
    fmpz_mpoly_clear(h, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(f, ctx);
    fmpz_mpoly_ctx_clear(ctx);

    fmpz_mpoly_ctx_init(ctx, 5, ORD_DEGREVLEX);

    fmpz_mpoly_init(f, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(h, ctx);
    fmpz_mpoly_init(k, ctx);

    fmpz_mpoly_set_str_pretty(f, "(1+x+y+2*z^2+3*t^3+5*u^5)", vars, ctx);

    printf("Timing sqrt(p^2) where p = (1+x+y+2z^2+3t^3+5u^5)^n\n\n");
    printf("DEGREVLEX ordering:\n\n");

    for (n = 1; n <= 20; n++)
    {
        fmpz_mpoly_pow_ui(k, f, n, ctx);

        fmpz_mpoly_pow_ui(g, f, 2*n, ctx);

        iters = iter_list[n];

        timeit_start(timer);

        for (j = 0; j < iters; j++)
            fmpz_mpoly_sqrt(h, g, ctx);

        timeit_stop(timer);

        if (!fmpz_mpoly_equal(h, k, ctx))
        {
            printf("ERROR\n");
            flint_abort();
        }

        flint_printf("n = %wd: %.10lf s", n, (((double)timer->wall)/iters)/1000);
        fflush(stdout);

        timeit_start(timer);

        for (j = 0; j < iters; j++)
            fmpz_mpoly_sqrt_heap(h, g, ctx, 0);

        timeit_stop(timer);

        if (!fmpz_mpoly_equal(h, k, ctx))
        {
            printf("ERROR\n");
            flint_abort();
        }

        flint_printf("  unchecked %.10lf s\n", n, (((double)timer->wall)/iters)/1000);
    }

    fmpz_mpoly_clear(k, ctx);
    fmpz_mpoly_clear(h, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(f, ctx);
    fmpz_mpoly_ctx_clear(ctx);

    return 0;
}

