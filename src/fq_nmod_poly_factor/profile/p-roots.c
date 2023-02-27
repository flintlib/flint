/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fq_nmod_poly.h"
#include "profiler.h"


int main(int argc, char *argv[])
{
    slong i;

    {
        fmpz_t p;
        fq_nmod_ctx_t ctx;
        fq_nmod_t lc;
        fq_nmod_poly_t f, g;
        fq_nmod_poly_factor_t r;
        flint_rand_t randstate;
        timeit_t timer;

        flint_randinit(randstate);

        fmpz_init_set_ui(p, UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX));
        fmpz_nextprime(p, p, 1);

        fq_nmod_ctx_init(ctx, p, 6, "a");

        fq_nmod_poly_init(f, ctx);
        fq_nmod_poly_init(g, ctx);
        fq_nmod_poly_factor_init(r, ctx);

        fq_nmod_poly_one(f, ctx);

        fq_nmod_init(lc, ctx);

        for (i = 1; i <= 32; i++)
        {
            fq_nmod_poly_fit_length(g, 2, ctx);
            fq_nmod_rand(g->coeffs + 0, randstate, ctx);
            fq_nmod_rand(g->coeffs + 1, randstate, ctx);
            if (fq_nmod_is_zero(g->coeffs + 1, ctx))
                fq_nmod_one(g->coeffs + 1, ctx);
            g->length = 2;

            fq_nmod_poly_mul(f, f, g, ctx);

            timeit_start(timer);
            fq_nmod_poly_roots(r, f, 0, ctx);
            timeit_stop(timer);
            flint_printf("degree %wd time: %wd %wd", i, timer->wall, timer->wall/i);

            timeit_start(timer);
            fq_nmod_poly_factor(r, lc, f, ctx);
            timeit_stop(timer);
            flint_printf("  factor time: %wd %wd", timer->wall, timer->wall/i);

            printf("\n");
        }

        fq_nmod_clear(lc, ctx);
        fq_nmod_poly_factor_clear(r, ctx);
        fq_nmod_poly_clear(g, ctx);
        fq_nmod_poly_clear(f, ctx);
        fq_nmod_ctx_clear(ctx);
        fmpz_clear(p);

        flint_randclear(randstate);
    }

    flint_cleanup_master();
    return 0;
}
