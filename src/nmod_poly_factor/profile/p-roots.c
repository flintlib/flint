/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "profiler.h"


int main(int argc, char *argv[])
{
    slong i;

    {
        nmod_poly_t f, g;
        nmod_poly_factor_t r;
        flint_rand_t state;
        timeit_t timer;
        ulong p = n_nextprime(UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX), 1);

        flint_randinit(state);

        nmod_poly_init(f, p);
        nmod_poly_init(g, p);
        nmod_poly_factor_init(r);

        nmod_poly_one(f);

        for (i = 1; i <= 200; i++)
        {
            nmod_poly_fit_length(g, 2);
            g->coeffs[0] = n_randint(state, p);
            g->coeffs[1] = n_randint(state, p);
            g->coeffs[1] = FLINT_MAX(g->coeffs[1], UWORD(1));
            g->length = 2;

            nmod_poly_mul(f, f, g);

            timeit_start(timer);
            nmod_poly_roots(r, f, 0);
            timeit_stop(timer);
            flint_printf("degree %wd time: %wd %wd", i, timer->wall, timer->wall/i);

            timeit_start(timer);
            nmod_poly_factor(r, f);
            timeit_stop(timer);
            flint_printf("  factor time: %wd %wd", timer->wall, timer->wall/i);

            printf("\n");
        }

        nmod_poly_factor_clear(r);
        nmod_poly_clear(g);
        nmod_poly_clear(f);

        flint_randclear(state);
    }

    flint_cleanup_master();
    return 0;
}
