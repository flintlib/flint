/*
    Copyright 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "fmpz_mpoly_factor.h"


slong check_omega(slong om, const fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_factor_t g;
    fmpz_t omega;
    timeit_t timer;

    fmpz_init(omega);
    fmpz_mpoly_factor_init(g, ctx);

    timeit_start(timer);
    if (!fmpz_mpoly_factor(g, p, ctx))
    {
        flint_printf("oops! could not factor\n");
        flint_abort();
    }
    timeit_stop(timer);

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, om) != 0)
    {
        flint_printf("factorization has wrong number of factors\n");
        flint_abort();        
    }

    fmpz_mpoly_factor_clear(g, ctx);
    fmpz_clear(omega);

    return timer->wall;
}


int main(int argc, char *argv[])
{
    slong i, j, k;
    slong time, total_time;

    flint_printf("\n------ 4 variables ------\n");
    total_time = 0;
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c;
        const char * vars[] = {"x", "y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);

        for (i = 0; i <= 20; i++)
        {
            fmpz_mpoly_set_str_pretty(a, "x", vars, ctx);
            fmpz_mpoly_set_str_pretty(b, "y", vars, ctx);
            fmpz_mpoly_set_str_pretty(c, "1+x+y+z+t", vars, ctx);
            fmpz_mpoly_pow_ui(c, c, i, ctx);
            fmpz_mpoly_add(a, a, c, ctx);
            fmpz_mpoly_add(b, b, c, ctx);
            fmpz_mpoly_mul(a, a, b, ctx);

            k = (i > 0);
            for (j = 1; j <= i; j++)
                if ((j%2) != 0 && (i%j) == 0)
                    k++;
            k = 2;

            time = check_omega(k, a, ctx);
            flint_printf("#%wd: %wd\n", i, time);
            total_time += time;
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }
    flint_printf("total_time: %wd\n", total_time);

    flint_printf("\n------ 5 variables ------\n");
    total_time = 0;
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a;
        const char * vars[] = {"s0", "s1", "s2", "s3", "s4"};
        const slong omegas[] = {2, 4, 1, 1};
        const char * polys[] = {"s1^6*s2^3-6*s0*s1^5*s2^2*s3+12*s0^2*s1^4*s2*s3^2-8*s0^3*s1^3*s3^3+3*s0^2*s1^4*s2^2*s4-12*s0^3*s1^3*s2*s3*s4+12*s0^4*s1^2*s3^2*s4+3*s0^4*s1^2*s2*s4^2-6*s0^5*s1*s3*s4^2+s0^6*s4^3+4*s1^4*s2^2*s3^2-16*s0*s1^3*s2*s3^3+16*s0^2*s1^2*s3^4-4*s1^4*s2^3*s4+16*s0*s1^3*s2^2*s3*s4-8*s0^2*s1^2*s2*s3^2*s4-16*s0^3*s1*s3^3*s4-8*s0^2*s1^2*s2^2*s4^2+16*s0^3*s1*s2*s3*s4^2+4*s0^4*s3^2*s4^2-4*s0^4*s2*s4^3+3*s1^2*s2*s3^4-6*s0*s1*s3^5-6*s1^2*s2^2*s3^2*s4+12*s0*s1*s2*s3^3*s4+3*s0^2*s3^4*s4+3*s1^2*s2^3*s4^2-6*s0*s1*s2^2*s3*s4^2-6*s0^2*s2*s3^2*s4^2+3*s0^2*s2^2*s4^3-2*s3^6+6*s2*s3^4*s4-6*s2^2*s3^2*s4^2+2*s2^3*s4^3",
                              "s3^8-4*s2*s3^6*s4+6*s2^2*s3^4*s4^2-4*s2^3*s3^2*s4^3+s2^4*s4^4",
                              "s1^4*s2^2-4*s0*s1^3*s2*s3+4*s0^2*s1^2*s3^2+2*s0^2*s1^2*s2*s4-4*s0^3*s1*s3*s4+s0^4*s4^2+2*s1^2*s2*s3^2-4*s0*s1*s3^3-2*s1^2*s2^2*s4+4*s0*s1*s2*s3*s4+2*s0^2*s3^2*s4-2*s0^2*s2*s4^2-s3^4+2*s2*s3^2*s4-s2^2*s4^2",
                              "s1^4-2*s1^2*s4-s4^2"};

        fmpz_mpoly_ctx_init(ctx, 5, ORD_LEX);
        fmpz_mpoly_init(a, ctx);

        for (i = 0; i < 4; i++)
        {
            fmpz_mpoly_set_str_pretty(a, polys[i], vars, ctx);
            time = check_omega(omegas[i], a, ctx);
            flint_printf("#%wd: %wd ms\n", i, time);
            total_time += time;
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }
    flint_printf("total_time: %wd\n", total_time);

    flint_cleanup_master();
    return 0;
}

