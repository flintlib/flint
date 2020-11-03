/*
    Copyright 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fq_nmod_mpoly_factor.h"
#include "profiler.h"


slong check_omega(slong om, const fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    fq_nmod_mpoly_t q;
    fq_nmod_mpoly_factor_t g;
    fmpz_t omega;
    timeit_t timer;

flint_printf("---------------------------\n");

    fmpz_init(omega);
    fq_nmod_mpoly_factor_init(g, ctx);
    fq_nmod_mpoly_init(q, ctx);

    timeit_start(timer);
    fq_nmod_mpoly_factor(g, p, ctx);
    timeit_stop(timer);

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, om) < 0)
    {
        flint_printf("factorization has wrong number of factors\n");
        flint_abort();        
    }

    fq_nmod_mpoly_factor_expand(q, g, ctx);
    if (!fq_nmod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        flint_abort();        
    }

    fq_nmod_mpoly_clear(q, ctx);
    fq_nmod_mpoly_factor_clear(g, ctx);
    fmpz_clear(omega);

    return timer->wall;
}


int main(int argc, char *argv[])
{
    slong i, time, total_time = 0;
    mp_limb_t p = UWORD(4611686018427388073);

    flint_printf("starting dense\n");
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, c;
        const char * vars[] = {"x", "y", "z", "t"};

        fq_nmod_mpoly_ctx_init_deg(ctx, 4, ORD_LEX, p, 2);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(c, ctx);

        for (i = 1; i <= 10; i++)
        {
            fq_nmod_mpoly_set_str_pretty(a, "x", vars, ctx);
            fq_nmod_mpoly_set_str_pretty(b, "y", vars, ctx);
            fq_nmod_mpoly_set_str_pretty(c, "#^4+#^3*x+#^2*y+#*z+0*t", vars, ctx);
            fq_nmod_mpoly_pow_ui(c, c, i, ctx);
            fq_nmod_mpoly_add(a, a, c, ctx);
            fq_nmod_mpoly_add(b, b, c, ctx);
            fq_nmod_mpoly_mul(a, a, b, ctx);

            time = check_omega(2, a, ctx);
            flint_printf("power %wd: %wd\n", i, time);
            total_time += time;
        }
        flint_printf("total_time: %wd\n", total_time);

        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(c, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a;
        const char * vars[] = {"x", "y", "z", "t" ,"u", "v", "w"};

    flint_printf("starting sparse\n");

        fq_nmod_mpoly_ctx_init_deg(ctx, 7, ORD_LEX, p, 2);
        fq_nmod_mpoly_init(a, ctx);

        fq_nmod_mpoly_set_str_pretty(a, "(x+(1+#*x+y+z+#*t+u+2*v+w)^7)"
                                       "*(y+(#+x+y+#*z-t+u+v+3*#*w)^7)", vars, ctx);

        time = check_omega(2, a, ctx);
        flint_printf("sparse: %wd\n", time);

        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    flint_cleanup_master();
    return 0;
}

