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
#include "fq_nmod_mpoly.h"
#include "profiler.h"


slong count = 0;
slong total_super = 0;

void profile_gcd(
    const fq_nmod_mpoly_t realG,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    fq_nmod_mpoly_t G;
    timeit_t timer;
    slong hensel = -1, brown = -1, zippel = -1, zippel2 = -1, super = -1;

    fq_nmod_mpoly_init(G, ctx);

    if (algo & MPOLY_GCD_USE_BROWN)
    {
        fq_nmod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        fq_nmod_mpoly_gcd_brown(G, A, B, ctx);
        timeit_stop(timer);
        brown = timer->wall;
        if (!fq_nmod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("brown is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", brown);
    fflush(stdout);

    if (algo & MPOLY_GCD_USE_HENSEL)
    {
        fq_nmod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        fq_nmod_mpoly_gcd_hensel(G, A, B, ctx);
        timeit_stop(timer);
        hensel = timer->wall;
        if (!fq_nmod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("hensel is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", hensel);
    fflush(stdout);

    if (algo & MPOLY_GCD_USE_ZIPPEL2)
    {
        fq_nmod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        fq_nmod_mpoly_gcd_zippel2(G, A, B, ctx);
        timeit_stop(timer);
        zippel2 = timer->wall;
        if (!fq_nmod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("zippel2 is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", zippel2);
    fflush(stdout);

    if (algo & MPOLY_GCD_USE_ZIPPEL)
    {
        fq_nmod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        fq_nmod_mpoly_gcd_zippel(G, A, B, ctx);
        timeit_stop(timer);
        zippel = timer->wall;
        if (!fq_nmod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("zippel is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", zippel);
    fflush(stdout);

    {
        fq_nmod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        fq_nmod_mpoly_gcd(G, A, B, ctx);
        timeit_stop(timer);
        super = timer->wall;
        if (!fq_nmod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("super is wrong\n");
            flint_abort();
        }
    }

    count++;

    flint_printf("%10wd     #%wd\n", super, count);
    fflush(stdout);

    total_super += super;

    fq_nmod_mpoly_clear(G, ctx);
}


void print_banner()
{
    flint_printf("|    brown |   hensel |  zippel2 |   zippel |    super |\n");
    flint_printf("+----------+----------+----------+----------+----------+\n");
}

int main(int argc, char *argv[])
{
    slong i;
    const char * vars[] = {"x", "y", "z", "t" ,"u", "v", "w", "s", "p"};
    mp_limb_t p = UWORD(4611686018427388073);

    print_banner();
    for (i = 50; i < 100; i += 4)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, t, m;

        fq_nmod_mpoly_ctx_init_deg(ctx, 2, ORD_LEX, p, 2);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);
        fq_nmod_mpoly_init(m, ctx);

        fq_nmod_mpoly_set_str_pretty(a, "1+#*x+y", vars, ctx);
        fq_nmod_mpoly_pow_ui(a, a, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "x", vars, ctx);
        fq_nmod_mpoly_add(a, a, m, ctx);

        fq_nmod_mpoly_set_str_pretty(b, "#-2*x-#*y", vars, ctx);
        fq_nmod_mpoly_pow_ui(b, b, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "y", vars, ctx);
        fq_nmod_mpoly_add(b, b, m, ctx);

        fq_nmod_mpoly_set_str_pretty(t, "2+#-x+y", vars, ctx);
        fq_nmod_mpoly_pow_ui(t, t, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "x-#*y", vars, ctx);
        fq_nmod_mpoly_add(t, t, m, ctx);

        fq_nmod_mpoly_mul(a, a, t, ctx);
        fq_nmod_mpoly_mul(b, b, t, ctx);
        fq_nmod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_BROWN | MPOLY_GCD_USE_HENSEL);

        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(m, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 3; i < 10; i += 1)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, t, m;

        fq_nmod_mpoly_ctx_init_deg(ctx, 2, ORD_LEX, p, 2);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);
        fq_nmod_mpoly_init(m, ctx);

        fq_nmod_mpoly_set_str_pretty(a, "#-#*x^31+y^51", vars, ctx);
        fq_nmod_mpoly_pow_ui(a, a, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "x^7", vars, ctx);
        fq_nmod_mpoly_add(a, a, m, ctx);

        fq_nmod_mpoly_set_str_pretty(b, "1-2*#*x^23+#*y^47", vars, ctx);
        fq_nmod_mpoly_pow_ui(b, b, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "y^9", vars, ctx);
        fq_nmod_mpoly_add(b, b, m, ctx);

        fq_nmod_mpoly_set_str_pretty(t, "2*#-x^39+#^2*y^24", vars, ctx);
        fq_nmod_mpoly_pow_ui(t, t, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "x^6*y^7", vars, ctx);
        fq_nmod_mpoly_add(t, t, m, ctx);

        fq_nmod_mpoly_mul(a, a, t, ctx);
        fq_nmod_mpoly_mul(b, b, t, ctx);
        fq_nmod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_BROWN | MPOLY_GCD_USE_HENSEL);

        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(m, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 15; i < 25; i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, t, m;

        fq_nmod_mpoly_ctx_init_deg(ctx, 3, ORD_LEX, p, 2);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);
        fq_nmod_mpoly_init(m, ctx);

        fq_nmod_mpoly_set_str_pretty(a, "#+x+#*y+z", vars, ctx);
        fq_nmod_mpoly_pow_ui(a, a, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "x", vars, ctx);
        fq_nmod_mpoly_add(a, a, m, ctx);

        fq_nmod_mpoly_set_str_pretty(b, "1-2*#^2*x-y+#*z", vars, ctx);
        fq_nmod_mpoly_pow_ui(b, b, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "y", vars, ctx);
        fq_nmod_mpoly_add(b, b, m, ctx);

        fq_nmod_mpoly_set_str_pretty(t, "3*#+x+#*y-2*z", vars, ctx);
        fq_nmod_mpoly_pow_ui(t, t, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "z", vars, ctx);
        fq_nmod_mpoly_add(t, t, m, ctx);

        fq_nmod_mpoly_mul(a, a, t, ctx);
        fq_nmod_mpoly_mul(b, b, t, ctx);
        fq_nmod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_ALL);

        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(m, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 1; i < 10; i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, t, m;

        fq_nmod_mpoly_ctx_init_deg(ctx, 7, ORD_LEX, p, 2);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);
        fq_nmod_mpoly_init(m, ctx);

        fq_nmod_mpoly_set_str_pretty(a, "1+#*x+y+z+#*t+#*u+2*v+w", vars, ctx);
        fq_nmod_mpoly_pow_ui(a, a, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "x", vars, ctx);
        fq_nmod_mpoly_add(a, a, m, ctx);

        fq_nmod_mpoly_set_str_pretty(b, "1+x+#*y+z+#*t+u+3*v+#*w", vars, ctx);
        fq_nmod_mpoly_pow_ui(b, b, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "y", vars, ctx);
        fq_nmod_mpoly_add(b, b, m, ctx);

        fq_nmod_mpoly_set_str_pretty(t, "1+#*x+y+#*z+#*t+u+4*#*v+#*w", vars, ctx);
        fq_nmod_mpoly_pow_ui(t, t, i, ctx);
        fq_nmod_mpoly_set_str_pretty(m, "z", vars, ctx);
        fq_nmod_mpoly_add(t, t, m, ctx);

        fq_nmod_mpoly_mul(a, a, t, ctx);
        fq_nmod_mpoly_mul(b, b, t, ctx);
        fq_nmod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_HENSEL | MPOLY_GCD_USE_ZIPPEL2);

        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(m, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

flint_printf("--------------------\ntotal time: %wd\n", total_super);

    flint_cleanup_master();
    return 0;
}

