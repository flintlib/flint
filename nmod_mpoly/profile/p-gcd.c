/*
    Copyright 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "nmod_mpoly.h"
#include "profiler.h"


slong count = 0;
slong total_super = 0;

void profile_gcd(
    const nmod_mpoly_t realG,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    nmod_mpoly_t G;
    timeit_t timer;
    slong hensel = -1, brown = -1, zippel = -1, zippel2 = -1, super = -1;

    nmod_mpoly_init(G, ctx);

    if (algo & MPOLY_GCD_USE_BROWN)
    {
        nmod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        nmod_mpoly_gcd_brown(G, A, B, ctx);
        timeit_stop(timer);
        brown = timer->wall;
        if (!nmod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("brown is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", brown);
    fflush(stdout);

    if (algo & MPOLY_GCD_USE_HENSEL)
    {
        nmod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        nmod_mpoly_gcd_hensel(G, A, B, ctx);
        timeit_stop(timer);
        hensel = timer->wall;
        if (!nmod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("hensel is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", hensel);
    fflush(stdout);

    if (algo & MPOLY_GCD_USE_ZIPPEL2)
    {
        nmod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        nmod_mpoly_gcd_zippel2(G, A, B, ctx);
        timeit_stop(timer);
        zippel2 = timer->wall;
        if (!nmod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("zippel2 is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", zippel2);
    fflush(stdout);

    if (algo & MPOLY_GCD_USE_ZIPPEL)
    {
        nmod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        nmod_mpoly_gcd_zippel(G, A, B, ctx);
        timeit_stop(timer);
        zippel = timer->wall;
        if (!nmod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("zippel is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", zippel);
    fflush(stdout);

    {
        nmod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        nmod_mpoly_gcd(G, A, B, ctx);
        timeit_stop(timer);
        super = timer->wall;
        if (!nmod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("super is wrong\n");
            flint_abort();
        }
    }

    count++;

    flint_printf("%10wd     #%wd\n", super, count);
    fflush(stdout);

    total_super += super;

    nmod_mpoly_clear(G, ctx);
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

    print_banner();
    for (i = 50; i < 100; i += 4)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, t, m;

        nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, 1073741827);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);
        nmod_mpoly_init(m, ctx);

        nmod_mpoly_set_str_pretty(a, "1+x+y", vars, ctx);
        nmod_mpoly_pow_ui(a, a, i, ctx);
        nmod_mpoly_set_str_pretty(m, "x", vars, ctx);
        nmod_mpoly_add(a, a, m, ctx);

        nmod_mpoly_set_str_pretty(b, "1-2*x-y", vars, ctx);
        nmod_mpoly_pow_ui(b, b, i, ctx);
        nmod_mpoly_set_str_pretty(m, "y", vars, ctx);
        nmod_mpoly_add(b, b, m, ctx);

        nmod_mpoly_set_str_pretty(t, "2-x+y", vars, ctx);
        nmod_mpoly_pow_ui(t, t, i, ctx);
        nmod_mpoly_set_str_pretty(m, "x-y", vars, ctx);
        nmod_mpoly_add(t, t, m, ctx);

        nmod_mpoly_mul(a, a, t, ctx);
        nmod_mpoly_mul(b, b, t, ctx);
        nmod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_BROWN | MPOLY_GCD_USE_HENSEL);

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(m, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 5; i < 15; i += 1)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, t, m;

        nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, 1073741827);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);
        nmod_mpoly_init(m, ctx);

        nmod_mpoly_set_str_pretty(a, "1+x^31+y^51", vars, ctx);
        nmod_mpoly_pow_ui(a, a, i, ctx);
        nmod_mpoly_set_str_pretty(m, "x^7", vars, ctx);
        nmod_mpoly_add(a, a, m, ctx);

        nmod_mpoly_set_str_pretty(b, "1-2*x^23-y^47", vars, ctx);
        nmod_mpoly_pow_ui(b, b, i, ctx);
        nmod_mpoly_set_str_pretty(m, "y^9", vars, ctx);
        nmod_mpoly_add(b, b, m, ctx);

        nmod_mpoly_set_str_pretty(t, "2-x^39+y^24", vars, ctx);
        nmod_mpoly_pow_ui(t, t, i, ctx);
        nmod_mpoly_set_str_pretty(m, "x^6*y^7", vars, ctx);
        nmod_mpoly_add(t, t, m, ctx);

        nmod_mpoly_mul(a, a, t, ctx);
        nmod_mpoly_mul(b, b, t, ctx);
        nmod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_BROWN | MPOLY_GCD_USE_HENSEL);

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(m, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 15; i < 30; i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, t, m;

        nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, 1073741827);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);
        nmod_mpoly_init(m, ctx);

        nmod_mpoly_set_str_pretty(a, "1+x+y+z", vars, ctx);
        nmod_mpoly_pow_ui(a, a, i, ctx);
        nmod_mpoly_set_str_pretty(m, "x", vars, ctx);
        nmod_mpoly_add(a, a, m, ctx);

        nmod_mpoly_set_str_pretty(b, "1-2*x-y+z", vars, ctx);
        nmod_mpoly_pow_ui(b, b, i, ctx);
        nmod_mpoly_set_str_pretty(m, "y", vars, ctx);
        nmod_mpoly_add(b, b, m, ctx);

        nmod_mpoly_set_str_pretty(t, "3+x+y-2*z", vars, ctx);
        nmod_mpoly_pow_ui(t, t, i, ctx);
        nmod_mpoly_set_str_pretty(m, "z", vars, ctx);
        nmod_mpoly_add(t, t, m, ctx);

        nmod_mpoly_mul(a, a, t, ctx);
        nmod_mpoly_mul(b, b, t, ctx);
        nmod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_ALL);

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(m, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 1; i < 10; i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, t, m;

        nmod_mpoly_ctx_init(ctx, 7, ORD_LEX, 1073741827);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);
        nmod_mpoly_init(m, ctx);

        nmod_mpoly_set_str_pretty(a, "1+x+y+z+t+u+2*v+w", vars, ctx);
        nmod_mpoly_pow_ui(a, a, i, ctx);
        nmod_mpoly_set_str_pretty(m, "x", vars, ctx);
        nmod_mpoly_add(a, a, m, ctx);

        nmod_mpoly_set_str_pretty(b, "1-2*x-y+z+t+u+v+w", vars, ctx);
        nmod_mpoly_pow_ui(b, b, i, ctx);
        nmod_mpoly_set_str_pretty(m, "y", vars, ctx);
        nmod_mpoly_add(b, b, m, ctx);

        nmod_mpoly_set_str_pretty(t, "1+x+y+z-t-u+v+3*w", vars, ctx);
        nmod_mpoly_pow_ui(t, t, i, ctx);
        nmod_mpoly_set_str_pretty(m, "z", vars, ctx);
        nmod_mpoly_add(t, t, m, ctx);

        nmod_mpoly_mul(a, a, t, ctx);
        nmod_mpoly_mul(b, b, t, ctx);
        nmod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_HENSEL | MPOLY_GCD_USE_ZIPPEL2);

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(m, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 1; i < 8; i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, t, m;

        nmod_mpoly_ctx_init(ctx, 9, ORD_LEX, 1073741827);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);
        nmod_mpoly_init(m, ctx);

        nmod_mpoly_set_str_pretty(a, "1+x+y+z+t+u+2*v+w+s-p", vars, ctx);
        nmod_mpoly_pow_ui(a, a, i, ctx);
        nmod_mpoly_set_str_pretty(m, "x", vars, ctx);
        nmod_mpoly_add(a, a, m, ctx);

        nmod_mpoly_set_str_pretty(b, "1-2*x-y+z+t+u+v+w-2*s+p", vars, ctx);
        nmod_mpoly_pow_ui(b, b, i, ctx);
        nmod_mpoly_set_str_pretty(m, "y", vars, ctx);
        nmod_mpoly_add(b, b, m, ctx);

        nmod_mpoly_set_str_pretty(t, "1+x+y+z-t-u+v+3*w+2*s-3*p", vars, ctx);
        nmod_mpoly_pow_ui(t, t, i, ctx);
        nmod_mpoly_set_str_pretty(m, "z", vars, ctx);
        nmod_mpoly_add(t, t, m, ctx);

        nmod_mpoly_mul(a, a, t, ctx);
        nmod_mpoly_mul(b, b, t, ctx);
        nmod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_HENSEL | MPOLY_GCD_USE_ZIPPEL2);

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(m, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

flint_printf("--------------------\ntotal time: %wd\n", total_super);

    flint_cleanup_master();
    return 0;
}

