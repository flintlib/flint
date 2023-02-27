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
#include "fmpz_mod_mpoly.h"
#include "profiler.h"


slong count = 0;
slong total_super = 0;

void profile_gcd(
    const fmpz_mod_mpoly_t realG,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    fmpz_mod_mpoly_t G;
    timeit_t timer;
    slong hensel = -1, brown = -1, zippel = -1, zippel2 = -1, super = -1;

    fmpz_mod_mpoly_init(G, ctx);

    if (algo & MPOLY_GCD_USE_BROWN)
    {
        fmpz_mod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        fmpz_mod_mpoly_gcd_brown(G, A, B, ctx);
        timeit_stop(timer);
        brown = timer->wall;
        if (!fmpz_mod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("brown is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", brown);
    fflush(stdout);

    if (algo & MPOLY_GCD_USE_HENSEL)
    {
        fmpz_mod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        fmpz_mod_mpoly_gcd_hensel(G, A, B, ctx);
        timeit_stop(timer);
        hensel = timer->wall;
        if (!fmpz_mod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("hensel is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", hensel);
    fflush(stdout);

    if (algo & MPOLY_GCD_USE_ZIPPEL2)
    {
        fmpz_mod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        fmpz_mod_mpoly_gcd_zippel2(G, A, B, ctx);
        timeit_stop(timer);
        zippel2 = timer->wall;
        if (!fmpz_mod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("zippel2 is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", zippel2);
    fflush(stdout);

    if (0 && algo & MPOLY_GCD_USE_ZIPPEL)
    {
        fmpz_mod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        fmpz_mod_mpoly_gcd_zippel(G, A, B, ctx);
        timeit_stop(timer);
        zippel = timer->wall;
        if (!fmpz_mod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("zippel is wrong\n");
            flint_abort();
        }
    }
    flint_printf("%10wd ", zippel);
    fflush(stdout);

    {
        fmpz_mod_mpoly_add(G, A, B, ctx);
        timeit_start(timer);
        fmpz_mod_mpoly_gcd(G, A, B, ctx);
        timeit_stop(timer);
        super = timer->wall;
        if (!fmpz_mod_mpoly_equal(G, realG, ctx))
        {
            flint_printf("super is wrong\n");
            flint_abort();
        }
    }

    count++;

    flint_printf("%10wd     #%wd\n", super, count);
    fflush(stdout);

    total_super += super;

    fmpz_mod_mpoly_clear(G, ctx);
}


void print_banner()
{
    flint_printf("|    brown |   hensel |  zippel2 |   zippel |    super |\n");
    flint_printf("+----------+----------+----------+----------+----------+\n");
}

int main(int argc, char *argv[])
{
    slong i, j;
    const char * vars[] = {"x", "y", "z", "t" ,"u", "v", "w", "s", "p"};
    fmpz_t p;

    fmpz_init(p);

    fmpz_set_ui(p, 3);
    fmpz_pow_ui(p, p, 70);
    fmpz_nextprime(p, p, 1);

    print_banner();
    for (i = 3; i <= 10; i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t g, a, b, t;

        fmpz_mod_mpoly_ctx_init(ctx, i, ORD_DEGREVLEX, p);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        fmpz_mod_mpoly_one(g, ctx);
        fmpz_mod_mpoly_one(a, ctx);
        fmpz_mod_mpoly_one(b, ctx);
        for (j = 0; j < i; j++)
        {
            fmpz_mod_mpoly_gen(t, j, ctx);
            fmpz_mod_mpoly_add_si(t, t, 1, ctx);
            fmpz_mod_mpoly_mul(g, g, t, ctx);
            fmpz_mod_mpoly_gen(t, j, ctx);
            fmpz_mod_mpoly_sub_si(t, t, 2, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_gen(t, j, ctx);
            fmpz_mod_mpoly_add_si(t, t, 2, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);
        }
        fmpz_mod_mpoly_sub_si(g, g, 2, ctx);
        fmpz_mod_mpoly_add_si(a, a, 2, ctx);
        fmpz_mod_mpoly_sub_si(b, b, 2, ctx);

        fmpz_mod_mpoly_mul(a, a, g, ctx);
        fmpz_mod_mpoly_mul(b, b, g, ctx);
        fmpz_mod_mpoly_make_monic(g, g, ctx);

        profile_gcd(g, a, b, ctx, MPOLY_GCD_USE_BROWN | MPOLY_GCD_USE_HENSEL);

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 50; i < 100; i += 4)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, t, m;

        fmpz_mod_mpoly_ctx_init(ctx, 2, ORD_LEX, p);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(m, ctx);

        fmpz_mod_mpoly_set_str_pretty(a, "1+x+y", vars, ctx);
        fmpz_mod_mpoly_pow_ui(a, a, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "x", vars, ctx);
        fmpz_mod_mpoly_add(a, a, m, ctx);

        fmpz_mod_mpoly_set_str_pretty(b, "1-2*x-y", vars, ctx);
        fmpz_mod_mpoly_pow_ui(b, b, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "y", vars, ctx);
        fmpz_mod_mpoly_add(b, b, m, ctx);

        fmpz_mod_mpoly_set_str_pretty(t, "2-x+y", vars, ctx);
        fmpz_mod_mpoly_pow_ui(t, t, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "x-y", vars, ctx);
        fmpz_mod_mpoly_add(t, t, m, ctx);

        fmpz_mod_mpoly_mul(a, a, t, ctx);
        fmpz_mod_mpoly_mul(b, b, t, ctx);
        fmpz_mod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_BROWN | MPOLY_GCD_USE_HENSEL);

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(m, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 5; i < 15; i += 1)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, t, m;

        fmpz_mod_mpoly_ctx_init(ctx, 2, ORD_LEX, p);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(m, ctx);

        fmpz_mod_mpoly_set_str_pretty(a, "1+x^31+y^51", vars, ctx);
        fmpz_mod_mpoly_pow_ui(a, a, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "x^7", vars, ctx);
        fmpz_mod_mpoly_add(a, a, m, ctx);

        fmpz_mod_mpoly_set_str_pretty(b, "1-2*x^23-y^47", vars, ctx);
        fmpz_mod_mpoly_pow_ui(b, b, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "y^9", vars, ctx);
        fmpz_mod_mpoly_add(b, b, m, ctx);

        fmpz_mod_mpoly_set_str_pretty(t, "2-x^39+y^24", vars, ctx);
        fmpz_mod_mpoly_pow_ui(t, t, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "x^6*y^7", vars, ctx);
        fmpz_mod_mpoly_add(t, t, m, ctx);

        fmpz_mod_mpoly_mul(a, a, t, ctx);
        fmpz_mod_mpoly_mul(b, b, t, ctx);
        fmpz_mod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_BROWN | MPOLY_GCD_USE_HENSEL);

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(m, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 15; i < 30; i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, t, m;

        fmpz_mod_mpoly_ctx_init(ctx, 3, ORD_LEX, p);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(m, ctx);

        fmpz_mod_mpoly_set_str_pretty(a, "1+x+y+z", vars, ctx);
        fmpz_mod_mpoly_pow_ui(a, a, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "x", vars, ctx);
        fmpz_mod_mpoly_add(a, a, m, ctx);

        fmpz_mod_mpoly_set_str_pretty(b, "1-2*x-y+z", vars, ctx);
        fmpz_mod_mpoly_pow_ui(b, b, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "y", vars, ctx);
        fmpz_mod_mpoly_add(b, b, m, ctx);

        fmpz_mod_mpoly_set_str_pretty(t, "3+x+y-2*z", vars, ctx);
        fmpz_mod_mpoly_pow_ui(t, t, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "z", vars, ctx);
        fmpz_mod_mpoly_add(t, t, m, ctx);

        fmpz_mod_mpoly_mul(a, a, t, ctx);
        fmpz_mod_mpoly_mul(b, b, t, ctx);
        fmpz_mod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_ALL);

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(m, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 1; i < 7; i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, t, m;

        fmpz_mod_mpoly_ctx_init(ctx, 7, ORD_LEX, p);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(m, ctx);

        fmpz_mod_mpoly_set_str_pretty(a, "1+x+y+z+t+u+2*v+w", vars, ctx);
        fmpz_mod_mpoly_pow_ui(a, a, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "x", vars, ctx);
        fmpz_mod_mpoly_add(a, a, m, ctx);

        fmpz_mod_mpoly_set_str_pretty(b, "1-2*x-y+z+t+u+v+w", vars, ctx);
        fmpz_mod_mpoly_pow_ui(b, b, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "y", vars, ctx);
        fmpz_mod_mpoly_add(b, b, m, ctx);

        fmpz_mod_mpoly_set_str_pretty(t, "1+x+y+z-t-u+v+3*w", vars, ctx);
        fmpz_mod_mpoly_pow_ui(t, t, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "z", vars, ctx);
        fmpz_mod_mpoly_add(t, t, m, ctx);

        fmpz_mod_mpoly_mul(a, a, t, ctx);
        fmpz_mod_mpoly_mul(b, b, t, ctx);
        fmpz_mod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_HENSEL | MPOLY_GCD_USE_ZIPPEL2);

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(m, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    print_banner();
    for (i = 1; i < 5; i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, t, m;

        fmpz_mod_mpoly_ctx_init(ctx, 9, ORD_LEX, p);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(m, ctx);

        fmpz_mod_mpoly_set_str_pretty(a, "1+x+y+z+t+u+2*v+w+s-p", vars, ctx);
        fmpz_mod_mpoly_pow_ui(a, a, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "x", vars, ctx);
        fmpz_mod_mpoly_add(a, a, m, ctx);

        fmpz_mod_mpoly_set_str_pretty(b, "1-2*x-y+z+t+u+v+w-2*s+p", vars, ctx);
        fmpz_mod_mpoly_pow_ui(b, b, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "y", vars, ctx);
        fmpz_mod_mpoly_add(b, b, m, ctx);

        fmpz_mod_mpoly_set_str_pretty(t, "1+x+y+z-t-u+v+3*w+2*s-3*p", vars, ctx);
        fmpz_mod_mpoly_pow_ui(t, t, i, ctx);
        fmpz_mod_mpoly_set_str_pretty(m, "z", vars, ctx);
        fmpz_mod_mpoly_add(t, t, m, ctx);

        fmpz_mod_mpoly_mul(a, a, t, ctx);
        fmpz_mod_mpoly_mul(b, b, t, ctx);
        fmpz_mod_mpoly_make_monic(t, t, ctx);

        profile_gcd(t, a, b, ctx, MPOLY_GCD_USE_HENSEL | MPOLY_GCD_USE_ZIPPEL2);

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(m, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    flint_printf("--------------------\n");
    flint_printf("total time: %wd\n", total_super);
    flint_printf("    record: 25312\n");

    fmpz_clear(p);

    flint_cleanup_master();
    return 0;
}

