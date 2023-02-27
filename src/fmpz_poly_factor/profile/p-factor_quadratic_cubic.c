/*
    Copyright (C) 2020 Daniel Schultz
    This file is part of FLINT.
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "profiler.h"


static slong count_omega(const fmpz_poly_factor_t fac)
{
    slong i, omega = 0;
    for (i = 0; i < fac->num; i++)
        omega += fac->exp[i];
    return omega;
}


int
main(void)
{
    flint_bitcnt_t max_bits;
    slong it, i, j, omega;
    ulong tab[11][20];
    timeit_t timer;
    FLINT_TEST_INIT(state);

    printf("        bits : ");

    omega = 0;
    for (max_bits = 10, it = 0; max_bits <= 3000; max_bits *= 2, it++)
    {
        fmpz_t d;
        fmpz_poly_t g;
        fmpz_poly_struct * f;
        fmpz_poly_factor_t fac;
        slong reps = 1000000/(max_bits + 100);
        slong outer_reps = 10;

        f = flint_malloc(reps * sizeof(fmpz_poly_struct));
        for (i = 0; i < reps; i++)
            fmpz_poly_init(f + i);
        fmpz_poly_init(g);
        fmpz_poly_factor_init(fac);
        fmpz_init(d);

        tab[0][it] = max_bits;
        flint_printf(" %06wu", tab[0][it]);
        fflush(stdout);

        /* (degree 1)^2 */
        for (i = 0; i < reps; i++)
        {
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_pow(f + i, g, 2);
        }
        timeit_start(timer);
        for (j = 0; j < outer_reps; j++)
        for (i = 0; i < reps; i++)
            fmpz_poly_factor(fac, f + i);
        timeit_stop(timer);
        tab[1][it] = timer->cpu;

        /* (degree 1)*(degree 1) */
        for (i = 0; i < reps; i++)
        {
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_set(f + i, g);
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_mul(f + i, f + i, g);
        }
        timeit_start(timer);
        for (i = 0; i < reps; i++)
        {
            fmpz_poly_factor(fac, f + i);
            omega += count_omega(fac);
        }
        for (j = 0; j < outer_reps; j++)
        for (i = 0; i < reps; i++)
            fmpz_poly_factor(fac, f + i);
        timeit_stop(timer);
        tab[2][it] = timer->cpu;

        /* (degree 2 +) */
        for (i = 0; i < reps; i++)
        {
            do {
               fmpz_poly_randtest(g, state, 3, n_randint(state, max_bits) + 2);
            } while (g->length != 3 || (fmpz_poly_discriminant(d, g), fmpz_sgn(d) < 0));
            fmpz_poly_set(f + i, g);
        }
        timeit_start(timer);
        for (i = 0; i < reps; i++)
        {
            fmpz_poly_factor(fac, f + i);
            omega += count_omega(fac);
        }
        for (j = 0; j < outer_reps; j++)
        for (i = 0; i < reps; i++)
            fmpz_poly_factor(fac, f + i);
        timeit_stop(timer);
        tab[3][it] = timer->cpu;

        /* (degree 2 -) */
        for (i = 0; i < reps; i++)
        {
            do {
               fmpz_poly_randtest(g, state, 3, n_randint(state, max_bits) + 2);
            } while (g->length != 3 || (fmpz_poly_discriminant(d, g), fmpz_sgn(d) > 0));
            fmpz_poly_set(f + i, g);
        }
        timeit_start(timer);
        for (i = 0; i < reps; i++)
        {
            fmpz_poly_factor(fac, f + i);
            omega += count_omega(fac);
        }
        for (j = 0; j < outer_reps; j++)
        for (i = 0; i < reps; i++)
            fmpz_poly_factor(fac, f + i);
        timeit_stop(timer);
        tab[4][it] = timer->cpu;

        /* (degree 1)^3 */
        for (i = 0; i< reps; i++)
        {
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_pow(f + i, g, 3);
        }
        timeit_start(timer);
        for (i = 0; i < reps; i++)
        {
            fmpz_poly_factor(fac, f + i);
            omega += count_omega(fac);
        }
        for (j = 0; j < outer_reps; j++)
        for (i = 0; i < reps; i++)
            fmpz_poly_factor(fac, f + i);
        timeit_stop(timer);
        tab[5][it] = timer->cpu;

        /* (degree 1)*(degree 1)^2 */
        for (i = 0; i < reps; i++)
        {
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_pow(f + i, g, 2);
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_mul(f + i, f + i, g);
        }
        timeit_start(timer);
        for (i = 0; i < reps; i++)
        {
            fmpz_poly_factor(fac, f + i);
            omega += count_omega(fac);
        }
        for (j = 0; j < outer_reps; j++)
        for (i = 0; i < reps; i++)
            fmpz_poly_factor(fac, f + i);
        timeit_stop(timer);
        tab[6][it] = timer->cpu;

        /* (degree 1)*(degree 1)*(degree 1) */
        for (i = 0; i < reps; i++)
        {
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_set(f + i, g);
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_mul(f + i, f + i, g);
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_mul(f + i, f + i, g);
        }
        timeit_start(timer);
        for (i = 0; i < reps; i++)
        {
            fmpz_poly_factor(fac, f + i);
            omega += count_omega(fac);
        }
        for (j = 0; j < outer_reps; j++)
        for (i = 0; i < reps; i++)
            fmpz_poly_factor(fac, f + i);
        timeit_stop(timer);
        tab[7][it] = timer->cpu;

        /* (degree 1)*(degree 2 +) */
        for (i = 0; i < reps; i++)
        {
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_set(f + i, g);
            do {
                fmpz_poly_randtest(g, state, 3, n_randint(state, max_bits) + 2);
            } while (g->length != 3 || (fmpz_poly_discriminant(d, g), fmpz_sgn(d) < 0));
            fmpz_poly_mul(f + i, f + i, g);
        }
        timeit_start(timer);
        for (i = 0; i < reps; i++)
        {
            fmpz_poly_factor(fac, f + i);
            omega += count_omega(fac);
        }
        for (j = 0; j < outer_reps; j++)
        for (i = 0; i < reps; i++)
            fmpz_poly_factor(fac, f + i);
        timeit_stop(timer);
        tab[8][it] = timer->cpu;

        /* (degree 1)*(degree 2 -) */
        for (i = 0; i < reps; i++)
        {
            do {
               fmpz_poly_randtest(g, state, 2, n_randint(state, max_bits) + 2);
            } while (g->length != 2);
            fmpz_poly_set(f + i, g);
            do {
                fmpz_poly_randtest(g, state, 3, n_randint(state, max_bits) + 2);
            } while (g->length != 3 || (fmpz_poly_discriminant(d, g), fmpz_sgn(d) > 0));
            fmpz_poly_mul(f + i, f + i, g);
        }
        timeit_start(timer);
        for (i = 0; i < reps; i++)
        {
            fmpz_poly_factor(fac, f + i);
            omega += count_omega(fac);
        }
        for (j = 0; j < outer_reps; j++)
        for (i = 0; i < reps; i++)
            fmpz_poly_factor(fac, f + i);
        timeit_stop(timer);
        tab[9][it] = timer->cpu;

        /* (degree 3) */
        for (i = 0; i < reps; i++)
        {
            do {
               fmpz_poly_randtest(g, state, 4, n_randint(state, max_bits) + 2);
            } while (g->length != 4);
            fmpz_poly_set(f + i, g);
        }
        timeit_start(timer);
        for (i = 0; i < reps; i++)
        {
            fmpz_poly_factor(fac, f + i);
            omega += count_omega(fac);
        }
        for (j = 0; j < outer_reps; j++)
        for (i = 0; i < reps; i++)
            fmpz_poly_factor(fac, f + i);
        timeit_stop(timer);
        tab[10][it] = timer->cpu;

        for (i = 0; i < reps; i++)
            fmpz_poly_clear(f + i);
        flint_free(f);
        fmpz_poly_clear(g);
        fmpz_poly_factor_clear(fac);
        fmpz_clear(d);
    }

    flint_printf("\n", max_bits);

    printf("---------------");
    for (i = 0; i < it; i++)
        flint_printf("-------");
    printf("\n");

    printf("      (d1)^2 : ");
    for (i = 0; i < it; i++)
        flint_printf(" %06wu", tab[1][i]);
    printf("\n");

    printf("    (d1)(d1) : ");
    for (i = 0; i < it; i++)
        flint_printf(" %06wu", tab[2][i]);
    printf("\n");

    printf("       (d2+) : ");
    for (i = 0; i < it; i++)
        flint_printf(" %06wu", tab[3][i]);
    printf("\n");

    printf("       (d2-) : ");
    for (i = 0; i < it; i++)
        flint_printf(" %06wu", tab[4][i]);
    printf("\n");

    printf("      (d1)^3 : ");
    for (i = 0; i < it; i++)
        flint_printf(" %06wu", tab[5][i]);
    printf("\n");

    printf("  (d1)(d1)^2 : ");
    for (i = 0; i < it; i++)
        flint_printf(" %06wu", tab[6][i]);
    printf("\n");

    printf("(d1)(d1)(d1) : ");
    for (i = 0; i < it; i++)
        flint_printf(" %06wu", tab[7][i]);
    printf("\n");

    printf("   (d1)(d2+) : ");
    for (i = 0; i < it; i++)
        flint_printf(" %06wu", tab[8][i]);
    printf("\n");

    printf("   (d1)(d2-) : ");
    for (i = 0; i < it; i++)
        flint_printf(" %06wu", tab[9][i]);
    printf("\n");

    printf("        (d3) : ");
    for (i = 0; i < it; i++)
        flint_printf(" %06wu", tab[10][i]);
    printf("\n");

    FLINT_TEST_CLEANUP(state);

    flint_printf("factors found: %wd\n", omega);
    return 0;
}
