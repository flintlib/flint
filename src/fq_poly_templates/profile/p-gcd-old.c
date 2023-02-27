/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include "profiler.h"

/* 
    Profiling and benchmarking code for GCD in fq_poly_templates.

    For three different prime moduli p[i], for a sequence of degrees degs[k], 
    we create 100 random polynomials A, B, C of degree degs[k]/2 and then 
    compute GCD(AC, BC) repeatedly, runs[i][k] times.
 */

#define N 50

int main(int argc, char** argv)
{
    FLINT_TEST_INIT(state);

    fmpz_t pp;
    slong ext = atoi(argv[1]);
    
    mp_limb_t p[] = {17ul, 2147483659ul, 9223372036854775837ul};
    const slong degs[]      = {   20,   40,  60,  80, 100, 120, 140, 160, 180, 200, 
                                220,  240, 260, 280, 300, 320, 340, 360, 380, 400, 
                                420,  440, 460, 480, 500, 520, 540, 560, 580, 600, 
                                620,  640, 660, 680, 700, 720, 740, 760, 780, 800, 
                                820,  840, 860, 880, 900, 920, 940, 960, 980, 1000};
    const slong runs[3][N] = {{ 2000, 1000, 500, 300, 200, 200, 200, 180, 140, 140, 
                                100,   80,  80,  80,  50,  50,  40,  30,  30,  20,  
                                 18,   16,  14,  12,  10,  10,  10,  10,  10,  10,  
                                  9,    9,   9,   9,   8,   8,   8,   8,   7,   7,
                                  7,    7,   6,   6,   6,   6,   5,   5,   5,   5},
                             { 1400,  800, 400, 260, 160, 140, 120, 100,  60,  60,  
                                 50,   50,  40,  40,  30,  30,  20,  20,  20,  15,  
                                 14,   13,  12,  11,  10,  10,  10,  10,  10,  10,  
                                  9,    9,   8,   8,   8,   7,   7,   7,   6,   6,
                                  6,    6,   6,   5,   5,   5,   5,   5,   4,   4},
                             { 1400,  800, 400, 260, 160, 120, 100,  80,  60,  50,  
                                 50,   40,  30,  20,  20,  20,  15,  15,  15,  12,  
                                 12,   11,  11,  10,  10,  10,  10,  10,  10,  10,  
                                  9,    9,   8,   8,   8,   7,   7,   7,   6,   6,
                                  6,    6,   6,   5,   5,   5,   5,   5,   4,   4}};

    long double cpu[3][2][N];
    slong i, k, c, n;

    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, poly_t) A, B, C, G;

    fmpz_init(pp);

    for (i = 0; i < 3; i++)
    {
        flint_printf("---[Modulus %wu^%wd]---\n", p[i], ext), fflush(stdout);
        fmpz_set_ui(pp, p[i]);
        TEMPLATE(T, ctx_init)(ctx, pp, ext, "a");

        for (k = 0; k < N; k++)
        {
            const slong d = degs[k];
            slong r = runs[i][k];
            /*r = 5;*/

            cpu[i][0][k] = (long double) 0.0;
            cpu[i][1][k] = (long double) 0.0;

            TEMPLATE(T, poly_init)(A, ctx);
            TEMPLATE(T, poly_init)(B, ctx);
            TEMPLATE(T, poly_init)(C, ctx);
            TEMPLATE(T, poly_init)(G, ctx);

            for (c = 0; c < 100; c++)
            {
                TEMPLATE(T, poly_randtest)(A, state, d/2, ctx);
                TEMPLATE(T, poly_randtest)(B, state, d/2, ctx);
                TEMPLATE(T, poly_randtest)(C, state, d/2, ctx);
                TEMPLATE(T, poly_mul)(A, A, C, ctx);
                TEMPLATE(T, poly_mul)(B, B, C, ctx);
                flint_printf("%4ld %10.8lf %10.8lf\n", A->length, cpu[i][0][k], cpu[i][1][k]);

                init_clock(0);
                prof_start();
                for (n = 0; n < r; n++)
                    TEMPLATE(T, poly_gcd_euclidean)(G, A, B, ctx);
                prof_stop();
                cpu[i][1][k] += get_clock(0);

                init_clock(0);
                prof_start();
                for (n = 0; n < r; n++)
                    TEMPLATE(T, poly_gcd_hgcd)(G, A, B, ctx);
                prof_stop();
                cpu[i][1][k] += get_clock(0);
            }

            cpu[i][0][k] = (long double) cpu[i][0][k] / (long double) CLOCKS_PER_SEC;
            cpu[i][1][k] = (long double) cpu[i][1][k] / (long double) CLOCKS_PER_SEC;

            cpu[i][0][k] = (long double) cpu[i][0][k] / (long double) (100*r);
            cpu[i][1][k] = (long double) cpu[i][1][k] / (long double) (100*r);

            flint_printf("%4ld %10.8llf %10.8llf\n", A->length, cpu[i][0][k], cpu[i][1][k]);
            fflush(stdout);

            TEMPLATE(T, poly_clear)(A, ctx);
            TEMPLATE(T, poly_clear)(B, ctx);
            TEMPLATE(T, poly_clear)(C, ctx);
            TEMPLATE(T, poly_clear)(G, ctx);

            TEMPLATE(T, ctx_clear)(ctx);
        }
    }

    flint_printf("cpu = [");
    for (i = 0; i < 3; i++)
    {
        flint_printf("[[");
        for (k = 0; k < N; k++)
            flint_printf("%.WORD(8)f,", cpu[i][0][k]);
        flint_printf("],");
        flint_printf("[");
        for (k = 0; k < N; k++)
            flint_printf("%.WORD(8)f,", cpu[i][1][k]);
        flint_printf("]],");
    }
    flint_printf("]\n");

    flint_randclear(state);
    return EXIT_SUCCESS;
}



#endif
