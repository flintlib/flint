/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <mpir.h>

#include "flint.h"
#include "nmod_poly.h"

/* 
    Profiling and benchmarking code for GCD in nmod_poly.

    For three different prime moduli p[i], for a sequence of degrees degs[k], 
    we create 100 random polynomials A, B, C of degree degs[k]/2 and then 
    compute GCD(AC, BC) repeatedly, runs[i][k] times.
 */

int main(void)
{
    flint_rand_t state;

    clock_t c0, c1;
    long double cpu[2];

    mp_limb_t p[] = {17ul, 2147483659ul, 9223372036854775837ul};
    const long degs[]      = {    20,   40,  60,  80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500};
    const long runs[3][25] = {{ 2000, 1000, 500, 300, 200, 200, 200, 180, 140, 140, 100,  80,  80,  80,  50,  50,  40,  30,  30,  20,  18,  16,  14,  12,  10},
                              { 1400,  800, 400, 260, 160, 140, 120, 100,  60,  60,  50,  50,  40,  40,  30,  30,  20,  20,  20,  15,  14,  13,  12,  11,  10},
                              { 1400,  800, 400, 260, 160, 120, 100,  80,  60,  50,  50,  40,  30,  20,  20,  20,  15,  15,  15,  12,  12,  11,  11,  10,  10}};
    long i, k, c, n;

    nmod_poly_t A, B, C, G;

    flint_randinit(state);

    for (i = 0; i < 3; i++)
    {
        printf("---[Modulus %lu]---\n", p[i]), fflush(stdout);

        for (k = 0; k < sizeof(degs)/sizeof(long); k++)
        {
            const long d = degs[k];
            const long r = runs[i][k];

            cpu[0] = 0;
            cpu[1] = 0;

            nmod_poly_init(A, p[i]);
            nmod_poly_init(B, p[i]);
            nmod_poly_init(C, p[i]);
            nmod_poly_init(G, p[i]);

            for (c = 0; c < 100; c++)
            {
                nmod_poly_randtest(A, state, d/2);
                nmod_poly_randtest(B, state, d/2);
                nmod_poly_randtest(C, state, d/2);
                nmod_poly_mul(A, A, C);
                nmod_poly_mul(B, B, C);

                c0 = clock();
                for (n = 0; n < r; n++)
                    nmod_poly_gcd_euclidean(G, A, B);
                c1 = clock();
                cpu[0] += (c1 - c0);

                c0 = clock();
                for (n = 0; n < r; n++)
                    nmod_poly_gcd_hgcd(G, A, B);
                c1 = clock();
                cpu[1] += (c1 - c0);
            }

            cpu[0] = (long double) cpu[0] / (long double) CLOCKS_PER_SEC;
            cpu[1] = (long double) cpu[1] / (long double) CLOCKS_PER_SEC;

            cpu[0] = (long double) cpu[0] / (long double) (100*r);
            cpu[1] = (long double) cpu[1] / (long double) (100*r);

            printf ("%4ld %10.8Lf %10.8Lf\n", A->length, cpu[0], cpu[1]);
            fflush(stdout);

            nmod_poly_clear(A);
            nmod_poly_clear(B);
            nmod_poly_clear(G);
        }
    }

    flint_randclear(state);
    return EXIT_SUCCESS;
}

