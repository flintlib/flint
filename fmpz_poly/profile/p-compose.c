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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "profiler.h"

/*
   Definitions for the parameters of the timing process.
   
   len1lo   Minimum length
   len1hi   Maximum length
   len1h    Step size for the length
   len2lo   Minimum length
   len2hi   Maximum length
   len2h    Step size for the length
   bits     Bit size of the coefficients
   cols     Number of different lengths
   rows     Number of different bit sizes
   cpumin   Minimum number of ms spent on each test case
   ncases   Number of test cases per point (length, bit size)
   nalgs    Number of algorithms
   img      Whether an RGB coloured image should be produced
   imgname  File name for image
 */

#define len1lo   1
#define len1hi   30
#define len1h    1
#define len2lo   1
#define len2hi   30
#define len2h    1
#define bits     112
#define cols     ((len1hi + 1 - len1lo + (len1h - 1)) / len1h)
#define rows     ((len2hi + 1 - len2lo + (len2h - 1)) / len2h)
#define cpumin   10
#define ncases   1
#define nalgs    2

int
main(void)
{
    int i, j, len1, len2;
    int X[rows][cols];
    double T[rows][cols][nalgs];
    fmpz_poly_t f, g, h;
    flint_rand_t state;
    flint_randinit(state);
       
    fmpz_poly_init2(f, len1hi);
    fmpz_poly_init2(g, len2hi);
    fmpz_poly_init2(h, (len1hi-1) * (len2hi-1) + 1);
    
    for (len1 = len1lo, j = 0; len1 <= len1hi; len1 += len1h, j++)
    {
        len_t s[nalgs];
        
        for (len2 = len2lo, i = 0; len2 <= len2hi; len2 += len2h, i++)
        {
            int c, n, reps = 0;
            
            for (c = 0; c < nalgs; c++)
                s[c] = 0L;
            
            for (n = 0; n < ncases; n++)
            {
                timeit_t t[nalgs];
                int l, loops = 1;
                
                /*
                   Construct random polynomials f and g
                 */
                {
                    len_t k;
                    for (k = 0; k < len1; k++)
                        fmpz_randbits(f->coeffs + k, state, bits);
                    if ((f->coeffs)[len1-1] == 0L)
                        fmpz_randtest_not_zero(f->coeffs + (len1 - 1), state, bits);
                    f->length = len1;
                }
                {
                    len_t k;
                    for (k = 0; k < len2; k++)
                        fmpz_randbits(g->coeffs + k, state, bits);
                    if ((g->coeffs)[len2-1] == 0L)
                        fmpz_randtest_not_zero(g->coeffs + (len2 - 1), state, bits);
                    g->length = len2;
                }
                
              loop:

                timeit_start(t[0]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_compose_horner(h, f, g);
                timeit_stop(t[0]);
                
                timeit_start(t[1]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_compose_divconquer(h, f, g);
                timeit_stop(t[1]);

                for (c = 0; c < nalgs; c++)
                    if (t[c]->cpu <= cpumin)
                    {
                        loops *= 10;
                        goto loop;
                    }
                
                for (c = 0; c < nalgs; c++)
                    s[c] += t[c]->cpu;
                reps += loops;
            }
            
            for (c = 0; c < nalgs; c++)
                T[i][j][c] = s[c] / (double) reps;
            
            if (s[0] <= s[1])
                X[i][j] = 0;
            else
                X[i][j] = 1;
        }
        printf("len1 = %d, time = %ldms\n", len1, s[0] + s[1]), fflush(stdout);
    }
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(h);
    
    /* 
       Print 2-D ASCII image of the winning algorithms
     */
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
            printf("%d", X[i][j]);
        printf("\n");
    }

    flint_randclear(state);
}
