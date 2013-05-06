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
#include <time.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_gmprand(fmpz_poly_t rop, gmp_randstate_t state, long len, long bits)
{
    long i;
    mpz_t c;
    fmpz_t d;
    mpz_init(c);
    fmpz_init(d);
    fmpz_poly_zero(rop);
    fmpz_poly_fit_length(rop, len);
    for (i = 0; i < len; i++)
    {
        mpz_urandomb(c, state, bits);
        if (rand() & 1L)
            mpz_neg(c, c);
        fmpz_set_mpz(d, c);
        fmpz_poly_set_coeff_fmpz(rop, i, d);
    }
    mpz_clear(c);
    fmpz_clear(d);
}

#define lenlo    1
#define lenhi    100
#define lenh     5
#define bitslo   16
#define bitshi   80
#define bitsh    32
#define cols     ((lenhi + 1 - lenlo + (lenh - 1)) / lenh)
#define rows     ((bitshi + 1 - bitslo + (bitsh - 1)) / bitsh)
#define cpumin   10
#define ncases   100

int
main(void)
{
    int i, j, len, bits;
    double X[rows][cols];
    fmpz_poly_t f, g, q, r;
    gmp_randstate_t state;
    
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 362436069UL);
    srand(521288629UL);
        
    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_init(q);
    fmpz_poly_init(r);

    fmpz_poly_fit_length(f, 2 * lenhi - 1);
    fmpz_poly_fit_length(g, lenhi);
    fmpz_poly_fit_length(q, lenhi);
    fmpz_poly_fit_length(r, 2 * lenhi - 1);
    
    printf("3 2 1");
    for (len = lenlo, j = 0; len <= lenhi; len += lenh, j++)
    {
        for (bits = bitslo, i = 0; bits <= bitshi; bits += bitsh, i++)
        {
            int c, n, reps = 0;
            long s = 0L;
            
            for (n = 0; n < ncases; n++)
            {
                clock_t c0, c1;
                int l, loops = 1;
                
                /*
                   Construct random polynomials f and g
                 */
                {
                    fmpz_poly_gmprand(f, state, 2*len - 1, bits);
                    fmpz_poly_gmprand(g, state, len, bits);
                }
                
              loop:

                c0 = clock();
                for (l = 0; l < loops; l++)
                {
                    fmpz_poly_div_divconquer(q, f, g);
                }
                c1 = clock();

                if (c1 - c0 <= cpumin)
                {
                    loops *= 10;
                    goto loop;
                }

                s += c1 - c0;
                reps += loops;
            }

            X[i][j] = (double) s / (double) reps;
            printf(" .");
            fflush(stdout);
        }
    }
    printf("\n");

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(q);
    fmpz_poly_clear(r);
    
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
            printf("%8.3f", X[i][j]);
        printf("\n");
    }

    gmp_randclear(state);
}
