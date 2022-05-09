/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void fmpz_poly_gmprand(fmpz_poly_t rop, gmp_randstate_t state, slong len, slong bits)
{
    slong i;
    mpz_t c;
    fmpz_t d;
    mpz_init(c);
    fmpz_init(d);
    fmpz_poly_zero(rop);
    fmpz_poly_fit_length(rop, len);
    for (i = 0; i < len; i++)
    {
        mpz_urandomb(c, state, bits);
        if (rand() & WORD(1))
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
    gmp_randseed_ui(state, UWORD(362436069));
    srand(UWORD(521288629));
        
    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_init(q);
    fmpz_poly_init(r);

    fmpz_poly_fit_length(f, 2 * lenhi - 1);
    fmpz_poly_fit_length(g, lenhi);
    fmpz_poly_fit_length(q, lenhi);
    fmpz_poly_fit_length(r, 2 * lenhi - 1);
    
    flint_printf("3 2 1");
    for (len = lenlo, j = 0; len <= lenhi; len += lenh, j++)
    {
        for (bits = bitslo, i = 0; bits <= bitshi; bits += bitsh, i++)
        {
            int n, reps = 0;
            slong s = WORD(0);
            
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
            flint_printf(" .");
            fflush(stdout);
        }
    }
    flint_printf("\n");

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(q);
    fmpz_poly_clear(r);
    
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
            flint_printf("%8.3f", X[i][j]);
        flint_printf("\n");
    }

    gmp_randclear(state);
}
