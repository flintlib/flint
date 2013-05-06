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
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "profiler.h"

#define lenlo  2
#define lenhi  2
#define lenh   1
#define bits   16
#define elo    26
#define ehi    50
#define eh     1
#define cpumin 100

int
main(void)
{
    int len, e;
    fmpz_poly_t f, g[1];
    
    flint_rand_t state;
    flint_randinit(state);
   
    fmpz_poly_init2(f, lenhi);
    fmpz_poly_init2(g[0], ehi * (lenhi - 1) + 1);
    fmpz_poly_init2(g[1], ehi * (lenhi - 1) + 1);
    
    printf("| len | exp | binomial |\n");
    
    for (len = lenlo; len <= lenhi; len += lenh)
    {
        /*
           Construct random polynomial f of length len
         */
        {
            long k;
            for (k = 0; k < len; k++)
                fmpz_randbits(f->coeffs + k, state, bits);
            if ((f->coeffs)[len-1] == 0L)
                fmpz_randtest_not_zero(f->coeffs + (len-1), state, bits);
            f->length = len;
        }
        
        for (e = elo; e <= ehi; e += eh)
        {
            timeit_t t[1];
            int l, loops = 1, r = 0;
            long s[1] = {0};
            
          loop:
            
            timeit_start(t[0]);
            for (l = 0; l < loops; l++)
                fmpz_poly_pow_binomial(g[0], f, e);
            timeit_stop(t[0]);
            
            if (t[0]->cpu <= cpumin)
            {
                loops *= 10;
                goto loop;
            }
            
            s[0] += t[0]->cpu;
            r    += loops;
            
            printf("%d %d %lf\n", len, e, s[0] / (double) r);
        }
    }
    
    fmpz_poly_clear(f);
    fmpz_poly_clear(g[0]);

    flint_randclear(state);
}
