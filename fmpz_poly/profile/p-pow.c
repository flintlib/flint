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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "profiler.h"

#define lenlo  50
#define lenhi  50
#define lenh   1
#define bits   16
#define elo    10
#define ehi    140
#define eh     10
#define cpumin 100

int
main(void)
{
    int len, e;
    fmpz_poly_t f, g[3];
    
    fmpz_poly_randinit();
    
    fmpz_poly_init2(f, lenhi);
    fmpz_poly_init2(g[0], ehi * (lenhi - 1) + 1);
    fmpz_poly_init2(g[1], ehi * (lenhi - 1) + 1);
    fmpz_poly_init2(g[2], ehi * (lenhi - 1) + 1);
    
    printf("| len | exp | binexp | addchains | multinomial |\n");
    
    for (len = lenlo; len <= lenhi; len += lenh)
    {
        
        /*
           Construct random polynomial f of length len
         */
        {
            long k;
            for (k = 0; k < len; k++)
                fmpz_randbits(f->coeffs + k, bits);
            while (k && !f->coeffs[k - 1])
                k--;
            f->length = k;
        }
        
        for (e = elo; e <= ehi; e += eh)
        {
            timeit_t t[3];
            int l, loops = 1, r = 0;
            long s[3] = {0, 0, 0};
            
          loop:
            
            timeit_start(t[0]);
            for (l = 0; l < loops; l++)
                fmpz_poly_pow_binexp(g[0], f, e);
            timeit_stop(t[0]);

            timeit_start(t[1]);
            for (l = 0; l < loops; l++)
                fmpz_poly_pow_addchains(g[1], f, e);
            timeit_stop(t[1]);
            
            timeit_start(t[2]);
            for (l = 0; l < loops; l++)
                fmpz_poly_pow_multinomial(g[2], f, e);
            timeit_stop(t[2]);
            
            if (t[0]->cpu <= cpumin || t[1]->cpu <= cpumin || t[2]->cpu <= cpumin)
            {
                loops *= 10;
                goto loop;
            }
            
            s[0] += t[0]->cpu;
            s[1] += t[1]->cpu;
            s[2] += t[2]->cpu;
            r    += loops;
            
            printf("%d %d %lf %lf %lf\n", len, e, s[0] / (double) r, s[1] / (double) r, s[2] / (double) r);
        }
    }
    
    fmpz_poly_clear(f);
    fmpz_poly_clear(g[0]);
    fmpz_poly_clear(g[1]);
    fmpz_poly_clear(g[2]);

    fmpz_poly_randclear();
}
