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

#define lenlo   1
#define lenhi   60
#define lenh    1
#define bitslo  16
#define bitshi  1024
#define bitsh   16
#define cols    ((lenhi + 1 - lenlo + (lenh - 1)) / lenh)
#define rows    ((bitshi + 1 - bitslo + (bitsh - 1)) / bitsh)
#define cpumin  10
#define N       10

int
main(void)
{
    int len, bits;
    int i, j, x[rows][cols];
    
    fmpz_poly_randinit();
    
    for (len = lenlo, j = 0; len <= lenhi; len += lenh, j++)
    {
        fmpz_poly_t f, g, h;
        fmpz_poly_init2(f, len);
        fmpz_poly_init2(g, len);
        fmpz_poly_init2(h, 2*len - 1);
        
        for (bits = bitslo, i = 0; bits <= bitshi; bits += bitsh, i++)
        {
            int n;
            long s[3] = {0, 0, 0};
            int reps = 0;
            
            for (n = 0; n < N; n++)
            {
                timeit_t t[3];
                int l, loops = 1;
                
                /*
                   Construct random polynomials f and g
                 */
                {
                    long k;
                    for (k = 0; k < len; k++)
                    {
                        fmpz_randbits(f->coeffs + k, bits);
                        fmpz_randbits(g->coeffs + k, bits);
                    }
                    k = len;
                    while (k && !f->coeffs[k - 1])
                        k--;
                    f->length = k;
                    k = len;
                    while (k && !g->coeffs[k - 1])
                        k--;
                    g->length = k;
                }
                
              loop:

                timeit_start(t[0]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_mul_classical(h, f, g);
                timeit_stop(t[0]);
                
                timeit_start(t[1]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_mul_karatsuba(h, f, g);
                timeit_stop(t[1]);

                timeit_start(t[2]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_mul_KS(h, f, g);
                timeit_stop(t[2]);
                
                if (t[0]->cpu <= cpumin || t[1]->cpu <= cpumin || t[2]->cpu <= cpumin)
                {
                    loops *= 10;
                    goto loop;
                }
                
                s[0] += t[0]->cpu;
                s[1] += t[1]->cpu;
                s[2] += t[2]->cpu;
                reps += loops;
            }
            
            printf("%d %d %lf %lf %lf\n", len, bits, s[0] / (double) reps, s[1] / (double) reps, s[2] / (double) reps);
            
            if (s[0] <= s[1] && s[0] <= s[2])
                x[i][j] = 0;
            else if (s[1] <= s[2])
                x[i][j] = 1;
            else
                x[i][j] = 2;
        }
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }
    
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
            printf("%d", x[i][j]);
        printf("\n");
    }

    fmpz_poly_randclear();
}
