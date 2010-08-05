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

int
main(void)
{
    long len, lenlo = 1, lenhi = 100, lenh = 1;
    long bits, bitslo = 16, bitshi = 2048, bitsh = 16;
    
    long cols = (lenhi + 1 - lenlo + (lenh - 1)) / lenh;     /* ceil((lenhi + 1 - lenlo) / lenh) */
    long rows = (bitshi + 1 - bitslo + (bitsh - 1)) / bitsh; /* ceil((bitshi + 1 - bitslo) / bitsh) */
    
    int i, j, x[rows][cols];
    long cpumin = 10;
    timeit_t tottime;
    
    timeit_start(tottime);
    
    fmpz_poly_randinit();
    
    for (len = lenlo, j = 0; len <= lenhi; len += lenh, j++)
    {
        for (bits = bitslo, i = 0; bits <= bitshi; bits += bitsh, i++)
        {
            /*
               Run the three methods and find out which is fastest.
             */
            fmpz_poly_t f, g, h;
            timeit_t t[3];
            long l, loops = 1;
            
            fmpz_poly_init2(f, len);
            fmpz_poly_init2(g, len);
            fmpz_poly_init2(h, 2*len - 1);
            
            fmpz_poly_randtest_not_zero(f, len, bits);
            fmpz_poly_randtest_not_zero(g, len, bits);

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
            
            if (t[0]->cpu <= 10 || t[1]->cpu <= 10 || t[2]->cpu <= 10)
            {
                loops *= 10;
                goto loop;
            }
            
            printf("%ld %ld %lf.3 %lf.3 %lf.3\n", len, bits, t[0]->cpu / (double) loops, t[1]->cpu / (double) loops, t[2]->cpu / (double) loops);

            if (t[0]->cpu <= t[1]->cpu && t[0]->cpu <= t[2]->cpu)
                x[i][j] = 0;
            else if (t[1]->cpu <= t[2]->cpu)
                x[i][j] = 1;
            else
                x[i][j] = 2;
            
            fmpz_poly_clear(f);
            fmpz_poly_clear(g);
            fmpz_poly_clear(h);
        }
    }
    
    timeit_stop(tottime);

    printf("Total time: %ldms\n", tottime->cpu);

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
            printf("%d", x[i][j]);
        printf("\n");
    }

    fmpz_poly_randclear();
}
