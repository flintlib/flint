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
#include <string.h>
#include <gmp.h>
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "profiler.h"

#define lenlo    1
#define lenhi    100
#define lenh     5
#define bitslo   16
#define bitshi   512
#define bitsh    32
#define explo    2
#define exphi    140
#define exph     5
#define cols     ((lenhi + 1 - lenlo + (lenh - 1)) / lenh)
#define rows     ((bitshi + 1 - bitslo + (bitsh - 1)) / bitsh)
#define height   ((exphi + 1 - explo + (exph - 1)) / exph)
#define cpumin   10
#define img      1

int
main(void)
{
    int k, exp;
    fmpz_poly_t f, g;
    
    FLINT_TEST_INIT(state);
    
   
    fmpz_poly_init2(f, lenhi);
    fmpz_poly_init2(g, exphi * (lenhi - 1) + 1);
    
    flint_printf("Comparative timing for fmpz_poly_pow\n");
    flint_printf("\n");
    flint_printf("Length:    [%d..%d] with step size %d\n", lenlo, lenhi, lenh);
    flint_printf("Bit size:  [%d..%d] with step size %d\n", bitslo, bitshi, bitsh);
    flint_printf("Exponents: [%d..%d] with step size %d\n", explo, exphi, exph);
    flint_printf("\n");
    flint_printf("exp len bits (Binary exponentiation) Multinomials\n");
    flint_printf("\n");
    
    for (exp = explo, k = 0; exp <= exphi; exp += exph, k++)
    {
        int i, j, len, bits;
    
        for (len = lenlo, j = 0; len <= lenhi; len += lenh, j++)
        {
            for (bits = bitslo, i = 0; bits <= bitshi; bits += bitsh, i++)
            {
            
                timeit_t t[2];
                int l, loops = 1, r = 0;
                slong s[2] = {WORD(0), WORD(0)};
                double T[2];
            
                /*
                   Construct random polynomial f of length len
                 */
                {
                    slong ell;
                    for (ell = 0; ell < len; ell++)
                        fmpz_randbits(f->coeffs + ell, state, bits);
                    if ((f->coeffs)[len - 1] == WORD(0))
                        fmpz_randtest_not_zero(f->coeffs + (len - 1), state, bits);
                    f->length = len;
                }
                
                /*
                   Run tests
                 */
                
              loop:
                
                timeit_start(t[0]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_pow_binexp(g, f, exp);
                timeit_stop(t[0]);

                timeit_start(t[1]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_pow_multinomial(g, f, exp);
                timeit_stop(t[1]);
                
                if (t[0]->cpu <= cpumin || t[1]->cpu <= cpumin)
                {
                    loops *= 10;
                    goto loop;
                }
                
                s[0] += t[0]->cpu;
                s[1] += t[1]->cpu;
                r    += loops;
                
                T[0] = s[0] / (double) r;
                T[1] = s[1] / (double) r;
                
                flint_printf("%d %d %d %f %f\n", exp, len, bits, T[0], T[1]);
                fflush(stdout);
            }
        }
    }
    
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);

    flint_randclear(state);  
}
