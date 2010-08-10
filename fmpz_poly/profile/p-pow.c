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
#include <string.h>
#include <mpir.h>
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
#define nalgs    3
#define img      1
/* #define imgname  "out" */

/*
   Write a binar 24-bit ppm image.
 */
int write_rgb_ppm(const char* file_name, unsigned char* pixels, 
                  unsigned int w, unsigned int h)
{
    FILE* file = fopen(file_name, "wb");
    if (file == NULL)
        return -1;
    fprintf(file, "P6\n%d %d\n255\n", w, h);
    fwrite(pixels, sizeof(unsigned char), w * h * 3, file);
    fclose(file);
    return 0;
}

int
main(void)
{
    int k, exp;
    int X[rows][cols];
    double T[rows][cols][nalgs];
    unsigned char PIXELS[3 * rows * cols * sizeof(unsigned char)];
    
    fmpz_poly_t f, g[3];
    
    fmpz_poly_randinit();
    
    fmpz_poly_init2(f, lenhi);
    fmpz_poly_init2(g[0], exphi * (lenhi - 1) + 1);
    fmpz_poly_init2(g[1], exphi * (lenhi - 1) + 1);
    fmpz_poly_init2(g[2], exphi * (lenhi - 1) + 1);
    
    printf("Comparative timing for fmpz_poly_pow\n");
    printf("\n");
    printf("Length:    [%d..%d] with step size %d\n", lenlo, lenhi, lenh);
    printf("Bit size:  [%d..%d] with step size %d\n", bitslo, bitshi, bitsh);
    printf("Exponents: [%d..%d] with step size %d\n", explo, exphi, exph);
    printf("\n");
    printf("exp len bits (Binary exponentiation) (Addition chains) Multinomials\n");
    printf("\n");
    
    for (exp = explo, k = 0; exp <= exphi; exp += exph, k++)
    {
        int i, j, len, bits;
    
        for (len = lenlo, j = 0; len <= lenhi; len += lenh, j++)
        {
            for (bits = bitslo, i = 0; bits <= bitshi; bits += bitsh, i++)
            {
            
                timeit_t t[nalgs];
                int l, loops = 1, r = 0;
                long s[nalgs] = {0L, 0L, 0L};
            
                /*
                   Construct random polynomial f of length len
                 */
                {
                    long ell;
                    for (ell = 0; ell < len; ell++)
                        fmpz_randbits(f->coeffs + ell, bits);
                    if ((f->coeffs)[len - 1] == 0L)
                        fmpz_randtest_not_zero(f->coeffs + (len - 1), bits);
                    f->length = len;
                }
                
                /*
                   Run tests
                 */
                
              loop:
                
                timeit_start(t[0]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_pow_binexp(g[0], f, exp);
                timeit_stop(t[0]);

                timeit_start(t[1]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_pow_addchains(g[1], f, exp);
                timeit_stop(t[1]);
                
                timeit_start(t[2]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_pow_multinomial(g[2], f, exp);
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
                
                T[i][j][0] = s[0] / (double) r;
                T[i][j][1] = s[1] / (double) r;
                T[i][j][2] = s[2] / (double) r;
                
                if (s[0] <= s[1] && s[0] <= s[2])
                    X[i][j] = 0;
                else if (s[1] <= s[2])
                    X[i][j] = 1;
                else
                    X[i][j] = 2;
            }
        }
        
        /*
           Print raw data
         */
        for (i = 0, len = lenlo; i < rows; i++, len += lenh)
        {
            for (j = 0, bits = bitslo; j < cols; j++, bits += bitsh)
            {
                printf("%d %d %d %f %f %f\n", exp, len, bits, 
                           T[i][j][0], T[i][j][1], T[i][j][2]);
            }
        }

        /*
           Print image
         */
        if (img)
        {
            int c = 0;
            char * fname = (char *) malloc(20 * sizeof(char));
            for (i = 0; i < rows; i++)
            {
                for (j = 0; j < cols; j++)
                {
                    double max = DBL_MIN, v[nalgs];
                    int m;
                    for (m = 0; m < nalgs; m++)
                    {
                        v[m] = T[i][j][m] - T[i][j][X[i][j]];
                        if (v[m] > max)
                            max = v[m];
                    }
                    for (m = 0; m < nalgs; m++)
                    {
                        v[m] = (max - v[m]) / max;
                        PIXELS[c++] = (unsigned char) (v[m] * 255);
                    }
                }
            }
            sprintf(fname, "out_%d.ppm", exp);
            c = write_rgb_ppm(fname, PIXELS, cols, rows);
            free(fname);
        }
        
    }
    
    fmpz_poly_clear(f);
    fmpz_poly_clear(g[0]);
    fmpz_poly_clear(g[1]);
    fmpz_poly_clear(g[2]);

    fmpz_poly_randclear();
}
