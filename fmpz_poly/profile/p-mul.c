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
#include <gmp.h>
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "profiler.h"

/*
   Definitions for the parameters of the timing process.
   
   lenlo    Minimum length
   lenhi    Maximum length
   lenh     Step size for the length
   bitslo   Minimum bit size
   bitshi   Maximum bit size
   bitsh    Step size for the bit size
   cols     Number of different lengths
   rows     Number of different bit sizes
   cpumin   Minimum number of ms spent on each test case
   ncases   Number of test cases per point (length, bit size)
   nalgs    Number of algorithms
   img      Whether an RGB coloured image should be produced
   imgname  File name for image
 */

#define lenlo    1
#define lenhi    60
#define lenh     1
#define bitslo   16
#define bitshi   2048
#define bitsh    32
#define cols     ((lenhi + 1 - lenlo + (lenh - 1)) / lenh)
#define rows     ((bitshi + 1 - bitslo + (bitsh - 1)) / bitsh)
#define cpumin   10
#define ncases   1
#define nalgs    3
#define img      1
#define imgname  "out.ppm"

/*
   Write a binary 24-bit ppm image.
 */
int write_rgb_ppm(const char* file_name, unsigned char* pixels, 
                   unsigned int width, unsigned int height)
{
    FILE* file = fopen(file_name, "wb");
    if (file == NULL)
        return -1;
    flint_fprintf(file, "P6\n%d %d\n255\n", width, height);
    fwrite(pixels, sizeof(unsigned char), width * height * 3, file);
    fclose(file);
    return 0;
}

int
main(void)
{
    int i, j, len, bits;
    int X[rows][cols];
    double T[rows][cols][nalgs];
    fmpz_poly_t f, g, h;
    
    FLINT_TEST_INIT(state);
    
       
    fmpz_poly_init2(f, lenhi);
    fmpz_poly_init2(g, lenhi);
    fmpz_poly_init2(h, 2*lenhi - 1);
    
    for (len = lenlo, j = 0; len <= lenhi; len += lenh, j++)
    {
        slong s[nalgs];
        
        for (bits = bitslo, i = 0; bits <= bitshi; bits += bitsh, i++)
        {
            int c, n, reps = 0;
            
            for (c = 0; c < nalgs; c++)
                s[c] = WORD(0);
            
            for (n = 0; n < ncases; n++)
            {
                timeit_t t[nalgs];
                int l, loops = 1;
                
                /*
                   Construct random polynomials f and g
                 */
                {
                    slong k;
                    for (k = 0; k < len; k++)
                    {
                        fmpz_randbits(f->coeffs + k, state, bits);
                        fmpz_randbits(g->coeffs + k, state, bits);
                    }
                    if ((f->coeffs)[len-1] == WORD(0))
                        fmpz_randtest_not_zero(f->coeffs + (len - 1), state, bits);
                    if ((g->coeffs)[len-1] == WORD(0))
                        fmpz_randtest_not_zero(g->coeffs + (len - 1), state, bits);
                    f->length = len;
                    g->length = len;
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
            
            if (s[0] <= s[1] && s[0] <= s[2])
                X[i][j] = 0;
            else if (s[1] <= s[2])
                X[i][j] = 1;
            else
                X[i][j] = 2;
        }
        {
           slong sum = 0, c;
           for (c = 0; c < nalgs; c++)
              sum += s[c];
           flint_printf("len = %d, time = %wdms\n", len, sum), fflush(stdout);
        }
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
            flint_printf("%d", X[i][j]);
        flint_printf("\n");
    }
    
    /*
       Print 2-D coloured image to file
     */
    if (img)
    {
        unsigned char * PIXELS;
        int k;
        
        PIXELS = (unsigned char *) flint_malloc(3 * rows * cols * sizeof(unsigned char));
        k = 0;
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
                    PIXELS[k++] = (unsigned char) (v[m] * 255);
                }
                for (; m < 3; m++)
                    PIXELS[k++] = (unsigned char) 0;
            }
        }

        k = write_rgb_ppm(imgname, PIXELS, cols, rows);
        flint_free(PIXELS);
        
        if (k)
        {
            flint_printf("Exception:  writing ppm image failed\n");
        }
    }

    flint_randclear(state);
}
