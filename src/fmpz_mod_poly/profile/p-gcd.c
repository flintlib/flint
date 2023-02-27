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
#include <math.h>
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"
#include "profiler.h"

/*
   Definitions for the parameters of the timing process.
   
   lenlo    Minimum length
   lenhi    Maximum length
   lenh     Ratio for the length
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
#define lenhi    2000
#define lenr     1.2
#define bitslo   16
#define bitshi   640
#define bitsr    1.2
#define cols     (slong) (ceil(log((double)(lenhi/lenlo))/log((double)lenr)))
#define rows     (slong) (ceil(log((double)(bitshi/bitslo))/log((double)bitsr)))
#define cpumin   10
#define ncases   1
#define nalgs    2
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
    int i, j, len, bits, maxcols, maxrows;
    int X[rows][cols];
    double T[rows][cols][nalgs];
    fmpz_mod_poly_t f, g, h;
    fmpz_t p;
    fmpz_mod_ctx_t ctx;
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_mod_ctx_init_ui(ctx, 2);
    
    for (len = lenlo, j = 0; len <= lenhi; len = ceil((double)len*lenr), j++)
    {
        slong s[nalgs];
        
        for (bits = bitslo, i = 0; bits <= bitshi; bits = ceil((double)bits*bitsr), i++)
        {
            int c, n, reps = 0;
            
            for (c = 0; c < nalgs; c++)
                s[c] = WORD(0);
            
            do {
               fmpz_randbits(p, state, bits);
            } while (!fmpz_is_probabprime(p));

            fmpz_mod_ctx_set_modulus(ctx, p);
            fmpz_mod_poly_init2(f, len, ctx);
            fmpz_mod_poly_init2(g, len, ctx);
            fmpz_mod_poly_init2(h, len, ctx);
    
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
                        fmpz_randm(f->coeffs + k, state, p);
                        fmpz_randm(g->coeffs + k, state, p);
                    }
                    while ((f->coeffs)[len-1] == WORD(0))
                        fmpz_randm(f->coeffs + (len - 1), state, p);
                    while ((g->coeffs)[len-1] == WORD(0))
                        fmpz_randm(g->coeffs + (len - 1), state, p);
                    f->length = len;
                    g->length = len;
                }
                
              loop:

                timeit_start(t[0]);
                for (l = 0; l < loops; l++)
                    fmpz_mod_poly_gcd_euclidean(h, f, g, ctx);
                timeit_stop(t[0]);
                
                timeit_start(t[1]);
                for (l = 0; l < loops; l++)
                    fmpz_mod_poly_gcd_hgcd(h, f, g, ctx);
                timeit_stop(t[1]);

                for (c = 0; c < nalgs; c++)
                    if (t[c]->cpu <= cpumin)
                    {
                        loops *= 2;
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

            fmpz_mod_poly_clear(f, ctx);
            fmpz_mod_poly_clear(g, ctx);
            fmpz_mod_poly_clear(h, ctx);
        }
        {
           slong sum = 0, c;
           for (c = 0; c < nalgs; c++)
              sum += s[c];
           flint_printf("len = %d, time = %wdms\n", len, sum), fflush(stdout);
        }
    }

    fmpz_clear(p);
    fmpz_mod_ctx_clear(ctx);
    
    maxcols = j;
    maxrows = i;
    
    /* 
       Print 2-D ASCII image of the winning algorithms
     */
    for (i = 0; i < maxrows; i++)
    {
        for (j = 0; j < maxcols; j++)
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
        
        PIXELS = (unsigned char *) flint_malloc(3 * maxrows * maxcols * sizeof(unsigned char));
        k = 0;
        for (i = 0; i < maxrows; i++)
        {
            for (j = 0; j < maxcols; j++)
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

        k = write_rgb_ppm(imgname, PIXELS, maxcols, maxrows);
        flint_free(PIXELS);
        
        if (k)
        {
            flint_printf("Exception:  writing ppm image failed\n");
        }
    }

    flint_randclear(state);
}
