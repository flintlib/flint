/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 William Hart

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
#include <math.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "profiler.h"

/*
   Definitions for the parameters of the timing process.
   
   lenlo    Minimum length
   lenhi    Maximum length
   lenh     Factor by which length increases at each step
   bitslo   Minimum bit size
   bitshi   Maximum bit size
   bitsh    Factor by which bitsize increases at each step
   cutoff   Maximum product of length*bits
   cols     Number of different lengths
   rows     Number of different bit sizes
   cpumin   Minimum number of ms spent on each test case
   ncases   Number of test cases per point (length, bit size)
   nalgs    Number of algorithms
   img      Whether an RGB coloured image should be produced
   imgname  File name for image
 */

#define lenlo    16
#define lenhi    131072
#define lenh     2
#define bitslo   256
#define bitshi   67108864
#define bitsh    2
#define cutoff   WORD(22000000000)
#define cols     ((slong)(log((double)lenhi/(double)lenlo)/log(lenh) + 1.5))
#define rows     ((slong)(log((double)bitshi/(double)bitslo)/log(bitsh) + 1.5))
#define cpumin   100
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
    slong i, j, k, len, bits;
    int X[rows][cols];
    double T[rows][cols][nalgs];
    fmpz_poly_t f, g, h;
    
    FLINT_TEST_INIT(state);
    
       
    fmpz_poly_init2(f, lenhi);
    fmpz_poly_init2(g, lenhi);
    fmpz_poly_init2(h, 2*lenhi - 1);
    
    for (len = lenlo, j = 0; len <= lenhi; len *= lenh, j++)
    {
        slong s[nalgs], sum;
        
        s[0] = 0;
        s[1] = 0;
        s[2] = 0;
        sum = 0;

        for (bits = bitslo, i = 0; bits <= bitshi; bits *= bitsh, i++)
        {
            int c, n, reps = 0, none = 0;
            
            for (c = 0; c < nalgs; c++)
                s[c] = WORD(0);
            
            for (n = 0; n < ncases; n++)
            {
                timeit_t t[nalgs];
                int l, loops = 1;
                
                if (bits*len <= cutoff)
                {    /*
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
                        fmpz_poly_mul_KS(h, f, g);
                    timeit_stop(t[0]);
                
                    timeit_start(t[1]);
                    for (l = 0; l < loops; l++)
                        fmpz_poly_mul_SS(h, f, g);
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
                } else
                    none = 1;
            }
            
            for (c = 0; c < nalgs; c++)
                T[i][j][c] = s[c] / (double) reps;
            
            if (none)
                X[i][j] = -1;
            else if (s[0] <= s[1] && s[0] <= s[2])
                X[i][j] = 0;
            else if (s[1] <= s[2])
                X[i][j] = 1;
            else
                X[i][j] = 2;

            if (s[0] != 0)
            {
                sum = 0;
                for (c = 0; c < nalgs; c++)
                   sum += s[c];
            }
        }
        flint_printf("len = %wd, time = %wdms\n", len, sum), fflush(stdout);
        for (k = 0; k < rows; k++)
        {
            if (X[k][j] == -1) 
                flint_printf(" ");
            else 
                flint_printf("%d", X[k][j]);
        }
        flint_printf("\n");
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
        {
            if (X[i][j] == -1) 
                flint_printf(" ");
            else 
                flint_printf("%d", X[i][j]);
        }
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
                
                if (X[i][j] == -1)
                {
                    for (m = 0; m < 3; m++)
                       PIXELS[k++] = (unsigned char) 0;
                } else
                {
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
