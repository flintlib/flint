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
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "profiler.h"
#include "NTL-interface.h"

/*
   Definitions for the parameters of the timing process.
   
   lenlo    Minimum extension size
   lenhi    Maximum extension size
   lenh     Step size for the length
   bitslo   Minimum prime size in bits
   bitshi   Maximum prime size in bitx
   bitsh    Step size for the extension size
   cols     Number of different lengths
   rows     Number of different bit sizes
   cpumin   Minimum number of ms spent on each test case
   ncases   Number of test cases per point (length, bit size)
   nalgs    Number of algorithms
   img      Whether an RGB coloured image should be produced
   imgname  File name for image
 */

#define lenlo    2
#define lenhi    4096
#define lenh     2
#define bitslo   2
#define bitshi   50
#define bitsh    2
#define cols     12 //((lenhi + 1 - lenlo + (lenh - 1)) / lenh)
#define rows     ((bitshi + 1 - bitslo + (bitsh - 1)) / bitsh)
#define cpumin   10
#define ncases   1
#define nalgs    3
#define img      1
#define imgname  "fq-ntl_vs_mul.ppm"

/*
   Write a binary 24-bit ppm image.
 */
int write_rgb_ppm(const char* file_name, unsigned char* pixels, 
                   unsigned int width, unsigned int height)
{
    FILE* file = fopen(file_name, "wb");
    if (file == NULL)
        return -1;
    fprintf(file, "P6\n%d %d\n255\n", width, height);
    fwrite(pixels, sizeof(unsigned char), width * height * 3, file);
    fclose(file);
    return 0;
}

int
main(int argc, char** argv)
{
    int i, j, len, bits;
    int X[rows][cols];
    double T[rows][cols][nalgs];
    fmpz_t p;
    ZZ prime;
    ZZ_pX mod;
    zz_pX nmod;
    fq_t f, g, h;
    fq_ctx_t ctx;
    fmpz_mod_poly_t fmod;
    
    flint_rand_t state;
    flint_randinit(state);

    fq_init(f, ctx);
    fq_init(g, ctx);
    fq_init(h, ctx);
    
    for (bits = bitslo, i = 0; bits <= bitshi; bits += bitsh, i++)
    {
        slong s[nalgs];
        
        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, bits, 1));
        fmpz_get_ZZ(prime, p);
        ZZ_p::init(prime);
        zz_p::init(fmpz_get_si(p));

        flint_printf("Row: p = "); fmpz_print(p); flint_printf("\n");

        for (len = lenlo, j = 0; len <= lenhi; len *= lenh, j++)
        {
            int c, n, reps = 0;

            fmpz_mod_poly_init(fmod, p);
            fmpz_mod_poly_randtest_not_zero(fmod, state, len + 1);
            fmpz_mod_poly_set_coeff_ui(fmod, len, 1);
            fmpz_mod_poly_get_ZZ_pX(mod, fmod);
            fmpz_mod_poly_get_zz_pX(nmod, fmod);

            ZZ_pE::init(mod);
            zz_pE::init(nmod);

            fq_ctx_init_modulus(ctx, p, len, fmod, "a");

            for (c = 0; c < nalgs; c++)
                s[c] = 0L;
            
            for (n = 0; n < ncases; n++)
            {
                timeit_t t[nalgs];
                int l, loops = 1;
                ZZ_pE nf, ng, nh;
                zz_pE nnf, nng, nnh;
                
                /*
                   Construct random elements of fq
                 */
                {
                    fq_randtest_dense(f, state, ctx);
                    fq_get_ZZ_pE(nf, f, ctx);
                    fq_get_zz_pE(nnf, f, ctx);
                    fq_randtest_dense(g, state, ctx);
                    fq_get_ZZ_pE(ng, g, ctx);
                    fq_get_zz_pE(nng, g, ctx);
                }
                
              loop:

                timeit_start(t[0]);
                for (l = 0; l < loops; l++)
                    fq_mul(h, f, g, ctx);
                timeit_stop(t[0]);
                
                timeit_start(t[1]);
                for (l = 0; l < loops; l++)
                    mul(nh, nf, ng);
                timeit_stop(t[1]);

                timeit_start(t[2]);
                for (l = 0; l < loops; l++)
                    mul(nnh, nnf, nng);
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
            
            flint_printf("%8d:", len);
            for (c = 0; c < nalgs; c++)
            {
                T[i][j][c] = s[c] / (double) reps;
                flint_printf("%10f", T[i][j][c]);
                fflush(stdout);
            }
            flint_printf("\n");
            
            if (s[0] <= s[1])
                X[i][j] = 0;
            else
                X[i][j] = 1;
            
            fmpz_mod_poly_clear(fmod);
            fq_ctx_clear(ctx);
        }
    }

    fmpz_clear(p);
    fq_clear(f, ctx);
    fq_clear(g, ctx);
    fq_clear(h, ctx);
    
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
