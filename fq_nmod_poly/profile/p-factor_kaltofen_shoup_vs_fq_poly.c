/*
    Copyright (C) 2013 Mike Hansen

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
#include "ulong_extras.h"
#include "profiler.h"

#include "fq_poly.h"
#include "fq_nmod_poly.h"

#define nalgs 2
#define ncases 1
#define cpumin 2

int
main(int argc, char** argv)
{
    slong s[nalgs];

    int c, n, len, ext, reps = 0;
    slong j;
    fmpz_t p, temp;
    fq_poly_t f, g;
    fq_nmod_poly_t fn;
    fq_ctx_t ctx;
    fq_nmod_ctx_t ctxn;
    fmpz_poly_t fpoly;
    nmod_poly_t nmod;
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);

    fmpz_init(temp);
       
    fmpz_set_str(temp, argv[2], 10);
    ext = fmpz_get_si(temp);

    fmpz_set_str(temp, argv[3], 10);
    len = fmpz_get_si(temp);

    fq_ctx_init(ctx, p, ext, "a");

    fmpz_poly_init(fpoly);
    nmod_poly_init(nmod, fmpz_get_ui(p));
    fmpz_poly_get_nmod_poly(nmod, fpoly);
    
    fq_nmod_ctx_init_modulus(ctxn, nmod, "a");
    
    fq_poly_init(f, ctx);
    fq_poly_init(g, ctx);
    fq_nmod_poly_init(fn, ctxn);

    for (c = 0; c < nalgs; c++)
    {
        s[c] = WORD(0);
    }
       
    for (n = 0; n < ncases; n++)
    {
        timeit_t t[nalgs];
        int l, loops = 1;
        fq_poly_factor_t res;
        fq_nmod_poly_factor_t resn;

        /*
           Construct random elements of fq
        */
        {
            fq_poly_randtest_irreducible(f, state, len + 1, ctx);
            fq_poly_randtest_irreducible(g, state, len + 2, ctx);
            fq_poly_mul(f, f, g, ctx);
            fq_poly_make_monic(f, f, ctx);

            fq_nmod_poly_fit_length(fn, f->length, ctxn);
            for (j = 0; j < f->length; j++)
            {
                fmpz_poly_get_nmod_poly(fn->coeffs + j, f->coeffs + j);
                
            }
            _fq_nmod_poly_set_length(fn, f->length, ctxn);
        }
        
    loop:
        fflush(stdout);
        timeit_start(t[0]);
        for (l = 0; l < loops; l++)
        {
            fq_poly_factor_init(res, ctx);
            fq_poly_factor_kaltofen_shoup(res, f, ctx);
            fq_poly_factor_clear(res, ctx);
        }
        timeit_stop(t[0]);

        timeit_start(t[1]);
        for (l = 0; l < loops; l++)
        {            
            fq_nmod_poly_factor_init(resn, ctxn);
            fq_nmod_poly_factor_kaltofen_shoup(resn, fn, ctxn);
            fq_nmod_poly_factor_clear(resn, ctxn);

        }
        timeit_stop(t[1]);

        for (c = 0; c < nalgs; c++)
            if (t[c]->cpu <= cpumin)
            {
                loops += 2;
                goto loop;
            }
        
        for (c = 0; c < nalgs; c++)
            s[c] += t[c]->cpu;
        reps += loops;
    }
        
    for (c = 0; c < nalgs; c++)
    {
        flint_printf("%20f ", s[c] / (double) reps);
        fflush(stdout);
    }
    printf("\n");
        
    fq_poly_clear(f, ctx);
    fq_poly_clear(g, ctx);
    fq_nmod_poly_clear(fn, ctxn);
    nmod_poly_clear(nmod);
    fq_ctx_clear(ctx);
    fq_nmod_ctx_clear(ctxn);
    fmpz_clear(p);
    fmpz_clear(temp);

    FLINT_TEST_CLEANUP(state);
    
    return 0;
}
