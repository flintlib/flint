/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fq_nmod_poly.h"
#include "fq_zech_poly.h"
#include "profiler.h"

#define nalgs 2
#define cpumin 2
#define ncases 10

int
main(int argc, char** argv)
{
    fmpz_t p;
    int c, n, reps = 0;
    slong d, len, i;
    fq_nmod_ctx_t ctx;
    fq_zech_ctx_t ctx_zech;
    fq_nmod_poly_t f, g;
    fq_zech_poly_t fz, gz;

    double s[nalgs];
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);

    d = atol(argv[2]);
    len = atol(argv[3]);

    fq_nmod_ctx_init(ctx, p, d, "a");
    fq_zech_ctx_init_fq_nmod_ctx(ctx_zech, ctx);

    fq_nmod_poly_init(f, ctx);
    fq_nmod_poly_init(g, ctx);
    fq_zech_poly_init(fz, ctx_zech);
    fq_zech_poly_init(gz, ctx_zech);

    
    for (c = 0; c < nalgs; c++)
        s[c] = 0.0;
            
    for (n = 0; n < ncases; n++)
    {
        double t[nalgs];
        int l, loops = 1;
        fq_nmod_poly_factor_t res;
        fq_zech_poly_factor_t resz;

        /*
          Construct random elements of fq[x]
        */
        {
            fq_nmod_poly_randtest_monic(f, state, len, ctx);
            fq_zech_poly_fit_length(fz, len, ctx_zech);
            for (i = 0; i < f->length; i++)
                fq_zech_set_fq_nmod(fz->coeffs + i, f->coeffs + i, ctx_zech);
            _fq_zech_poly_set_length(fz, len, ctx_zech);
            _fq_zech_poly_normalise(fz, ctx_zech);
        }
                
    loop:

        t[0] = 0.0;
        init_clock(0);
        prof_start();
        for (l = 0; l < loops; l++)
        {
            fq_nmod_poly_factor_init(res, ctx);
            fq_nmod_poly_factor_kaltofen_shoup(res, f, ctx);
            fq_nmod_poly_factor_clear(res, ctx);
        }
        prof_stop();
        t[0] += get_clock(0);
                
        t[1] = 0.0;
        init_clock(0);
        prof_start();
        for (l = 0; l < loops; l++)
        {
            fq_zech_poly_factor_init(resz, ctx_zech);
            fq_zech_poly_factor_kaltofen_shoup(resz, fz, ctx_zech);
            fq_zech_poly_factor_clear(resz, ctx_zech);
        }
        prof_stop();
        t[1] += get_clock(0);

        for (c = 0; c < nalgs; c++)
            if (t[c] * FLINT_CLOCK_SCALE_FACTOR <= cpumin)
            {
                loops *= 10;
                goto loop;
            }
                
        for (c = 0; c < nalgs; c++)
            s[c] += t[c];
        reps += loops;
    }
            
    for (c = 0; c < nalgs; c++)
    {
        printf("%20f", s[c] / (double) reps);
        fflush(stdout);
    }
    printf("\n");

    
    fq_nmod_poly_clear(f, ctx);
    fq_nmod_poly_clear(g, ctx);
    fq_zech_poly_clear(fz, ctx_zech);
    fq_zech_poly_clear(gz, ctx_zech);

    fq_nmod_ctx_clear(ctx);
    fq_zech_ctx_clear(ctx_zech);
    fmpz_clear(p);

    FLINT_TEST_CLEANUP(state);
    
    return 0;
}
