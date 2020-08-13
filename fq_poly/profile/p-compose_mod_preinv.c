/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fq_poly.h"
#include "profiler.h"

#define nalgs 2
#define cpumin 2
#define ncases 1

int
main(int argc, char** argv)
{
    fmpz_t p;
    int c, n, reps = 0;
    slong d, lenf, leng, lenh;
    fq_ctx_t ctx;
    fq_poly_t f, g, h, hinv, res;

    double s[nalgs];
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);

    d = atol(argv[2]);
    lenf = atol(argv[3]);
    leng = atol(argv[4]);
    lenh = atol(argv[5]);

    fq_ctx_init(ctx, p, d, "a");

    fq_poly_init(f, ctx);
    fq_poly_init(g, ctx);
    fq_poly_init(h, ctx);
    fq_poly_init(hinv, ctx);

    fq_poly_init(res, ctx);
    
    for (c = 0; c < nalgs; c++)
        s[c] = 0.0;
            
    for (n = 0; n < ncases; n++)
    {
        double t[nalgs];
        int l, loops = 1;

        /*
          Construct random elements of fq
        */
        {
            fq_poly_randtest_monic(f, state, lenf, ctx);
            fq_poly_randtest_monic(g, state, leng, ctx);
            fq_poly_randtest_monic(h, state, lenh, ctx);
            fq_poly_reverse(hinv, h, h->length, ctx);
            fq_poly_inv_series_newton(hinv, hinv, h->length, ctx);
        }
                
    loop:

        t[0] = 0.0;
        init_clock(0);
        prof_start();
        for (l = 0; l < loops; l++)
            fq_poly_compose_mod_horner_preinv(res, f, g, h, hinv, ctx);
        prof_stop();
        t[0] += get_clock(0);
                
        t[1] = 0.0;
        init_clock(0);
        prof_start();
        for (l = 0; l < loops; l++)
            fq_poly_compose_mod_brent_kung_preinv(res, f, g, h, hinv, ctx);
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

    
    fq_poly_clear(f, ctx);
    fq_poly_clear(g, ctx);
    fq_poly_clear(h, ctx);
    fq_poly_clear(res, ctx);

    fq_ctx_clear(ctx);
    fmpz_clear(p);

    FLINT_TEST_CLEANUP(state);
    
    return 0;
}
