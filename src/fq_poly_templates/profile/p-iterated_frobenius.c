/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "flint.h"
#include "templates.h"

#include <math.h>
#include "profiler.h"

#define nalgs 2
#define cpumin 2
#define ncases 1

int
main(int argc, char** argv)
{
    fmpz_t p, q;
    int c, n, reps = 0;
    slong d, lenf;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, poly_t) f, *h, finv;
    TEMPLATE(T, mat_t) HH;
    slong i, l;
    double beta;

    double s[nalgs];
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);

    d = atol(argv[2]);
    lenf = atol(argv[3]);

    beta = 0.5 * (1. - (log(2) / log(lenf)));
    l = ceil(pow(lenf, beta));

    TEMPLATE(T, ctx_init)(ctx, p, d, "a");

    TEMPLATE(T, poly_init)(f, ctx);
    TEMPLATE(T, poly_init)(finv, ctx);
    
    fmpz_init(q);
    TEMPLATE(T, ctx_order)(q, ctx);

    if (l < 2)
    {
        printf("l < 2!\n");
    }

    if (!(h = flint_malloc((l + 1) * sizeof(TEMPLATE(T, poly_struct)))))
    {
        printf("Exception (p-iterated_frobenius):\n");
        printf("Not enough memory.\n");
        flint_abort();
    }

    for (i = 0; i < l + 1; i++)
        TEMPLATE(T, poly_init)(h[i], ctx);

    for (c = 0; c < nalgs; c++)
        s[c] = 0.0;
            
    for (n = 0; n < ncases; n++)
    {
        double t[nalgs];
        int lo, loops = 1;

        /*
          Construct random elements of fq
        */
        {
            TEMPLATE(T, poly_randtest_monic)(f, state, lenf, ctx);
            TEMPLATE(T, poly_reverse)(finv, f, f->length, ctx);
            TEMPLATE(T, poly_inv_series_newton)(finv, finv, f->length, ctx);
        }
                
    loop:

        t[0] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_gen)(h[0], ctx);
            TEMPLATE(T, mat_init)(HH, n_sqrt(f->length - 1) + 1, f->length - 1, ctx);
            TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(h[1], h[0], q, 0, f, finv, ctx);
            TEMPLATE(T, poly_precompute_matrix)(HH, h[1], f, finv, ctx);
            for (i = 2; i < l + 1; i++)
                TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)(h[i], h[i - 1],
                                                                        HH, f, finv, ctx);
            TEMPLATE(T, mat_clear)(HH, ctx);
        }
        prof_stop();
        t[0] += get_clock(0);

        
        t[1] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_gen)(h[0], ctx);
            TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(h[1], h[0], q, 0, f, finv, ctx);
            for (i = 2; i < l + 1; i++)
                TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(h[i], h[i-1], q, 0, f, finv, ctx);
        }
        prof_stop();
        t[1] += get_clock(0);

        /*t[2] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_gen)(h[0], ctx);
            TEMPLATE(T, poly_powmod_x_fmpz_preinv)(h[1], q, f, finv, ctx);
            for (i = 2; i < l + 1; i++)
                TEMPLATE(T, poly_compose_mod_preinv)(h[i], h[i-1], h[1], f, finv, ctx);
        }
        prof_stop();
        t[2] += get_clock(0);
               
        t[3] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_gen)(h[0], ctx);
            TEMPLATE(T, poly_powmod_x_fmpz_preinv)(h[1], q, f, finv, ctx);
            for (i = 2; i < l + 1; i++)
                TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(h[i], h[i-1], q, 0, f, finv, ctx);
        }
        prof_stop();
        t[3] += get_clock(0);
        */

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
    
    TEMPLATE(T, poly_clear)(f, ctx);
    TEMPLATE(T, poly_clear)(finv, ctx);
    for (i = 0; i < l + 1; i++)
        TEMPLATE(T, poly_clear)(h[i], ctx);
    flint_free(h);
    TEMPLATE(T, ctx_clear)(ctx);
    fmpz_clear(p);
    fmpz_clear(q);

    FLINT_TEST_CLEANUP(state);
    
    return 0;
}

#endif
