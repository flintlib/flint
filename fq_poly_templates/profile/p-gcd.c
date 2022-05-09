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
#define ncases 50

int
main(int argc, char** argv)
{
    fmpz_t p, q;
    int l, n, reps = 0;
    slong d, len;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, poly_t) a, b, c, g;

    double s[nalgs];
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);
    d = atol(argv[2]);
    len = atol(argv[3]);

    TEMPLATE(T, ctx_init)(ctx, p, d, "a");

    TEMPLATE(T, poly_init)(a, ctx);
    TEMPLATE(T, poly_init)(b, ctx);
    TEMPLATE(T, poly_init)(c, ctx);
    TEMPLATE(T, poly_init)(g, ctx);
    
    fmpz_init(q);
    TEMPLATE(T, ctx_order)(q, ctx);

    for (l = 0; l < nalgs; l++)
        s[l] = 0.0;
            
    for (n = 0; n < ncases; n++)
    {
        double t[nalgs];
        int lo, loops = 1;

        /*
          Construct random elements of fq
        */
        {
            TEMPLATE(T, poly_randtest_monic)(a, state, len / 2, ctx);
            TEMPLATE(T, poly_randtest_monic)(b, state, len / 2, ctx);
            TEMPLATE(T, poly_randtest_monic)(c, state, len / 2, ctx);
            TEMPLATE(T, poly_mul)(a, a, c, ctx);
            TEMPLATE(T, poly_mul)(b, b, c, ctx);
        }
                
    loop:

        t[0] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_gcd_euclidean)(g, a, b, ctx);
        }
        prof_stop();
        t[0] += get_clock(0);

        
        t[1] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_gcd_hgcd)(g, a, b, ctx);
        }
        prof_stop();
        t[1] += get_clock(0);

        for (l = 0; l < nalgs; l++)
            if (t[l] * FLINT_CLOCK_SCALE_FACTOR <= cpumin)
            {
                loops *= 10;
                goto loop;
            }
                
        for (l = 0; l < nalgs; l++)
            s[l] += t[l];
        reps += loops;
    }
            
    for (l = 0; l < nalgs; l++)
    {
        printf("%20f", s[l] / (double) reps);
        fflush(stdout);
    }
    printf("\n");
    
    TEMPLATE(T, poly_clear)(a, ctx);
    TEMPLATE(T, poly_clear)(b, ctx);
    TEMPLATE(T, poly_clear)(c, ctx);
    TEMPLATE(T, poly_clear)(g, ctx);
    TEMPLATE(T, ctx_clear)(ctx);
    fmpz_clear(p);
    fmpz_clear(q);

    FLINT_TEST_CLEANUP(state);
    
    return 0;
}

#endif
