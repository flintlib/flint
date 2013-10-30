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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"
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
    long d, lenf;
    fq_ctx_t ctx;
    fq_poly_t f, g, *h, finv, res;
    slong i, l, m, index;
    double beta;

    double s[nalgs];
    
    flint_rand_t state;
    flint_randinit(state);

    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);

    d = atol(argv[2]);
    lenf = atol(argv[3]);

    beta = 0.5 * (1. - (log(2) / log(lenf)));
    l = ceil(pow(lenf, beta));
    m = ceil(0.5 * lenf / l);

    fq_ctx_init(ctx, p, d, "a");

    fq_poly_init(f, ctx);
    fq_poly_init(finv, ctx);
    
    fmpz_init(q);
    fq_ctx_order(q, ctx);

    if (l < 2)
    {
        printf("l < 2!\n");
    }

    if (!(h = flint_malloc((l + 1) * sizeof(fq_poly_struct))))
    {
        printf("Exception (p-iterated_frobenius):\n");
        printf("Not enough memory.\n");
        abort();
    }

    for (i = 0; i < l + 1; i++)
        fq_poly_init(h[i], ctx);

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
            fq_poly_randtest_monic(f, state, lenf, ctx);
            fq_poly_reverse(finv, f, f->length, ctx);
            fq_poly_inv_series_newton(finv, finv, f->length, ctx);
        }
                
    loop:

        t[0] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            fq_poly_gen(h[0], ctx);
            fq_poly_powmod_fmpz_binexp_preinv(h[1], h[0], q, f, finv, ctx);
            for (i = 2; i < l + 1; i++)
                fq_poly_compose_mod_preinv(h[i], h[i-1], h[1], f, finv, ctx);
        }
        prof_stop();
        t[0] += get_clock(0);
                
        t[1] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            fq_poly_gen(h[0], ctx);
            for (i = 1; i < l + 1; i++)
                fq_poly_powmod_fmpz_binexp_preinv(h[i], h[i-1], q, f, finv, ctx);
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
    
    fq_poly_clear(f, ctx);
    fq_poly_clear(finv, ctx);
    for (i = 0; i < l + 1; i++)
        fq_poly_clear(h[i], ctx);
    flint_free(h);
    fq_ctx_clear(ctx);
    fmpz_clear(p);
    fmpz_clear(q);

    return 0;
}
