/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "ulong_extras.h"
#include "gr.h"
#include "gr_poly.h"
#include "profiler.h"

/* the largest prime p = m*2^16 + 1 of `bits` bits, for which the evaluation
   can use a DFT rather than the geometric progression */
static ulong
_fft_prime(int bits)
{
    ulong m;

    for (m = (UWORD(1) << (bits - 16)) - 1; m > (UWORD(1) << (bits - 17)); m--)
    {
        ulong p = (m << 16) + 1;

        if (n_is_prime(p))
            return p;
    }

    return 0;
}

/* random bivariate polynomial of length leny in y whose coefficients have
   length lenx in x */
static void
_randtest_bivariate(gr_poly_t f, flint_rand_t state, slong leny, slong lenx,
                    gr_ctx_t ctx, gr_ctx_t cctx)
{
    gr_poly_t c;
    slong i;

    gr_poly_init(c, cctx);
    GR_MUST_SUCCEED(gr_poly_zero(f, ctx));

    for (i = 0; i < leny; i++)
    {
        GR_MUST_SUCCEED(gr_poly_randtest(c, state, lenx, cctx));

        /* keep the announced degrees, so that the number of evaluation
           points really is the one reported */
        GR_MUST_SUCCEED(gr_poly_set_coeff_ui(c, lenx - 1,
            1 + n_randint(state, 1000), cctx));

        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, i, c, ctx));
    }

    gr_poly_clear(c, cctx);
}

int main(int argc, char * argv[])
{
    flint_rand_t state;
    slong sizes[] = {64, 128, 256, 512};
    slong nsizes = 4;
    slong threads[] = {1, 2, 4, 12};
    slong nthreads = 4;
    slong si, ti, leny, lenx, custom[1];
    int kind;

    /* p-resultant_multipoint_threaded [n], for f and g of length n in y with
       coefficients of length n in x; the sizes below are used otherwise */
    if (argc > 1)
    {
        custom[0] = atol(argv[1]);
        sizes[0] = custom[0];
        nsizes = 1;
    }

    flint_rand_init(state);

    flint_printf("res_y(f, g) over (Z/pZ)[x][y], with f and g of length n in y\n");
    flint_printf("and coefficients of length n in x.\n\n");
    flint_printf("Both the evaluation of the coefficients and the univariate\n");
    flint_printf("resultants at the points are split over the thread pool.\n");
    flint_printf("The two halves of the table are the two evaluations: a modulus\n");
    flint_printf("admitting no DFT, and one admitting one.\n\n");

    flint_printf("                geometric evaluation                DFT evaluation\n");
    flint_printf("    n  npoints    1 thr");
    for (ti = 1; ti < nthreads; ti++)
        flint_printf("  %2wd thr", threads[ti]);
    flint_printf("      1 thr");
    for (ti = 1; ti < nthreads; ti++)
        flint_printf("  %2wd thr", threads[ti]);
    flint_printf("\n");

    for (si = 0; si < nsizes; si++)
    {
        leny = lenx = sizes[si];

        flint_printf("%5wd %8wd", leny,
            (leny - 1) * (lenx - 1) + (leny - 1) * (lenx - 1) + 1);
        fflush(stdout);

        for (kind = 0; kind < 2; kind++)
        {
            gr_ctx_t cctx, ctx;
            gr_poly_t f, g;
            gr_ptr x, y;
            double t1 = 0.0;
            ulong p;

            /* the same size of modulus either way, so that only the
               evaluation the algorithm can use differs */
            if (kind == 0)
                p = n_nextprime(UWORD(1) << 49, 1);
            else
                p = _fft_prime(50);

            gr_ctx_init_nmod(cctx, p);
            gr_ctx_init_gr_poly(ctx, cctx);

            gr_poly_init(f, ctx);
            gr_poly_init(g, ctx);
            x = gr_heap_init(ctx);
            y = gr_heap_init(ctx);

            _randtest_bivariate(f, state, leny, lenx, ctx, cctx);
            _randtest_bivariate(g, state, leny, lenx, ctx, cctx);

            for (ti = 0; ti < nthreads; ti++)
            {
                timeit_t timer;
                double t;

                flint_set_num_threads(threads[ti]);

                timeit_start(timer);
                GR_MUST_SUCCEED(gr_poly_resultant_multipoint(x, f, g, ctx));
                timeit_stop(timer);

                t = timer->wall / 1000.0;

                if (ti == 0)
                {
                    t1 = t;
                    flint_printf("  %7.2fs", t);
                }
                else
                {
                    flint_printf(" %6.2fx", t > 0.0 ? t1 / t : 0.0);
                }

                fflush(stdout);

                /* the answer must not depend on the number of threads */
                if (ti == 0)
                    GR_MUST_SUCCEED(gr_set(y, x, ctx));
                else if (gr_equal(x, y, ctx) == T_FALSE)
                {
                    flint_printf("\nFAIL: the result depends on the thread count\n");
                    flint_abort();
                }
            }

            gr_poly_clear(f, ctx);
            gr_poly_clear(g, ctx);
            gr_heap_clear(x, ctx);
            gr_heap_clear(y, ctx);
            gr_ctx_clear(ctx);
            gr_ctx_clear(cctx);
        }

        flint_printf("\n");
    }

    flint_set_num_threads(1);
    flint_rand_clear(state);
    return 0;
}
