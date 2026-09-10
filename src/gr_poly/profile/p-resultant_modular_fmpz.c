/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_poly.h"
#include "profiler.h"

/* random bivariate polynomial of length leny in y whose coefficients have
   length lenx in x and entries of `bits` bits */
static void
_randtest_bivariate(gr_poly_t f, flint_rand_t state, slong leny, slong lenx,
                    flint_bitcnt_t bits, gr_ctx_t ctx, gr_ctx_t cctx)
{
    gr_poly_t c;
    gr_ptr v;
    fmpz_t t;
    slong i, k;

    gr_poly_init(c, cctx);
    v = gr_heap_init(cctx);
    fmpz_init(t);

    GR_MUST_SUCCEED(gr_poly_zero(f, ctx));

    for (i = 0; i < leny; i++)
    {
        GR_MUST_SUCCEED(gr_poly_zero(c, cctx));

        for (k = 0; k < lenx; k++)
        {
            /* nonzero throughout, so that the degrees, and with them the
               number of evaluation points, really are the ones reported */
            fmpz_randtest_not_zero(t, state, bits);
            GR_MUST_SUCCEED(gr_set_fmpz(v, t, cctx));
            GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(c, k, v, cctx));
        }

        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, i, c, ctx));
    }

    fmpz_clear(t);
    gr_heap_clear(v, cctx);
    gr_poly_clear(c, cctx);
}

int main(void)
{
    flint_rand_t state;
    slong leny, lenx;
    flint_bitcnt_t bits;

    flint_rand_init(state);

    flint_printf("res_y(f, g) over Z[x][y], with f and g of length leny in y,\n");
    flint_printf("coefficients of length lenx in x and entries of `bits` bits.\n\n");
    flint_printf("The modular algorithm reduces modulo word-size primes, calls the\n");
    flint_printf("multipoint algorithm on each of them and reconstructs by CRT; the\n");
    flint_printf("subresultant PRS is what the generic resultant uses otherwise.\n");
    flint_printf("`proved` runs to the coefficient bound, `heuristic` stops once the\n");
    flint_printf("reconstruction has been stable over 100 bits of primes.\n\n");

    flint_printf("leny  lenx  bits     subres      proved   ratio"
                 "   heuristic   ratio\n");

    for (leny = 4; leny <= 16; leny *= 2)
    {
        for (lenx = 4; lenx <= 16; lenx *= 2)
        {
            for (bits = 8; bits <= 512; bits *= 8)
            {
                gr_ctx_t cctx, ctx;
                gr_poly_t f, g;
                gr_ptr x, y, z;
                double t1, t2, t3, FLINT_SET_BUT_UNUSED(tt);

                gr_ctx_init_fmpz(cctx);
                gr_ctx_init_gr_poly(ctx, cctx);

                gr_poly_init(f, ctx);
                gr_poly_init(g, ctx);
                x = gr_heap_init(ctx);
                y = gr_heap_init(ctx);
                z = gr_heap_init(ctx);

                _randtest_bivariate(f, state, leny, lenx, bits, ctx, cctx);
                _randtest_bivariate(g, state, leny, lenx, bits, ctx, cctx);

                TIMEIT_START;
                GR_MUST_SUCCEED(gr_poly_resultant_subresultant(y, f, g, ctx));
                TIMEIT_STOP_VALUES(tt, t1);

                TIMEIT_START;
                GR_MUST_SUCCEED(gr_poly_resultant_modular(x, f, g, 1, ctx));
                TIMEIT_STOP_VALUES(tt, t2);

                TIMEIT_START;
                GR_MUST_SUCCEED(gr_poly_resultant_modular(z, f, g, 0, ctx));
                TIMEIT_STOP_VALUES(tt, t3);

                // for theoratical purposes, the proved version is used by default
                // but in practice this test must never fail since the probability 
                // of failing is bounded above by 2^-50
                if (gr_equal(x, y, ctx) == T_FALSE || gr_equal(z, y, ctx) == T_FALSE)
                {
                    flint_printf("\nFAIL: the algorithms disagree\n");
                    flint_abort();
                }

                flint_printf("%4wd  %4wd  %4wu   %8.2es  %8.2es  %6.2fx"
                             "   %8.2es  %6.2fx\n",
                    leny, lenx, bits, t1, t2, t1 / t2, t3, t1 / t3);
                fflush(stdout);

                gr_poly_clear(f, ctx);
                gr_poly_clear(g, ctx);
                gr_heap_clear(x, ctx);
                gr_heap_clear(y, ctx);
                gr_heap_clear(z, ctx);
                gr_ctx_clear(ctx);
                gr_ctx_clear(cctx);
            }
        }
    }

    flint_rand_clear(state);
    return 0;
}
