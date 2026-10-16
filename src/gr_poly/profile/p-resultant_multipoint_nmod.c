/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr.h"
#include "gr_poly.h"
#include "profiler.h"

/* the largest prime p = m*2^20 + 1 of `bits` bits, for which the evaluation
   can use a DFT rather than the geometric progression */
static ulong
_fft_prime(int bits)
{
    ulong m;

    for (m = (UWORD(1) << (bits - 20)) - 1; m > (UWORD(1) << (bits - 21)); m--)
    {
        ulong p = (m << 20) + 1;

        if (n_is_prime(p))
            return p;
    }

    return 0;
}

/* random bivariate polynomial of length leny in y whose coefficients have
   length at most lenx in x */
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

int main(void)
{
    flint_rand_t state;
    slong leny, lenx, npoints;
    int kind;

    flint_rand_init(state);

    flint_printf("res_y(f, g) over (Z/pZ)[x][y], with f and g of length leny\n");
    flint_printf("in y and coefficients of length lenx in x.\n\n");
    flint_printf("The multipoint algorithm evaluates at a geometric progression,\n");
    flint_printf("or at roots of unity when the modulus admits a long enough DFT,\n");
    flint_printf("which is what the two halves of the table compare.\n\n");

    flint_printf("                        p without a DFT              p with a DFT\n");
    flint_printf("leny  lenx  npoints   subres   multipt  ratio      subres   multipt  ratio\n");

    for (leny = 4; leny <= 1024; leny *= 4)
    {
        for (lenx = 4; lenx <= 1024; lenx *= 8)
        {
            npoints = (leny - 1) * (lenx - 1) + (leny - 1) * (lenx - 1) + 1;

            flint_printf("%4wd  %4wd  %7wd", leny, lenx, npoints);
            fflush(stdout);

            for (kind = 0; kind < 2; kind++)
            {
                gr_ctx_t cctx, ctx;
                gr_poly_t f, g;
                gr_ptr x, y;
                double t1, t2, FLINT_SET_BUT_UNUSED(tt);
                ulong p;

                /* the same size of modulus either way, so that the timings
                   differ only in the evaluation that the algorithm can use */
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

                TIMEIT_START;
                GR_MUST_SUCCEED(gr_poly_resultant_subresultant(y, f, g, ctx));
                TIMEIT_STOP_VALUES(tt, t1);

                TIMEIT_START;
                GR_MUST_SUCCEED(gr_poly_resultant_multipoint(x, f, g, ctx));
                TIMEIT_STOP_VALUES(tt, t2);

                if (gr_equal(x, y, ctx) == T_FALSE)
                {
                    flint_printf("\nFAIL: the two algorithms disagree\n");
                    flint_abort();
                }

                flint_printf("   %8.2es %8.2es %6.2fx", t1, t2, t1 / t2);
                fflush(stdout);

                gr_poly_clear(f, ctx);
                gr_poly_clear(g, ctx);
                gr_heap_clear(x, ctx);
                gr_heap_clear(y, ctx);
                gr_ctx_clear(ctx);
                gr_ctx_clear(cctx);
            }

            flint_printf("\n");
        }
    }

    flint_rand_clear(state);
    return 0;
}
