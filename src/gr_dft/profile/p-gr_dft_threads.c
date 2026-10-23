/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "profiler.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_dft.h"

/* Parallel speedup of the threaded transforms, per algorithm.

   Usage: p-gr_dft_threads [acb|nfixed|nmod] [max_n] [prec] [max_threads] [-v]

   The ring defaults to acb; prec (in bits, for acb and nfixed)
   defaults to 128; max_n to 1048576; max_threads to 8. With -v, each
   plan is printed (gr_dft_precomp_print) before its row. */

static double
time_dft(gr_dft_pre_t P, gr_ptr y, gr_srcptr x, gr_ctx_t ctx)
{
    timeit_t tm;
    slong r, reps = 1;
    int status = GR_SUCCESS;

    status |= gr_dft_precomp(y, x, P, ctx);

    for (;;)
    {
        timeit_start(tm);
        for (r = 0; r < reps; r++)
            status |= gr_dft_precomp(y, x, P, ctx);
        timeit_stop(tm);
        if (tm->wall >= 100 || reps >= 65536)
            return (status == GR_SUCCESS) ? (double) tm->wall / reps : -1.0;
        reps *= 4;
    }
}

/* uniform random data of full working precision (randtest generates
   special values and mixed magnitudes, which misrepresent the cost
   over rings with data-dependent arithmetic) */
static void
rand_vec(gr_ptr x, ulong n, int ring, slong prec, flint_rand_t state,
        gr_ctx_t ctx, gr_ctx_t rctx)
{
    ulong j;

    if (ring == 0)          /* acb */
    {
        acb_ptr v = (acb_ptr) x;
        for (j = 0; j < n; j++)
            acb_urandom(v + j, state, prec);
    }
    else if (ring == 1)     /* nfixed: uniform components scaled below
                               the transform magnitude bound */
    {
        slong rsz = rctx->sizeof_elem;
        slong e = FLINT_BIT_COUNT(n) + 2;
        arb_t u;
        arb_init(u);
        for (j = 0; j < 2 * n; j++)
        {
            gr_ptr comp = (gr_ptr) ((char *) x + j * rsz);
            arb_urandom(u, state, prec + 64);
            arf_mul_2exp_si(arb_midref(u), arb_midref(u), -e);
            GR_IGNORE(gr_dft_nfixed_set_arf(comp, arb_midref(u), rctx));
        }
        arb_clear(u);
    }
    else                    /* nmod: arithmetic cost is data
                               independent; randtest is uniform */
    {
        GR_IGNORE(_gr_vec_randtest(x, state, n, ctx));
    }
}

int
main(int argc, char * argv[])
{
    flint_rand_t state;
    slong max_n = 1048576, prec = 128, max_threads = 8;
    int verbose = 0, ring = 0;
    const int algs[] = { GR_DFT_ALG_CT, GR_DFT_ALG_SPLIT,
        GR_DFT_ALG_BAILEY, GR_DFT_ALG_AUTO };
    const char * alg_names[] = { "ct", "split", "bailey", "auto" };
    ulong n;
    slong ai, nt, i;
    slong nums[3];
    slong num_count = 0;
    gr_ctx_t ctx, rctx;
    int have_rctx = 1;

    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-v") == 0)
            verbose = 1;
        else if (strcmp(argv[i], "acb") == 0)
            ring = 0;
        else if (strcmp(argv[i], "nfixed") == 0)
            ring = 1;
        else if (strcmp(argv[i], "nmod") == 0)
            ring = 2;
        else if (atol(argv[i]) > 0 && num_count < 3)
            nums[num_count++] = atol(argv[i]);
    }

    /* positional numbers: [max_n] [prec] [max_threads] for acb and
       nfixed; [max_n] [max_threads] for nmod (no precision) */
    if (num_count > 0) max_n = nums[0];
    if (ring == 2)
    {
        if (num_count > 1) max_threads = nums[1];
    }
    else
    {
        if (num_count > 1) prec = nums[1];
        if (num_count > 2) max_threads = nums[2];
    }

    flint_rand_init(state);

    if (ring == 0)
    {
        gr_ctx_init_complex_acb(ctx, prec);
        gr_ctx_init_real_arb(rctx, prec);
        flint_printf("acb, prec = %wd", prec);
    }
    else if (ring == 1)
    {
        slong nl = (prec + FLINT_BITS - 1) / FLINT_BITS;
        GR_MUST_SUCCEED(gr_dft_ctx_init_nfixed(rctx, nl));
        GR_MUST_SUCCEED(gr_dft_ctx_init_nfixed_complex(ctx, nl));
        flint_printf("nfixed, %wd limbs (%wd bits)", nl, nl * FLINT_BITS);
    }
    else
    {
        gr_ctx_init_nmod(ctx, UWORD(998244353));
        have_rctx = 0;
        flint_printf("nmod, q = 998244353");
    }
    flint_printf("; per-transform ms by thread count\n\n");

    flint_printf("%8s %8s", "n", "alg");
    for (nt = 1; nt <= max_threads; nt *= 2)
        flint_printf("  %5wd thr", nt);
    flint_printf("   speedup\n");

    for (n = 64; n <= (ulong) max_n; n *= 4)
    {
        for (ai = 0; ai < 4; ai++)
        {
            gr_dft_pre_t P;
            gr_ptr x, y;
            double t1 = 0.0, t = 0.0;
            int status;

            if (ring == 2)
            {
                /* q - 1 = 119 * 2^23: use w = 3^((q-1)/n) */
                gr_ptr w = gr_heap_init(ctx);
                status = gr_set_ui(w, 3, ctx);
                status |= gr_pow_ui(w, w, UWORD(998244352) / n, ctx);
                status |= gr_dft_precomp_init_root(P, w, n, algs[ai], 0, ctx);
                gr_heap_clear(w, ctx);
            }
            else
            {
                status = gr_dft_precomp_init_karatsuba(P, n, algs[ai], 0,
                        rctx, ctx);
            }
            if (status != GR_SUCCESS)
                continue;

            if (verbose)
                gr_dft_precomp_print(P);

            x = gr_heap_init_vec(n, ctx);
            y = gr_heap_init_vec(n, ctx);
            rand_vec(x, n, ring, prec, state, ctx,
                    have_rctx ? rctx : NULL);

            flint_printf("%8wu %8s", n, alg_names[ai]);
            for (nt = 1; nt <= max_threads; nt *= 2)
            {
                flint_set_num_threads(nt);
                t = time_dft(P, y, x, ctx);
                if (nt == 1)
                    t1 = t;
                flint_printf("  %9.3f", t);
            }
            flint_printf("  %7.2fx\n", (t > 0.0) ? t1 / t : 0.0);
            flint_set_num_threads(1);

            gr_heap_clear_vec(x, n, ctx);
            gr_heap_clear_vec(y, n, ctx);
            gr_dft_precomp_clear(P);
        }
    }

    gr_ctx_clear(ctx);
    if (have_rctx)
        gr_ctx_clear(rctx);
    flint_rand_clear(state);
    return 0;
}
