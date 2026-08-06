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
#include "profiler.h"
#include "acb.h"
#include "gr.h"
#include "gr_dft.h"

/* Compares the two internal paths of the gr_dft_acb drop-in: ball
   arithmetic (a Karatsuba plan over arb/acb contexts, as used by the
   fallback) against the default fixed-point path, reporting the
   precomputation and per-transform times separately. The path column
   shows which arithmetic the automatic plan selected at this length
   and precision (it falls back to ball when the fixed-point limb
   count would be excessive). Optional arguments: max_n max_prec. */

static double
time_ms(void (*f)(void *), void * arg, slong min_reps)
{
    timeit_t tm;
    slong r, reps = min_reps;

    f(arg);   /* warm up */

    for (;;)
    {
        timeit_start(tm);
        for (r = 0; r < reps; r++)
            f(arg);
        timeit_stop(tm);
        if (tm->wall >= 100 || reps >= 100000)
            return (double) tm->wall / reps;
        reps *= 4;
    }
}

typedef struct
{
    acb_ptr w;
    acb_srcptr v;
    slong n;
    slong prec;
    gr_ctx_struct * rctx;
    gr_ctx_struct * cctx;
    gr_dft_pre_struct * bpre;
    gr_dft_acb_pre_struct * qpre;
}
args_t;

/* ball path: contexts and Karatsuba plan, as in the which = 1 branch
   of gr_dft_acb_precomp_init */
static void f_ball_pre(void * p)
{
    args_t * a = p;
    gr_ctx_t rctx, cctx;
    gr_dft_pre_t pre;

    gr_ctx_init_real_arb(rctx, a->prec);
    gr_ctx_init_complex_acb(cctx, a->prec);
    GR_IGNORE(gr_dft_precomp_init_karatsuba(pre, a->n, GR_DFT_ALG_AUTO, 0,
            rctx, cctx));
    gr_dft_precomp_clear(pre);
    gr_ctx_clear(rctx);
    gr_ctx_clear(cctx);
}

static void f_ball_dft(void * p)
{
    args_t * a = p;
    GR_IGNORE(gr_dft_precomp((gr_ptr) a->w, (gr_srcptr) a->v, a->bpre,
            a->cctx));
}

/* fixed-point path: the automatic drop-in plan */
static void f_nfixed_pre(void * p)
{
    args_t * a = p;
    gr_dft_acb_pre_t Q;
    GR_IGNORE(gr_dft_acb_precomp_init(Q, a->n, a->prec));
    gr_dft_acb_precomp_clear(Q);
}

static void f_nfixed_dft(void * p)
{
    args_t * a = p;
    gr_dft_acb_precomp(a->w, a->v, (const gr_dft_acb_pre_struct *) a->qpre,
            a->prec);
}

int
main(int argc, char * argv[])
{
    flint_rand_t state;
    slong max_n = 16384, max_prec = 65536;
    const slong precs[] = { 128, 256, 1024, 4096, 65536, 0 };
    const ulong ns[] = { 64, 100, 101, 256, 720, 729, 1000, 1009, 1024,
        2310, 3125, 4096, 5040, 10000, 10007, 16384, 0 };
    slong pi, ni;

    if (argc > 1) max_n = atol(argv[1]);
    if (argc > 2) max_prec = atol(argv[2]);

    flint_set_num_threads(8);

    flint_rand_init(state);

    flint_printf("wall-clock ms; precomputation and transform reported "
            "separately\n\n");

    for (pi = 0; precs[pi] != 0; pi++)
    {
        slong prec = precs[pi];

        if (prec > max_prec)
            continue;

        flint_printf("prec = %wd\n", prec);
        flint_printf("%8s %12s %12s %12s %12s %10s %7s\n", "n",
                "ball pre", "ball dft", "auto pre", "auto dft", "dft ratio",
                "path");

        for (ni = 0; ns[ni] != 0; ni++)
        {
            ulong n = ns[ni];
            slong j;
            acb_ptr v, w;
            gr_ctx_t rctx, cctx;
            gr_dft_pre_t bpre;
            gr_dft_acb_pre_t qpre;
            args_t a;
            double t_ball_pre, t_ball_dft, t_nf_pre, t_nf_dft;

            if ((slong) n > max_n)
                continue;

            v = _acb_vec_init(n);
            w = _acb_vec_init(n);
            for (j = 0; j < (slong) n; j++)
                acb_urandom(v + j, state, prec);

            gr_ctx_init_real_arb(rctx, prec);
            gr_ctx_init_complex_acb(cctx, prec);
            GR_IGNORE(gr_dft_precomp_init_karatsuba(bpre, n,
                    GR_DFT_ALG_AUTO, 0, rctx, cctx));
            GR_IGNORE(gr_dft_acb_precomp_init(qpre, n, prec));

            a.w = w;
            a.v = v;
            a.n = n;
            a.prec = prec;
            a.rctx = rctx;
            a.cctx = cctx;
            a.bpre = bpre;
            a.qpre = qpre;

            t_ball_pre = time_ms(f_ball_pre, &a, 1);
            t_ball_dft = time_ms(f_ball_dft, &a, 1);
            t_nf_pre = time_ms(f_nfixed_pre, &a, 1);
            t_nf_dft = time_ms(f_nfixed_dft, &a, 1);

            flint_printf("%8wu %12.4g %12.4g %12.4g %12.4g %9.2fx %7s\n",
                    n, t_ball_pre, t_ball_dft, t_nf_pre, t_nf_dft,
                    t_ball_dft / t_nf_dft,
                    (qpre->which == 2) ? "nfixed" : "ball");

            gr_dft_precomp_clear(bpre);
            gr_ctx_clear(rctx);
            gr_ctx_clear(cctx);
            gr_dft_acb_precomp_clear(qpre);
            _acb_vec_clear(v, n);
            _acb_vec_clear(w, n);
        }

        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
