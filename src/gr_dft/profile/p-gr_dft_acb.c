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
#include "acb_dft.h"
#include "gr.h"
#include "gr_dft.h"

/* Compares acb_dft against the drop-in gr_dft_acb, reporting the
   precomputation and per-transform times separately. Optional
   arguments: max_n max_prec. */

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
    acb_dft_pre_struct * apre;
    gr_dft_acb_pre_struct * qpre;
}
args_t;

static void f_acb_pre(void * p)
{
    args_t * a = p;
    acb_dft_pre_t pre;
    acb_dft_precomp_init(pre, a->n, a->prec);
    acb_dft_precomp_clear(pre);
}

static void f_acb_dft(void * p)
{
    args_t * a = p;
    acb_dft_precomp(a->w, a->v, a->apre, a->prec);
}

static void f_gr_pre(void * p)
{
    args_t * a = p;
    gr_dft_acb_pre_t Q;
    GR_IGNORE(gr_dft_acb_precomp_init(Q, a->n, a->prec));
    gr_dft_acb_precomp_clear(Q);
}

static void f_gr_dft(void * p)
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
        flint_printf("%8s %12s %12s %12s %12s %10s\n", "n",
                "acb pre", "acb dft", "gr pre", "gr dft", "dft ratio");

        for (ni = 0; ns[ni] != 0; ni++)
        {
            ulong n = ns[ni];
            slong j;
            acb_ptr v, w;
            acb_dft_pre_t apre;
            gr_dft_acb_pre_t qpre;
            args_t a;
            double t_acb_pre, t_acb_dft, t_gr_pre, t_gr_dft;

            if ((slong) n > max_n)
                continue;

            v = _acb_vec_init(n);
            w = _acb_vec_init(n);
            for (j = 0; j < (slong) n; j++)
                acb_urandom(v + j, state, prec);

            acb_dft_precomp_init(apre, n, prec);
            GR_IGNORE(gr_dft_acb_precomp_init(qpre, n, prec));

            a.w = w;
            a.v = v;
            a.n = n;
            a.prec = prec;
            a.apre = apre;
            a.qpre = qpre;

            t_acb_pre = time_ms(f_acb_pre, &a, 1);
            t_acb_dft = time_ms(f_acb_dft, &a, 1);
            t_gr_pre = time_ms(f_gr_pre, &a, 1);
            t_gr_dft = time_ms(f_gr_dft, &a, 1);

            flint_printf("%8wu %12.4g %12.4g %12.4g %12.4g %9.2fx\n",
                    n, t_acb_pre, t_acb_dft, t_gr_pre, t_gr_dft,
                    t_acb_dft / t_gr_dft);

            acb_dft_precomp_clear(apre);
            gr_dft_acb_precomp_clear(qpre);
            _acb_vec_clear(v, n);
            _acb_vec_clear(w, n);
        }

        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
