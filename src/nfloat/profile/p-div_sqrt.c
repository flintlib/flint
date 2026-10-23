/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Timings (ns per operation) of nfloat division and square roots against
   MPFR at the same precision.  Usage: p-div_sqrt [flags] where flags is
   the nfloat context flags (e.g. 8 for NFLOAT_RND_FLOOR). */

#include <stdlib.h>
#include <mpfr.h>
#include "gr.h"
#include "gr_vec.h"
#include "nfloat.h"
#include "profiler.h"

#define N 256

static double
now(void)
{
    timeit_t t;
    timeit_start_us(t);
    return -1e-6 * (double) t->wall;
}

/* Best of 3 runs of body, each repeated until it takes >= 20 ms; result
   in ns per element of a length-N loop. */
#define TIME(result, body) \
    do { \
        double __best = 1e300; \
        int __k; \
        for (__k = 0; __k < 7; __k++) \
        { \
            slong __reps = 1, __r; \
            double __t; \
            for (;;) \
            { \
                __t = now(); \
                for (__r = 0; __r < __reps; __r++) { body; } \
                __t = now() - __t; \
                if (__t > 0.02) break; \
                __reps *= 2; \
            } \
            __t = __t / __reps; \
            if (__t < __best) __best = __t; \
        } \
        (result) = __best * 1e9 / N; \
    } while (0)

int main(int argc, char * argv[])
{
    slong nlimbs, i;
    int flags = (argc > 1) ? atoi(argv[1]) : 0;
    slong nlist[] = { 1, 2, 3, 4, 5, 6, 8, 12, 16, 24, 32, 48, 66, 0 };
    mpfr_rnd_t rnd;
    int status = GR_SUCCESS;

    if (flags & NFLOAT_RND_FLOOR)
        rnd = MPFR_RNDD;
    else if (flags & NFLOAT_RND_CEIL)
        rnd = MPFR_RNDU;
    else
        rnd = MPFR_RNDZ;

    flint_printf("flags = %d; ns per op, nfloat / mpfr (speedup)\n", flags);
    flint_printf("limbs |        div           |        inv           |        sqrt          |       rsqrt          |    vec_div           |  vec_div_scalar      |  vec_div_scalar_ui\n");

    for (i = 0; nlist[i] != 0; i++)
    {
        gr_ctx_t ctx;
        gr_ptr X, Y, Z, A;
        mpfr_ptr fx, fy, fz, fa;
        flint_rand_t state;
        slong j, prec, sz;
        double t1, t2;
        double tn[7], tm[7];
        ulong c;

        nlimbs = nlist[i];
        prec = nlimbs * FLINT_BITS;
        flint_rand_init(state);

        nfloat_ctx_init(ctx, prec, flags);
        sz = ctx->sizeof_elem;
        X = gr_heap_init_vec(N, ctx);
        Y = gr_heap_init_vec(N, ctx);
        Z = gr_heap_init_vec(N, ctx);
        A = gr_heap_init_vec(N, ctx);

        fx = flint_malloc(N * sizeof(__mpfr_struct));
        fy = flint_malloc(N * sizeof(__mpfr_struct));
        fz = flint_malloc(N * sizeof(__mpfr_struct));
        fa = flint_malloc(N * sizeof(__mpfr_struct));

        for (j = 0; j < N; j++)
        {
            mpfr_init2(fx + j, prec);
            mpfr_init2(fy + j, prec);
            mpfr_init2(fz + j, prec);
            mpfr_init2(fa + j, prec);

            {
                gr_ptr xj = GR_ENTRY(X, j, sz), yj = GR_ENTRY(Y, j, sz);
                slong k;

                for (k = 0; k < nlimbs; k++)
                {
                    NFLOAT_D(xj)[k] = n_randlimb(state);
                    NFLOAT_D(yj)[k] = n_randlimb(state);
                }
                NFLOAT_D(xj)[nlimbs - 1] |= UWORD(1) << (FLINT_BITS - 1);
                NFLOAT_D(yj)[nlimbs - 1] |= UWORD(1) << (FLINT_BITS - 1);
                NFLOAT_EXP(xj) = (slong) n_randint(state, 20) - 10;
                NFLOAT_EXP(yj) = (slong) n_randint(state, 20) - 10;
                NFLOAT_SGNBIT(xj) = n_randint(state, 2);
                NFLOAT_SGNBIT(yj) = n_randint(state, 2);

                /* mpfr_t has the same mantissa layout as nfloat */
                mpfr_set_ui(fx + j, 1, MPFR_RNDN);
                mpfr_set_ui(fy + j, 1, MPFR_RNDN);
                for (k = 0; k < nlimbs; k++)
                {
                    fx[j]._mpfr_d[k] = NFLOAT_D(xj)[k];
                    fy[j]._mpfr_d[k] = NFLOAT_D(yj)[k];
                }
                mpfr_set_exp(fx + j, NFLOAT_EXP(xj));
                mpfr_set_exp(fy + j, NFLOAT_EXP(yj));
                mpfr_setsign(fx + j, fx + j, NFLOAT_SGNBIT(xj), MPFR_RNDN);
                mpfr_setsign(fy + j, fy + j, NFLOAT_SGNBIT(yj), MPFR_RNDN);
                mpfr_abs(fa + j, fx + j, MPFR_RNDN);
            }
            status |= nfloat_abs(GR_ENTRY(A, j, sz), GR_ENTRY(X, j, sz), ctx);
        }

        c = n_randtest_not_zero(state) | 1;

        TIME(tn[0], for (j = 0; j < N; j++) status |= nfloat_div(GR_ENTRY(Z, j, sz), GR_ENTRY(X, j, sz), GR_ENTRY(Y, j, sz), ctx));
        TIME(tm[0], for (j = 0; j < N; j++) mpfr_div(fz + j, fx + j, fy + j, rnd));
        TIME(tn[1], for (j = 0; j < N; j++) status |= nfloat_inv(GR_ENTRY(Z, j, sz), GR_ENTRY(Y, j, sz), ctx));
        TIME(tm[1], for (j = 0; j < N; j++) mpfr_ui_div(fz + j, 1, fy + j, rnd));
        TIME(tn[2], for (j = 0; j < N; j++) status |= nfloat_sqrt(GR_ENTRY(Z, j, sz), GR_ENTRY(A, j, sz), ctx));
        TIME(tm[2], for (j = 0; j < N; j++) mpfr_sqrt(fz + j, fa + j, rnd));
        TIME(tn[3], for (j = 0; j < N; j++) status |= nfloat_rsqrt(GR_ENTRY(Z, j, sz), GR_ENTRY(A, j, sz), ctx));
        TIME(tm[3], for (j = 0; j < N; j++) mpfr_rec_sqrt(fz + j, fa + j, rnd));
        TIME(tn[4], status |= _gr_vec_div(Z, X, Y, N, ctx));
        TIME(tm[4], for (j = 0; j < N; j++) mpfr_div(fz + j, fx + j, fy + j, rnd));
        TIME(tn[5], status |= _gr_vec_div_scalar(Z, X, N, Y, ctx));
        TIME(tm[5], for (j = 0; j < N; j++) mpfr_div(fz + j, fx + j, fy, rnd));
        TIME(tn[6], status |= _gr_vec_div_scalar_ui(Z, X, N, c, ctx));
        TIME(tm[6], for (j = 0; j < N; j++) mpfr_div_ui(fz + j, fx + j, c, rnd));

        flint_printf("%5wd ", nlimbs);
        for (j = 0; j < 7; j++)
            flint_printf("| %6.1f %6.1f (%4.2fx) ", tn[j], tm[j], tm[j] / tn[j]);
        flint_printf("\n");

        (void) t1; (void) t2;

        for (j = 0; j < N; j++)
        {
            mpfr_clear(fx + j);
            mpfr_clear(fy + j);
            mpfr_clear(fz + j);
            mpfr_clear(fa + j);
        }
        flint_free(fx);
        flint_free(fy);
        flint_free(fz);
        flint_free(fa);
        gr_heap_clear_vec(X, N, ctx);
        gr_heap_clear_vec(Y, N, ctx);
        gr_heap_clear_vec(Z, N, ctx);
        gr_heap_clear_vec(A, N, ctx);
        gr_ctx_clear(ctx);
        flint_rand_clear(state);
    }

    if (status != GR_SUCCESS)
        flint_printf("(some operations did not return GR_SUCCESS)\n");

    flint_cleanup_master();
    return 0;
}
