/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz.h"
#include "gr.h"
#include "gr_special.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "arf.h"
#include "nfloat.h"
#include "profiler.h"
#include "double_extras.h"

#define TABN (NFLOAT_MAX_LIMBS + 1)

#if 1
#undef TIMEIT_END_REPEAT
#define TIMEIT_END_REPEAT(__timer, __reps) \
            } \
            timeit_stop(__timer); \
            if (__timer->cpu >= 100) \
                break; \
            __reps *= 10; \
        } \
    } while (0);
#endif

void
randmat(gr_mat_t mat, flint_rand_t state, gr_ctx_t ctx)
{
    slong m = gr_mat_nrows(mat, ctx);
    slong n = gr_mat_ncols(mat, ctx);

    slong i, j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            gr_ptr v = gr_mat_entry_ptr(mat, i, j, ctx);

            GR_MUST_SUCCEED(gr_set_si(v, 1 + n_randint(state, 1000), ctx));
            GR_MUST_SUCCEED(gr_div_ui(v, v, 1 + n_randint(state, 1000), ctx));
            if (n_randint(state, 2))
                GR_MUST_SUCCEED(gr_neg(v, v, ctx));
        }
    }
}

void tune_classical_vs_fixed(int * cutoffs)
{
    gr_ctx_t ctx;
    gr_mat_t A, B, C;
    slong i, n, nn;
    slong prec;
    double FLINT_SET_BUT_UNUSED(__), t1, t2;

    flint_rand_t state;
    flint_rand_init(state);

    for (i = 0; i < TABN; i++)
        cutoffs[i] = -1;

    for (prec = 64; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += 64)
    {
        flint_printf("prec = %wd\n", prec);

        nfloat_ctx_init(ctx, prec, 0);

        for (nn = 1; nn <= 128; nn++)
        {
            n = (nn == 1) ? 128 : nn;

            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            randmat(A, state, ctx);
            randmat(B, state, ctx);

            TIMEIT_START
            GR_MUST_SUCCEED(gr_mat_mul_classical(C, A, B, ctx));
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul_fixed(C, A, B, 1000, ctx));
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf("%wd  %wd   %e   %e   %.3f\n", prec, n, t1, t2, t1 / t2);

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);

            if (nn == 1)
            {
                if (t2 < t1)
                    continue;
                else
                    break;
            }

            if (t2 < t1)
            {
                cutoffs[prec / 64] = n;

                flint_printf("short tab_classical_vs_fixed[] = {\n");
                for (i = 0; i <= prec / 64; i++)
                    flint_printf("    %d, /* prec = %wd */\n", cutoffs[i], i * 64);
                flint_printf("}\n");

                break;
            }
        }
    }

    flint_rand_clear(state);
}

slong ns[] = { 2, 3, 4, 8, 16, 24, 32, 48, 64, 80, 96, 128, 144, 256, 512, 1024, 0 };

void prof_classical_vs_fixed()
{
    gr_ctx_t ctx;
    gr_mat_t A, B, C;
    slong ni, n;
    slong prec;
    double FLINT_SET_BUT_UNUSED(__), t1, t2;

    flint_rand_t state;
    flint_rand_init(state);

    for (prec = 64; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec = (prec < 1024) ? prec + 64 : prec + 256)
    {
        flint_printf("%wd     ", prec);

        nfloat_ctx_init(ctx, prec, 0);

        for (ni = 8; (n = ns[ni]) != 0; ni++)
        {
            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            randmat(A, state, ctx);
            randmat(B, state, ctx);

            TIMEIT_START
            GR_MUST_SUCCEED(gr_mat_mul_classical(C, A, B, ctx));
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul_fixed(C, A, B, 1000, ctx));
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf("%.3f  ", t1 / t2);
            fflush(stdout);

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);
        }

        flint_printf("\n");
    }

    flint_rand_clear(state);
}

void tune_fixed_vs_block(int * cutoffs)
{
    gr_ctx_t ctx;
    gr_mat_t A, B, C;
    slong i, n;
    slong prec;
    double FLINT_SET_BUT_UNUSED(__), t1, t2;

    flint_rand_t state;
    flint_rand_init(state);

    for (i = 0; i < TABN; i++)
        cutoffs[i] = -1;

    for (prec = 64; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += 64)
    {
        flint_printf("prec = %wd\n", prec);

        nfloat_ctx_init(ctx, prec, 0);

        for (n = 16; ; n = FLINT_MAX(n + 1, n * 1.1))
        {
            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            randmat(A, state, ctx);
            randmat(B, state, ctx);

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul_fixed(C, A, B, 1000, ctx));
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul_block(C, A, B, 1, ctx));
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf("%wd  %wd   %e   %e   %.3f\n", prec, n, t1, t2, t1 / t2);

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);

            if (t2 < t1)
            {
                cutoffs[prec / 64] = n;

                flint_printf("short tab_fixed_vs_block[] = {\n");
                for (i = 0; i <= prec / 64; i++)
                    flint_printf("    %d, /* prec = %wd */\n", cutoffs[i], i * 64);
                flint_printf("}\n");

                break;
            }
        }
    }

    flint_rand_clear(state);
}

void prof_fixed_vs_block()
{
    gr_ctx_t ctx;
    gr_mat_t A, B, C;
    slong ni, n;
    slong prec;
    double FLINT_SET_BUT_UNUSED(__), t1, t2;

    flint_rand_t state;
    flint_rand_init(state);

    flint_printf("        ");
    for (ni = 8; (n = ns[ni]) != 0; ni++)
        flint_printf("%5wd  ", n);
    flint_printf("\n");

    for (prec = 64; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec = (prec < 1024) ? prec + 64 : prec + 256)
    {
        flint_printf("%4wd     ", prec);

        nfloat_ctx_init(ctx, prec, 0);

        for (ni = 8; (n = ns[ni]) != 0; ni++)
        {
            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            randmat(A, state, ctx);
            randmat(B, state, ctx);

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul_fixed(C, A, B, 1000, ctx));
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul_block(C, A, B, 1, ctx));
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf("%.3f  ", t1 / t2);
            fflush(stdout);

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);
        }

        flint_printf("\n");
    }

    flint_rand_clear(state);
}

void tune_classical_vs_block(int * cutoffs)
{
    gr_ctx_t ctx;
    gr_mat_t A, B, C;
    slong i, n;
    slong prec;
    double FLINT_SET_BUT_UNUSED(__), t1, t2;

    flint_rand_t state;
    flint_rand_init(state);

    for (i = 0; i < TABN; i++)
        cutoffs[i] = -1;

    for (prec = 64; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += 64)
    {
        flint_printf("prec = %wd\n", prec);

        nfloat_ctx_init(ctx, prec, 0);

        for (n = 16; ; n = FLINT_MAX(n + 1, n * 1.1))
        {
            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            randmat(A, state, ctx);
            randmat(B, state, ctx);

            TIMEIT_START
            GR_MUST_SUCCEED(gr_mat_mul_classical(C, A, B, ctx));
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul_block(C, A, B, 1, ctx));
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf("%wd  %wd   %e   %e   %.3f\n", prec, n, t1, t2, t1 / t2);

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);

            if (t2 < t1)
            {
                cutoffs[prec / 64] = n;

                flint_printf("short tab_classical_vs_block[] = {\n");
                for (i = 0; i <= prec / 64; i++)
                    flint_printf("    %d, /* prec = %wd */\n", cutoffs[i], i * 64);
                flint_printf("}\n");

                break;
            }
        }
    }

    flint_rand_clear(state);
}

void prof_mul()
{
    gr_ctx_t ctx;
    gr_mat_t A, B, C;
    slong ni, n;
    slong prec;
    double FLINT_SET_BUT_UNUSED(__), t1, t2;

    flint_rand_t state;
    flint_rand_init(state);

    flint_printf("        ");
    for (ni = 0; (n = ns[ni]) != 0; ni++)
        flint_printf("%5wd ", n);
    flint_printf("\n");

    for (prec = 64; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec = (prec < 1024) ? prec + 64 : prec + 256)
    {
        flint_printf("%4wd     ", prec);

        nfloat_ctx_init(ctx, prec, 0);

        for (ni = 0; (n = ns[ni]) != 0; ni++)
        {
            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            randmat(A, state, ctx);
            randmat(B, state, ctx);

            TIMEIT_START
            GR_MUST_SUCCEED(gr_mat_mul_classical(C, A, B, ctx));
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul(C, A, B, ctx));
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf("%.3f ", t1 / t2);
            fflush(stdout);

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);

            if (t1 > 3.0)
                break;
        }

        flint_printf("\n");
    }

    flint_rand_clear(state);
}

int main()
{
    //int tab[TABN];

    //tune_classical_vs_fixed(tab);
    //tune_fixed_vs_block(tab);
    //prof_classical_vs_fixed();
    //prof_fixed_vs_block();
    //tune_classical_vs_block(tab);
    prof_mul();
}
