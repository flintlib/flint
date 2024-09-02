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
#define WAKSMAN_MIN_PREC 320

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

void tune_fixed_vs_waksman(int * cutoffs)
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

    for (prec = WAKSMAN_MIN_PREC; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += 64)
    {
        flint_printf("prec = %wd\n", prec);

        nfloat_ctx_init(ctx, prec, 0);

        for (n = 2; n <= 128; n++)
        {
            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            randmat(A, state, ctx);
            randmat(B, state, ctx);

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul_fixed_classical(C, A, B, ctx));
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul_waksman(C, A, B, ctx));
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf("%wd  %wd   %e   %e   %.3f\n", prec, n, t1, t2, t1 / t2);

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);

            if (t2 < t1 * 0.99)
            {
                cutoffs[prec / 64] = n;

                flint_printf("short tab_fixed_classical_vs_waksman[] = {\n");
                for (i = 0; i <= prec / 64; i++)
                    flint_printf("    %d, /* prec = %wd */\n", cutoffs[i], i * 64);
                flint_printf("}\n");

                break;
            }
        }
    }

    flint_rand_clear(state);
}

void tune_strassen(int * cutoffs)
{
    gr_ctx_t ctx;
    gr_mat_t A, B, C;
    slong i, n;
    slong prec;
    double FLINT_SET_BUT_UNUSED(__), t1, t2;
    int prev_ok = 0;

    flint_rand_t state;
    flint_rand_init(state);

    for (i = 0; i < TABN; i++)
        cutoffs[i] = -1;

    for (prec = 64; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += 64)
    {
        flint_printf("prec = %wd\n", prec);

        prev_ok = 0;

        nfloat_ctx_init(ctx, prec, 0);

        for (n = 2; n <= 128; n++)
        {
            gr_mat_init(A, n, n, ctx);
            gr_mat_init(B, n, n, ctx);
            gr_mat_init(C, n, n, ctx);

            randmat(A, state, ctx);
            randmat(B, state, ctx);

            TIMEIT_START
            if (prec < WAKSMAN_MIN_PREC)
                GR_MUST_SUCCEED(nfloat_mat_mul_fixed_classical(C, A, B, ctx));
            else
                GR_MUST_SUCCEED(nfloat_mat_mul_waksman(C, A, B, ctx));
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            GR_MUST_SUCCEED(nfloat_mat_mul_strassen(C, A, B, n, ctx));
            TIMEIT_STOP_VALUES(__, t2)

            flint_printf("%wd  %wd   %e   %e   %.3f\n", prec, n, t1, t2, t1 / t2);

            gr_mat_clear(A, ctx);
            gr_mat_clear(B, ctx);
            gr_mat_clear(C, ctx);

            if (t2 < t1 * 0.99)
            {
                if (prev_ok)
                {
                    cutoffs[prec / 64] = n;

                    flint_printf("short tab_strassen[] = {\n");
                    for (i = 0; i <= prec / 64; i++)
                        flint_printf("    %d, /* prec = %wd */\n", cutoffs[i], i * 64);
                    flint_printf("}\n");

                    break;
                }
                else
                {
                    prev_ok = 1;
                }
            }
            else
            {
                prev_ok = 0;
            }
        }
    }

    flint_rand_clear(state);
}


short tab_fixed_classical_vs_waksman[] = {
    -1, /* prec = 0 */
    -1, /* prec = 64 */
    -1, /* prec = 128 */
    -1, /* prec = 192 */
    -1, /* prec = 256 */
    16, /* prec = 320 */
    10, /* prec = 384 */
    7, /* prec = 448 */
    7, /* prec = 512 */
    6, /* prec = 576 */
    5, /* prec = 640 */
    4, /* prec = 704 */
    4, /* prec = 768 */
    4, /* prec = 832 */
    4, /* prec = 896 */
    4, /* prec = 960 */
    4, /* prec = 1024 */
    4, /* prec = 1088 */
    4, /* prec = 1152 */
    4, /* prec = 1216 */
    4, /* prec = 1280 */
    3, /* prec = 1344 */
    4, /* prec = 1408 */
    3, /* prec = 1472 */
    3, /* prec = 1536 */
    3, /* prec = 1600 */
    3, /* prec = 1664 */
    3, /* prec = 1728 */
    3, /* prec = 1792 */
    3, /* prec = 1856 */
    3, /* prec = 1920 */
    3, /* prec = 1984 */
    3, /* prec = 2048 */
    3, /* prec = 2112 */
    3, /* prec = 2176 */
    3, /* prec = 2240 */
    3, /* prec = 2304 */
    3, /* prec = 2368 */
    3, /* prec = 2432 */
    3, /* prec = 2496 */
    3, /* prec = 2560 */
    3, /* prec = 2624 */
    3, /* prec = 2688 */
    3, /* prec = 2752 */
    3, /* prec = 2816 */
    3, /* prec = 2880 */
    3, /* prec = 2944 */
    3, /* prec = 3008 */
    2, /* prec = 3072 */
    3, /* prec = 3136 */
    3, /* prec = 3200 */
    2, /* prec = 3264 */
    2, /* prec = 3328 */
    2, /* prec = 3392 */
    2, /* prec = 3456 */
    2, /* prec = 3520 */
    3, /* prec = 3584 */
    2, /* prec = 3648 */
    2, /* prec = 3712 */
    2, /* prec = 3776 */
    2, /* prec = 3840 */
    2, /* prec = 3904 */
    2, /* prec = 3968 */
    2, /* prec = 4032 */
    2, /* prec = 4096 */
    2, /* prec = 4160 */
    2, /* prec = 4224 */
};



int main()
{
    int tab_fixed_classical_vs_waksman[TABN];
    int tab_strassen[TABN];

    tune_strassen(tab_strassen);


    tune_fixed_vs_waksman(tab_fixed_classical_vs_waksman);

}
