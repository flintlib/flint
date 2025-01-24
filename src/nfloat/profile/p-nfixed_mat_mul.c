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

void
nfixed_rand(nn_ptr a, flint_rand_t state, slong nlimbs)
{
    a[0] = n_randint(state, 2);
    flint_mpn_rrandom(a + 1, state, nlimbs);
    a[nlimbs] >>= 10;
}

void
nfixed_randmat(nn_ptr a, slong m, slong n, flint_rand_t state, slong nlimbs)
{
    slong i;
    for (i = 0; i < m * n; i++)
        nfixed_rand(a + i * (nlimbs + 1), state, nlimbs);
}

void tune_fixed_vs_waksman(int * cutoffs)
{
    nn_ptr A, B, C;
    slong i, n, nlimbs, nn;
    double FLINT_SET_BUT_UNUSED(__), t1, t2;

    flint_rand_t state;
    flint_rand_init(state);

    for (i = 0; i < TABN; i++)
        cutoffs[i] = -1;

    for (nlimbs = 2; nlimbs <= NFLOAT_MAX_LIMBS; nlimbs++)
    {
        flint_printf("nlimbs = %wd\n", nlimbs);

        for (nn = 1; nn <= 64; nn++)
        {
            n = (nn == 1) ? 128 : nn;

            A = flint_malloc((n * n) * (nlimbs + 1) * sizeof(ulong));
            B = flint_malloc((n * n) * (nlimbs + 1) * sizeof(ulong));
            C = flint_malloc((n * n) * (nlimbs + 1) * sizeof(ulong));

            nfixed_randmat(A, n, n, state, nlimbs);
            nfixed_randmat(B, n, n, state, nlimbs);

            TIMEIT_START
            _nfixed_mat_mul_classical(C, A, B, n, n, n, nlimbs);
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            _nfixed_mat_mul_waksman(C, A, B, n, n, n, nlimbs);
            TIMEIT_STOP_VALUES(__, t2)

            flint_free(A);
            flint_free(B);
            flint_free(C);

            flint_printf("%wd  %wd   %e   %e   %.3f\n", nlimbs * FLINT_BITS, n, t1, t2, t1 / t2);

            if (nn == 1)
            {
                if (t2 < t1)
                    continue;
                else
                    break;
            }

            if (t2 < t1)
            {
                cutoffs[nlimbs] = n;

                flint_printf("short tab_fixed_classical_vs_waksman[] = {\n");
                for (i = 0; i <= nlimbs; i++)
                    flint_printf("    %d, /* nlimbs = %wd */\n", cutoffs[i], i);
                flint_printf("}\n");

                break;
            }
        }
    }
}

void tune_strassen(int * cutoffs)
{
    nn_ptr A, B, C;
    slong i, n, nlimbs;
    int parity;
    int last_ok;
    double FLINT_SET_BUT_UNUSED(__), t1, t2;

    flint_rand_t state;
    flint_rand_init(state);

    for (parity = 0; parity < 2; parity++)
    {
        for (i = 0; i < TABN; i++)
            cutoffs[i] = -1;

        for (nlimbs = 2; nlimbs <= NFLOAT_MAX_LIMBS; nlimbs++)
        {
            flint_printf("nlimbs = %wd\n", nlimbs);

            last_ok = 0;

            for (n = parity ? 1 : 2; ; n += 2)
            {
                A = flint_malloc((n * n) * (nlimbs + 1) * sizeof(ulong));
                B = flint_malloc((n * n) * (nlimbs + 1) * sizeof(ulong));
                C = flint_malloc((n * n) * (nlimbs + 1) * sizeof(ulong));

                nfixed_randmat(A, n, n, state, nlimbs);
                nfixed_randmat(B, n, n, state, nlimbs);

                TIMEIT_START
                _nfixed_mat_mul_strassen(C, A, B, n, n, n, n + 1, nlimbs);
                TIMEIT_STOP_VALUES(__, t1)

                TIMEIT_START
                _nfixed_mat_mul_strassen(C, A, B, n, n, n, n, nlimbs);
                TIMEIT_STOP_VALUES(__, t2)

                flint_free(A);
                flint_free(B);
                flint_free(C);

                flint_printf("%wd  %wd   %e   %e   %.3f\n", nlimbs * FLINT_BITS, n, t1, t2, t1 / t2);

                if (t2 < t1)
                {
                    if (!last_ok)
                    {
                        last_ok = 1;
                        continue;
                    }

                    cutoffs[nlimbs] = n;

                    if (parity)
                        flint_printf("short tab_strassen_odd[] = {\n");
                    else
                        flint_printf("short tab_strassen_even[] = {\n");
                    for (i = 0; i <= nlimbs; i++)
                        flint_printf("    %d, /* nlimbs = %wd */\n", cutoffs[i], i);
                    flint_printf("}\n");

                    break;
                }
                else
                {
                    last_ok = 0;
                }
            }
        }
    }
}

void prof_strassen_1()
{
    nn_ptr A, B, C;
    slong n, nlimbs;
    int parity;
    double FLINT_SET_BUT_UNUSED(__), t1, t2;

    flint_rand_t state;
    flint_rand_init(state);

    for (nlimbs = 2; nlimbs <= NFLOAT_MAX_LIMBS; nlimbs++)
    {
        for (parity = 0; parity < 2; parity++)
        {
            flint_printf("nlimbs = %wd     ", nlimbs);

            if (nlimbs <= 3)
                n = parity ? 57 : 50;
            else
                n = parity ? 37 : 26;

            A = flint_malloc((n * n) * (nlimbs + 1) * sizeof(ulong));
            B = flint_malloc((n * n) * (nlimbs + 1) * sizeof(ulong));
            C = flint_malloc((n * n) * (nlimbs + 1) * sizeof(ulong));

            nfixed_randmat(A, n, n, state, nlimbs);
            nfixed_randmat(B, n, n, state, nlimbs);

            TIMEIT_START
            _nfixed_mat_mul_strassen(C, A, B, n, n, n, n + 1, nlimbs);
            TIMEIT_STOP_VALUES(__, t1)

            TIMEIT_START
            _nfixed_mat_mul_strassen(C, A, B, n, n, n, n, nlimbs);
            TIMEIT_STOP_VALUES(__, t2)

            flint_free(A);
            flint_free(B);
            flint_free(C);

            flint_printf("%wd  %e   %e   %.3f    ", n, t1, t2, t1 / t2);
        }

        flint_printf("\n");
    }
}

int main()
{
    //int tab_fixed_classical_vs_waksman[TABN];
    //int tab_strassen[TABN];

    //tune_fixed_vs_waksman(tab_fixed_classical_vs_waksman);
    //tune_strassen(tab_strassen);
    prof_strassen_1();
}
