/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Times nmod_mat_transpose, out of place (m x n -> n x m) and in place
    (square), against the entrywise loops it replaces, in nanoseconds per
    entry. With arguments "m n", times that shape only.
*/

#include <stdlib.h>
#include "profiler.h"
#include "nmod_mat.h"

typedef struct
{
    slong m, n;
    int algo;   /* 0: nmod_mat_transpose, 1: entrywise */
    int inplace;
} info_t;

/* the loops of nmod_mat_transpose in FLINT 3.4, not inlined so as to
   compare calls with calls */
FLINT_STATIC_NOINLINE void
transpose_entrywise(nmod_mat_t B, const nmod_mat_t A)
{
    slong i, j;

    if (A == B)
    {
        for (i = 0; i < A->r - 1; i++)
            for (j = i + 1; j < A->c; j++)
                FLINT_SWAP(ulong, nmod_mat_entry(A, i, j), nmod_mat_entry(A, j, i));
    }
    else
    {
        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
                nmod_mat_entry(B, i, j) = nmod_mat_entry(A, j, i);
    }
}

static void
sample(void * arg, ulong count)
{
    info_t * info = (info_t *) arg;
    nmod_mat_t A, B;
    ulong i;
    FLINT_TEST_INIT(state);

    nmod_mat_init(A, info->m, info->n, UWORD(1) << 62);
    nmod_mat_init(B, info->n, info->m, UWORD(1) << 62);
    nmod_mat_randfull(A, state);
    nmod_mat_randfull(B, state);

    prof_start();
    for (i = 0; i < count; i++)
    {
        if (info->inplace)
        {
            if (info->algo == 0)
                nmod_mat_transpose(A, A);
            else
                transpose_entrywise(A, A);
        }
        else
        {
            if (info->algo == 0)
                nmod_mat_transpose(B, A);
            else
                transpose_entrywise(B, A);
        }
    }
    prof_stop();

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    FLINT_TEST_CLEAR(state);
}

/* ns per entry */
static double
time_ns(slong m, slong n, int algo, int inplace)
{
    info_t info;
    double min, max;

    info.m = m;
    info.n = n;
    info.algo = algo;
    info.inplace = inplace;
    prof_repeat(&min, &max, sample, &info);

    return 1000 * min / (m * n);
}

static void
profile_shape(slong m, slong n)
{
    double t, t0;

    flint_printf("%5wd x %-5wd", m, n);

    t = time_ns(m, n, 0, 0);
    t0 = time_ns(m, n, 1, 0);
    flint_printf("  %8.3f %8.3f  %5.2fx", t0, t, t0 / t);

    if (m == n)
    {
        t = time_ns(m, n, 0, 1);
        t0 = time_ns(m, n, 1, 1);
        flint_printf("  %8.3f %8.3f  %5.2fx", t0, t, t0 / t);
    }

    flint_printf("\n");
}

int
main(int argc, char ** argv)
{
    slong dims[] = {4, 8, 16, 30, 32, 64, 100, 128, 250, 256, 500, 512,
                    1000, 1024, 2000, 2048, 4000, 4096};
    slong rect[][2] = {{1, 10000}, {10000, 1}, {4, 10000}, {10000, 4},
                       {10, 10000}, {10000, 10}, {100, 10000}, {10000, 100},
                       {1000, 4000}, {4000, 1000}};
    slong i;

    flint_printf("nmod_mat_transpose, ns per entry: entrywise loop, nmod_mat_transpose, ratio\n");
    flint_printf("%-13s  %-26s  %s\n", "m x n", "  out of place", "  in place (square)");

    if (argc == 3)
    {
        profile_shape(atol(argv[1]), atol(argv[2]));
        return 0;
    }

    for (i = 0; i < (slong) (sizeof(dims) / sizeof(dims[0])); i++)
        profile_shape(dims[i], dims[i]);

    for (i = 0; i < (slong) (sizeof(rect) / sizeof(rect[0])); i++)
        profile_shape(rect[i][0], rect[i][1]);

    return 0;
}
