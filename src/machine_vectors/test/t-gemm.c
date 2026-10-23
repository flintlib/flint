/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "machine_vectors.h"
#include "ulong_extras.h"
#include "thread_support.h"

/*
    Entries are small integers, so every product is exact in both
    precisions and the reference comparison is an equality test. This
    covers all the internal paths: packed, small (all dims <= 96), thin-m
    (m <= MR), thin-n (n <= NR), the k = 0 case, and every combination
    with leading dimensions larger than the matrices.

    flint_sgemm_fallback/flint_dgemm_fallback are tested rather than
    flint_sgemm/flint_dgemm so that our own kernels are exercised even
    in a BLAS build; the dispatch itself is tested separately below.
*/

TEST_FUNCTION_START(machine_vectors_gemm, state)
{
    slong iter;

    /* the dispatcher must agree with both implementations */
    {
        slong m = 40, k = 33, n = 51, i, j, l;
        double * A = flint_malloc(m * k * sizeof(double));
        double * B = flint_malloc(k * n * sizeof(double));
        double * C = flint_malloc(m * n * sizeof(double));
        int saved = flint_gemm_use_blas;
        /* FLINT_USES_BLAS is undefined rather than 0 without BLAS */
#if FLINT_USES_BLAS
        int max_use_blas = 1;
#else
        int max_use_blas = 0;
#endif

        for (i = 0; i < m * k; i++)
            A[i] = (double) (slong) (n_randint(state, 21) - 10);
        for (i = 0; i < k * n; i++)
            B[i] = (double) (slong) (n_randint(state, 21) - 10);

        for (flint_gemm_use_blas = 0; flint_gemm_use_blas <= max_use_blas;
             flint_gemm_use_blas++)
        {
            flint_dgemm(m, k, n, A, k, B, n, C, n);

            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                {
                    double s = 0.0;

                    for (l = 0; l < k; l++)
                        s += A[i * k + l] * B[l * n + j];

                    if (C[i * n + j] != s)
                        TEST_FUNCTION_FAIL("dispatch (use_blas = %d) "
                            "wrong at %wd, %wd: got %g, want %g\n",
                            flint_gemm_use_blas, i, j, C[i * n + j], s);
                }
        }

        flint_gemm_use_blas = saved;
        flint_free(A); flint_free(B); flint_free(C);
    }

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong m, k, n, lda, ldb, ldc, i, j, l;
        float *Af, *Bf, *Cf;
        double *Ad, *Bd, *Cd;
        int nthreads;

        if (n_randint(state, 4) == 0)
        {
            /* occasionally exercise sizes that reach the packed path */
            m = 1 + n_randint(state, 200);
            k = 1 + n_randint(state, 200);
            n = 1 + n_randint(state, 200);
        }
        else
        {
            m = n_randint(state, 40);
            k = n_randint(state, 40);
            n = n_randint(state, 40);
        }

        lda = k + n_randint(state, 3);
        ldb = n + n_randint(state, 3);
        ldc = n + n_randint(state, 3);

        nthreads = 1 + n_randint(state, 4);
        flint_set_num_threads(nthreads);

        Af = flint_malloc(FLINT_MAX(1, m * lda) * sizeof(float));
        Bf = flint_malloc(FLINT_MAX(1, k * ldb) * sizeof(float));
        Cf = flint_malloc(FLINT_MAX(1, m * ldc) * sizeof(float));
        Ad = flint_malloc(FLINT_MAX(1, m * lda) * sizeof(double));
        Bd = flint_malloc(FLINT_MAX(1, k * ldb) * sizeof(double));
        Cd = flint_malloc(FLINT_MAX(1, m * ldc) * sizeof(double));

        for (i = 0; i < m * lda; i++)
        {
            Ad[i] = (double) (slong) (n_randint(state, 21) - 10);
            Af[i] = (float) Ad[i];
        }

        for (i = 0; i < k * ldb; i++)
        {
            Bd[i] = (double) (slong) (n_randint(state, 21) - 10);
            Bf[i] = (float) Bd[i];
        }

        for (i = 0; i < m * ldc; i++)
        {
            Cd[i] = 1234.0;
            Cf[i] = 1234.0f;
        }

        /* always exercise our own kernels, whatever the BLAS default is */
        flint_sgemm_fallback(m, k, n, Af, lda, Bf, ldb, Cf, ldc);
        flint_dgemm_fallback(m, k, n, Ad, lda, Bd, ldb, Cd, ldc);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                double s = 0.0;

                for (l = 0; l < k; l++)
                    s += Ad[i * lda + l] * Bd[l * ldb + j];

                if (Cd[i * ldc + j] != s || Cf[i * ldc + j] != (float) s)
                    TEST_FUNCTION_FAIL(
                        "m = %wd, k = %wd, n = %wd\n"
                        "lda = %wd, ldb = %wd, ldc = %wd\n"
                        "threads = %d\n"
                        "i = %wd, j = %wd\n"
                        "sgemm gave %g, dgemm gave %g, want %g\n",
                        m, k, n, lda, ldb, ldc, nthreads, i, j,
                        (double) Cf[i * ldc + j], Cd[i * ldc + j], s);
            }

            /* padding columns of C must be untouched */
            for (j = n; j < ldc; j++)
                if (Cd[i * ldc + j] != 1234.0 || Cf[i * ldc + j] != 1234.0f)
                    TEST_FUNCTION_FAIL(
                        "wrote past column n\n"
                        "m = %wd, k = %wd, n = %wd, ldc = %wd, i = %wd, "
                        "j = %wd\n", m, k, n, ldc, i, j);
        }

        flint_free(Af); flint_free(Bf); flint_free(Cf);
        flint_free(Ad); flint_free(Bd); flint_free(Cd);
    }

    /*
        Cleanup and reuse cycles: the packing scratch registers a cleanup
        function with FLINT, and must re-register correctly after
        flint_cleanup has emptied the registry.
    */
    {
        slong d = 200, i;
        double * A = flint_malloc(d * d * sizeof(double));
        double * B = flint_malloc(d * d * sizeof(double));
        double * C = flint_malloc(d * d * sizeof(double));

        for (i = 0; i < d * d; i++)
        {
            A[i] = (double) (i % 5);
            B[i] = (double) (i % 3);
        }

        for (i = 0; i < 3; i++)
        {
            flint_set_num_threads(4);
            flint_dgemm(d, d, d, A, d, B, d, C, d);
            flint_cleanup();
        }

        flint_free(A); flint_free(B); flint_free(C);
    }

    flint_set_num_threads(1);

    TEST_FUNCTION_END(state);
}
