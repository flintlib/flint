/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_factor.h"

typedef struct
{
    slong rows;
    slong cols;
    slong bits;
    slong rank;
} snf_params_t;

/*
    Section 0: SNF dispatcher timing.
    Generate matrix with given rank, apply random ops, time fmpz_mat_snf.
*/
void sample_snf(void * arg, ulong count)
{
    snf_params_t * params = (snf_params_t *) arg;
    slong m = params->rows, n = params->cols;
    slong bits = params->bits, rank = params->rank;
    ulong i;
    fmpz_mat_t A, S;
    FLINT_TEST_INIT(state);

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(S, m, n);

    fmpz_mat_randrank(A, state, rank, bits);
    fmpz_mat_randops(A, state, 2 * m * n);

    prof_start();
    for (i = 0; i < count; i++)
        fmpz_mat_snf(S, A);
    prof_stop();

    fmpz_mat_clear(S);
    fmpz_mat_clear(A);
    FLINT_TEST_CLEAR(state);
}

/*
    Section 1: HNF with and without transform.
*/
void sample_hnf(void * arg, ulong count)
{
    snf_params_t * params = (snf_params_t *) arg;
    slong m = params->rows, n = params->cols;
    slong bits = params->bits, rank = params->rank;
    ulong i;
    fmpz_mat_t A, H;
    FLINT_TEST_INIT(state);

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(H, m, n);

    fmpz_mat_randrank(A, state, rank, bits);
    fmpz_mat_randops(A, state, 2 * m * n);

    prof_start();
    for (i = 0; i < count; i++)
        fmpz_mat_hnf(H, A);
    prof_stop();

    fmpz_mat_clear(H);
    fmpz_mat_clear(A);
    FLINT_TEST_CLEAR(state);
}

void sample_hnf_transform(void * arg, ulong count)
{
    snf_params_t * params = (snf_params_t *) arg;
    slong m = params->rows, n = params->cols;
    slong bits = params->bits, rank = params->rank;
    ulong i;
    fmpz_mat_t A, H, U;
    FLINT_TEST_INIT(state);

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(H, m, n);
    fmpz_mat_init(U, m, m);

    fmpz_mat_randrank(A, state, rank, bits);
    fmpz_mat_randops(A, state, 2 * m * n);

    prof_start();
    for (i = 0; i < count; i++)
        fmpz_mat_hnf_transform(H, U, A);
    prof_stop();

    fmpz_mat_clear(U);
    fmpz_mat_clear(H);
    fmpz_mat_clear(A);
    FLINT_TEST_CLEAR(state);
}

/*
    Section 2: det + factor feasibility.
*/
void sample_det(void * arg, ulong count)
{
    snf_params_t * params = (snf_params_t *) arg;
    slong n = params->rows, bits = params->bits;
    ulong i;
    fmpz_mat_t A;
    fmpz_t d;
    FLINT_TEST_INIT(state);

    fmpz_mat_init(A, n, n);
    fmpz_init(d);

    fmpz_mat_randrank(A, state, n, bits);
    fmpz_mat_randops(A, state, 2 * n * n);

    prof_start();
    for (i = 0; i < count; i++)
        fmpz_mat_det(d, A);
    prof_stop();

    fmpz_clear(d);
    fmpz_mat_clear(A);
    FLINT_TEST_CLEAR(state);
}

void sample_det_factor(void * arg, ulong count)
{
    snf_params_t * params = (snf_params_t *) arg;
    slong n = params->rows, bits = params->bits;
    ulong i;
    fmpz_mat_t A;
    fmpz_t d;
    fmpz_factor_t fac;
    FLINT_TEST_INIT(state);

    fmpz_mat_init(A, n, n);
    fmpz_init(d);
    fmpz_factor_init(fac);

    fmpz_mat_randrank(A, state, n, bits);
    fmpz_mat_randops(A, state, 2 * n * n);
    fmpz_mat_det(d, A);
    fmpz_abs(d, d);

    prof_start();
    for (i = 0; i < count; i++)
        fmpz_factor(fac, d);
    prof_stop();

    fmpz_factor_clear(fac);
    fmpz_clear(d);
    fmpz_mat_clear(A);
    FLINT_TEST_CLEAR(state);
}

/*
    Section 3: Iterative Hermite building blocks.
    Time one iteration: hnf_transform + transpose + mul.
*/
void sample_hermite_iter(void * arg, ulong count)
{
    snf_params_t * params = (snf_params_t *) arg;
    slong n = params->rows, bits = params->bits;
    ulong i;
    fmpz_mat_t X, M, U;
    FLINT_TEST_INIT(state);

    fmpz_mat_init(X, n, n);
    fmpz_mat_init(M, n, n);
    fmpz_mat_init(U, n, n);

    fmpz_mat_randrank(X, state, n, bits);
    fmpz_mat_randops(X, state, 2 * n * n);
    fmpz_mat_one(U);

    prof_start();
    for (i = 0; i < count; i++)
    {
        fmpz_mat_hnf_transform(X, M, X);
        fmpz_mat_mul(U, M, U);
        fmpz_mat_transpose(X, X);
    }
    prof_stop();

    fmpz_mat_clear(U);
    fmpz_mat_clear(M);
    fmpz_mat_clear(X);
    FLINT_TEST_CLEAR(state);
}

int main(void)
{
    double min, max;
    snf_params_t params;
    slong dim, bits;
    slong dims[] = {5, 10, 20, 30, 50, 75, 100};
    slong bitvals[] = {10, 50, 200};
    slong ndims = sizeof(dims) / sizeof(dims[0]);
    slong nbits = sizeof(bitvals) / sizeof(bitvals[0]);
    slong i, j;

    /* ================================================================ */
    flint_printf("=== Section 0: SNF dispatcher timing ===\n\n");

    /* Square nonsingular */
    flint_printf("--- Square nonsingular ---\n");
    flint_printf("%5s", "dim");
    for (j = 0; j < nbits; j++)
        flint_printf("  %8wd-bit", bitvals[j]);
    flint_printf("\n");

    for (i = 0; i < ndims; i++)
    {
        dim = dims[i];
        flint_printf("%5wd", dim);
        for (j = 0; j < nbits; j++)
        {
            params.rows = dim;
            params.cols = dim;
            params.bits = bitvals[j];
            params.rank = dim;
            prof_repeat(&min, &max, sample_snf, &params);
            flint_printf("  %12.1f", min);
        }
        flint_printf("  (us)\n");
    }
    flint_printf("\n");

    /* Square singular (rank n/2) */
    flint_printf("--- Square singular (rank n/2) ---\n");
    flint_printf("%5s", "dim");
    for (j = 0; j < nbits; j++)
        flint_printf("  %8wd-bit", bitvals[j]);
    flint_printf("\n");

    for (i = 0; i < ndims; i++)
    {
        dim = dims[i];
        flint_printf("%5wd", dim);
        for (j = 0; j < nbits; j++)
        {
            params.rows = dim;
            params.cols = dim;
            params.bits = bitvals[j];
            params.rank = dim / 2;
            prof_repeat(&min, &max, sample_snf, &params);
            flint_printf("  %12.1f", min);
        }
        flint_printf("  (us)\n");
    }
    flint_printf("\n");

    /* Non-square m x 2m */
    flint_printf("--- Non-square m x 2m ---\n");
    flint_printf("%5s", "m");
    for (j = 0; j < nbits; j++)
        flint_printf("  %8wd-bit", bitvals[j]);
    flint_printf("\n");

    for (i = 0; i < ndims - 1; i++)
    {
        dim = dims[i];
        flint_printf("%5wd", dim);
        for (j = 0; j < nbits; j++)
        {
            params.rows = dim;
            params.cols = 2 * dim;
            params.bits = bitvals[j];
            params.rank = dim;
            prof_repeat(&min, &max, sample_snf, &params);
            flint_printf("  %12.1f", min);
        }
        flint_printf("  (us)\n");
    }
    flint_printf("\n");

    /* ================================================================ */
    flint_printf("=== Section 1: HNF with/without transform ===\n\n");
    flint_printf("%5s", "dim");
    for (j = 0; j < nbits; j++)
        flint_printf("  %8wd-bit  %8wd-bit", bitvals[j], bitvals[j]);
    flint_printf("\n");
    flint_printf("%5s", "");
    for (j = 0; j < nbits; j++)
        flint_printf("  %10s  %10s", "hnf", "hnf+U");
    flint_printf("\n");

    for (i = 0; i < ndims; i++)
    {
        double min_hnf, min_hnf_t;
        dim = dims[i];
        flint_printf("%5wd", dim);
        for (j = 0; j < nbits; j++)
        {
            params.rows = dim;
            params.cols = dim;
            params.bits = bitvals[j];
            params.rank = dim;
            prof_repeat(&min_hnf, &max, sample_hnf, &params);
            prof_repeat(&min_hnf_t, &max, sample_hnf_transform, &params);
            flint_printf("  %10.1f  %10.1f", min_hnf, min_hnf_t);
        }
        flint_printf("  (us)\n");
    }
    flint_printf("\n");

    /* ================================================================ */
    flint_printf("=== Section 2: det + factor feasibility ===\n\n");
    flint_printf("%5s", "dim");
    for (j = 0; j < nbits; j++)
        flint_printf("  %8wd-bit  %8wd-bit", bitvals[j], bitvals[j]);
    flint_printf("\n");
    flint_printf("%5s", "");
    for (j = 0; j < nbits; j++)
        flint_printf("  %10s  %10s", "det", "factor");
    flint_printf("\n");

    for (i = 0; i < ndims; i++)
    {
        double min_det, min_fac;
        dim = dims[i];
        flint_printf("%5wd", dim);
        for (j = 0; j < nbits; j++)
        {
            params.rows = dim;
            params.cols = dim;
            params.bits = bitvals[j];
            params.rank = dim;
            prof_repeat(&min_det, &max, sample_det, &params);
            prof_repeat(&min_fac, &max, sample_det_factor, &params);
            flint_printf("  %10.1f  %10.1f", min_det, min_fac);
        }
        flint_printf("  (us)\n");
    }
    flint_printf("\n");

    /* ================================================================ */
    flint_printf("=== Section 3: Iterative Hermite (one iteration) ===\n\n");
    flint_printf("%5s", "dim");
    for (j = 0; j < nbits; j++)
        flint_printf("  %8wd-bit", bitvals[j]);
    flint_printf("\n");

    for (i = 0; i < ndims; i++)
    {
        dim = dims[i];
        flint_printf("%5wd", dim);
        for (j = 0; j < nbits; j++)
        {
            params.rows = dim;
            params.cols = dim;
            params.bits = bitvals[j];
            params.rank = dim;
            prof_repeat(&min, &max, sample_hermite_iter, &params);
            flint_printf("  %12.1f", min);
        }
        flint_printf("  (us)\n");
    }
    flint_printf("\n");

    return 0;
}
