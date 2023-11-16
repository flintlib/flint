/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

#ifdef __GNUC__
# define frexp __builtin_frexp
# define ldexp __builtin_ldexp
#else
# include <math.h>
#endif

/* Block size for better cache locality. */
#define BLOCK_SIZE 32

/* Don't convert to doubles when smaller than this block size. */
#define MIN_D_BLOCK_SIZE 5

/* With doubles, we can have an exponent range of about 1024, minus some
   slack for accumulated sums. */
#define DOUBLE_MAX_OFFSET 900

static inline double dot8(const double * A, const double * B)
{
    return ((A[0] * B[0] + A[1] * B[1]) + (A[2] * B[2] + A[3] * B[3])) +
           ((A[4] * B[4] + A[5] * B[5]) + (A[6] * B[6] + A[7] * B[7]));
}

/* Upper bound of matrix product, assuming nonnegative entries and
   no overflow/underflow. B is pre-transposed. Straightforward blocked
   implementation; could use BLAS, but this matrix product is rarely going
   to be the bottleneck. */

static void
_d_mat_addmul(double * C, const double * A, const double * B, slong ar, slong ac, slong bc)
{
    slong ii, jj, kk, i, j, k;
    double t, eps;

    eps = ldexp(1.0, -52);

    for (ii = 0; ii < ar; ii += BLOCK_SIZE)
    {
        for (jj = 0; jj < bc; jj += BLOCK_SIZE)
        {
            for (kk = 0; kk < ac; kk += BLOCK_SIZE)
            {
                for (i = ii; i < FLINT_MIN(ii + BLOCK_SIZE, ar); i++)
                {
                    for (j = jj; j < FLINT_MIN(jj + BLOCK_SIZE, bc); j++)
                    {
                        if (BLOCK_SIZE == 32 && kk + BLOCK_SIZE <= ac)
                        {
                            double t0, t1, t2, t3;

                            t0 = dot8(A + i * ac + kk + 0, B + j * ac + kk + 0);
                            t1 = dot8(A + i * ac + kk + 8, B + j * ac + kk + 8);
                            t2 = dot8(A + i * ac + kk + 16, B + j * ac + kk + 16);
                            t3 = dot8(A + i * ac + kk + 24, B + j * ac + kk + 24);

                            t = (t0 + t1) + (t2 + t3);
                        }
                        else
                        {
                            t = 0.0;

                            for (k = kk; k < FLINT_MIN(kk + BLOCK_SIZE, ac); k++)
                                t += A[i * ac + k] * B[j * ac + k];
                        }

                        C[i * bc + j] += t;
                    }
                }
            }
        }
    }

    /* Compensate for possible rounding errors */
    for (i = 0; i < ar; i++)
        for (j = 0; j < bc; j++)
            C[i * bc + j] *= (1.0 + 2.01 * (ac + 1) * eps);
}

/* We use WORD_MIN to represent zero here. */
static inline slong _mag_get_exp(const mag_t x)
{
    if (mag_is_special(x))
        return WORD_MIN;
    else
        return MAG_EXP(x);
}

static double
mag_get_d_fixed_si(const mag_t x, slong e)
{
    return ldexp(MAG_MAN(x), MAG_EXP(x) - e - MAG_BITS);
}

void
_arb_mat_addmul_rad_mag_fast(arb_mat_t C, mag_srcptr A, mag_srcptr B,
    slong ar, slong ac, slong bc)
{
    slong i, j, k, M, N, P, top, n, block_start, block_end;
    slong *A_min, *A_max, *B_min, *B_max, max_offset;
    double *CC, *AA, *BB;

    M = ar;
    N = ac;
    P = bc;

    /* todo: could use TMP_ALLOC */
    A_min = flint_malloc(sizeof(slong) * M);
    A_max = flint_malloc(sizeof(slong) * M);
    B_min = flint_malloc(sizeof(slong) * P);
    B_max = flint_malloc(sizeof(slong) * P);

    AA = flint_malloc(ar * ac * sizeof(double));
    BB = flint_malloc(ac * bc * sizeof(double));
    CC = flint_malloc(ar * bc * sizeof(double));

    max_offset = DOUBLE_MAX_OFFSET;

    block_start = 0;
    while (block_start < N)
    {
        block_end = block_start + 1;  /* index is exclusive block_end */

        /* begin with this column of A and row of B */
        for (i = 0; i < M; i++)
            A_max[i] = A_min[i] = _mag_get_exp(A + i * N + block_start);
        for (i = 0; i < P; i++)
            B_max[i] = B_min[i] = _mag_get_exp(B + i * N + block_start);

        while (block_end < N)
        {
            /* check if we can extend with column [block_end] of A */
            for (i = 0; i < M; i++)
            {
                top = _mag_get_exp(A + i * N + block_end);
                /* zeros are irrelevant */
                if (top == WORD_MIN || A_max[i] == WORD_MIN)
                    continue;
                /* jump will be too big */
                if (top > A_min[i] + max_offset || top < A_max[i] - max_offset)
                    goto mblocks_built;
            }

            /* check if we can extend with row [block_end] of B */
            for (i = 0; i < P; i++)
            {
                top = _mag_get_exp(B + i * N + block_end);
                if (top == WORD_MIN || B_max[i] == WORD_MIN)
                    continue;
                if (top > B_min[i] + max_offset || top < B_max[i] - max_offset)
                    goto mblocks_built;
            }

            /* second pass to update the extreme values */
            for (i = 0; i < M; i++)
            {
                top = _mag_get_exp(A + i * N + block_end);
                if (A_max[i] == WORD_MIN)
                {
                    A_max[i] = top;
                    A_min[i] = top;
                }
                else if (top != WORD_MIN)
                {
                    if (top < A_min[i]) A_min[i] = top;
                    if (top > A_max[i]) A_max[i] = top;
                }
            }

            for (i = 0; i < P; i++)
            {
                top = _mag_get_exp(B + i * N + block_end);
                if (B_max[i] == WORD_MIN)
                {
                    B_max[i] = top;
                    B_min[i] = top;
                }
                else if (top != WORD_MIN)
                {
                    if (top < B_min[i]) B_min[i] = top;
                    if (top > B_max[i]) B_max[i] = top;
                }
            }

            block_end++;
        }

    mblocks_built:

        n = block_end - block_start;

        if (n <= MIN_D_BLOCK_SIZE)
        {
            /* increment so we don't just do steps of 1 in degenerate cases */
            block_end = FLINT_MIN(block_start + MIN_D_BLOCK_SIZE, N);
            n = block_end - block_start;

            for (i = 0; i < ar; i++)
            {
                for (j = 0; j < bc; j++)
                {
                    for (k = 0; k < n; k++)
                    {
                        mag_fast_addmul(arb_radref(arb_mat_entry(C, i, j)),
                            A + i * ac + block_start + k,
                            B + j * ac + block_start + k);
                    }
                }
            }
        }
        else
        {
            for (i = 0; i < ar; i++)
            {
                if (A_min[i] == WORD_MIN)  /* only zeros in this row */
                    continue;

                A_min[i] = (A_min[i] + A_max[i]) / 2;

                for (j = 0; j < n; j++)
                    AA[i * n + j] = mag_get_d_fixed_si(A + i * ac + block_start + j, A_min[i]);
            }

            /* Note: B and BB are both transposed in memory */
            for (i = 0; i < bc; i++)
            {
                if (B_min[i] == WORD_MIN)  /* only zeros in this column */
                    continue;

                B_min[i] = (B_min[i] + B_max[i]) / 2;

                for (j = 0; j < n; j++)
                    BB[i * n + j] = mag_get_d_fixed_si(B + i * ac + block_start + j, B_min[i]);
            }

            for (i = 0; i < ar * bc; i++)
                CC[i] = 0.0;

            _d_mat_addmul(CC, AA, BB, ar, n, bc);

            for (i = 0; i < ar; i++)
            {
                if (A_min[i] == WORD_MIN)
                    continue;

                for (j = 0; j < bc; j++)
                {
                    if (B_min[j] == WORD_MIN)
                        continue;

                    if (CC[i * bc + j] != 0.0)
                    {
                        mag_t t;
                        MAG_SET_D_2EXP(MAG_MAN(t), MAG_EXP(t), CC[i * bc + j], A_min[i] + B_min[j]);
                        mag_add(arb_radref(arb_mat_entry(C, i, j)),
                                arb_radref(arb_mat_entry(C, i, j)), t);
                    }
                }
            }
        }

        block_start = block_end;
    }

    flint_free(A_max);
    flint_free(A_min);
    flint_free(B_max);
    flint_free(B_min);

    flint_free(AA);
    flint_free(BB);
    flint_free(CC);
}

