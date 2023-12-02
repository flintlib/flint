/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "arb_mat.h"

int arb_mat_is_lagom(const arb_mat_t A)
{
    slong i, j, M, N;

    M = arb_mat_nrows(A);
    N = arb_mat_ncols(A);

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (!ARB_IS_LAGOM(arb_mat_entry(A, i, j)))
                return 0;
        }
    }

    return 1;
}

/* allow changing this from the test code */
FLINT_DLL slong arb_mat_mul_block_min_block_size = 0;

void
arb_mat_mid_addmul_block_fallback(arb_mat_t C,
    const arb_mat_t A, const arb_mat_t B,
    slong block_start,
    slong block_end,
    slong prec)
{
    slong M, P, n;
    slong i, j;
    arb_ptr tmpA, tmpB;

    M = arb_mat_nrows(A);
    P = arb_mat_ncols(B);

    n = block_end - block_start;

    tmpA = flint_malloc(sizeof(arb_struct) * (M * n + P * n));
    tmpB = tmpA + M * n;

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < n; j++)
        {
            *arb_midref(tmpA + i * n + j) = *arb_midref(arb_mat_entry(A, i, block_start + j));
            mag_init(arb_radref(tmpA + i * n + j));
        }
    }

    for (i = 0; i < P; i++)
    {
        for (j = 0; j < n; j++)
        {
            *arb_midref(tmpB + i * n + j) = *arb_midref(arb_mat_entry(B, block_start + j, i));
            mag_init(arb_radref(tmpB + i * n + j));
        }
    }

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < P; j++)
        {
            arb_dot(arb_mat_entry(C, i, j),
                (block_start == 0) ? NULL : arb_mat_entry(C, i, j), 0,
                tmpA + i * n, 1, tmpB + j * n, 1, n, prec);
        }
    }

    flint_free(tmpA);
}

void
arb_mat_mid_addmul_block_prescaled(arb_mat_t C,
    const arb_mat_t A, const arb_mat_t B,
    slong block_start,
    slong block_end,
    const slong * A_min,  /* A per-row bottom exponent */
    const slong * B_min,  /* B per-row bottom exponent */
    slong prec)
{
    slong M, P, n;
    slong i, j;
    slong M0, M1, P0, P1, Mstep, Pstep;
    int inexact;

    /* flint_printf("block mul from %wd to %wd\n", block_start, block_end); */

    M = arb_mat_nrows(A);
    P = arb_mat_ncols(B);

    n = block_end - block_start;

    /* Create sub-blocks to keep matrices nearly square. Necessary? */
#if 1
    Mstep = (M < 2 * n) ? M : n;
    Pstep = (P < 2 * n) ? P : n;
#else
    Mstep = M;
    Pstep = P;
#endif

    for (M0 = 0; M0 < M; M0 += Mstep)
    {
        for (P0 = 0; P0 < P; P0 += Pstep)
        {
            fmpz_mat_t AA, BB, CC;
            arb_t t;
            fmpz_t e;

            M1 = FLINT_MIN(M0 + Mstep, M);
            P1 = FLINT_MIN(P0 + Pstep, P);

            fmpz_mat_init(AA, M1 - M0, n);
            fmpz_mat_init(BB, n, P1 - P0);
            fmpz_mat_init(CC, M1 - M0, P1 - P0);

            /* Convert to fixed-point matrices. */
            for (i = M0; i < M1; i++)
            {
                if (A_min[i] == WORD_MIN)  /* only zeros in this row */
                    continue;

                for (j = 0; j < n; j++)
                {
                    inexact = arf_get_fmpz_fixed_si(fmpz_mat_entry(AA, i - M0, j),
                        arb_midref(arb_mat_entry(A, i, block_start + j)), A_min[i]);

                    if (inexact)
                    {
                        flint_throw(FLINT_ERROR, "matrix multiplication: bad exponent!\n");
                    }
                }
            }

            for (i = P0; i < P1; i++)
            {
                if (B_min[i] == WORD_MIN)  /* only zeros in this column */
                    continue;

                for (j = 0; j < n; j++)
                {
                    inexact = arf_get_fmpz_fixed_si(fmpz_mat_entry(BB, j, i - P0),
                        arb_midref(arb_mat_entry(B, block_start + j, i)), B_min[i]);

                    if (inexact)
                    {
                        flint_throw(FLINT_ERROR, "matrix multiplication: bad exponent!\n");
                    }
                }
            }

            /* The main multiplication */
            fmpz_mat_mul(CC, AA, BB);
            /* flint_printf("bits %wd %wd %wd\n", fmpz_mat_max_bits(CC),
                        fmpz_mat_max_bits(AA), fmpz_mat_max_bits(BB)); */

            fmpz_mat_clear(AA);
            fmpz_mat_clear(BB);

            arb_init(t);

            /* Add to the result matrix */
            for (i = M0; i < M1; i++)
            {
                for (j = P0; j < P1; j++)
                {
                    *e = A_min[i] + B_min[j];

                    /* The first time we write this Cij */
                    if (block_start == 0)
                    {
                        arb_set_round_fmpz_2exp(arb_mat_entry(C, i, j),
                            fmpz_mat_entry(CC, i - M0, j - P0), e, prec);
                    }
                    else
                    {
                        arb_set_round_fmpz_2exp(t, fmpz_mat_entry(CC, i - M0, j - P0), e, prec);
                        arb_add(arb_mat_entry(C, i, j), arb_mat_entry(C, i, j), t, prec);
                    }
                }
            }
            arb_clear(t);

            fmpz_mat_clear(CC);
        }
    }
}

/* todo: squaring optimizations */
void
arb_mat_mul_block(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong M, N, P;
    slong *A_min, *A_max, *B_min, *B_max;
    short *A_bits, *B_bits;
    slong *A_bot, *B_bot;
    slong block_start, block_end, i, j, bot, top, max_height;
    slong b, A_max_bits, B_max_bits;
    slong min_block_size;
    arb_srcptr t;
    int A_exact, B_exact;
    double A_density, B_density;

    M = arb_mat_nrows(A);
    N = arb_mat_ncols(A);
    P = arb_mat_ncols(B);

    if (N != arb_mat_nrows(B) || M != arb_mat_nrows(C) || P != arb_mat_ncols(C))
    {
        flint_throw(FLINT_ERROR, "arb_mat_mul_block: incompatible dimensions\n");
    }

    if (M == 0 || N == 0 || P == 0 || arb_mat_is_zero(A) || arb_mat_is_zero(B))
    {
        arb_mat_zero(C);
        return;
    }

    if (A == C || B == C)
    {
        arb_mat_t T;
        arb_mat_init(T, M, P);
        arb_mat_mul_block(T, A, B, prec);
        arb_mat_swap_entrywise(T, C);
        arb_mat_clear(T);
        return;
    }

    /* We assume everywhere below that exponents cannot overflow/underflow
       the small fmpz value range. */
    if (!arb_mat_is_lagom(A) || !arb_mat_is_lagom(B))
    {
        arb_mat_mul_classical(C, A, B, prec);
        return;
    }

    /* bottom exponents of A */
    A_bot = flint_malloc(sizeof(slong) * M * N);
    /* minimum bottom exponent in current row */
    A_min = flint_malloc(sizeof(slong) * M);
    /* maximum top exponent in current row */
    A_max = flint_malloc(sizeof(slong) * M);

    B_bot = flint_malloc(sizeof(slong) * N * P);
    B_min = flint_malloc(sizeof(slong) * P);
    B_max = flint_malloc(sizeof(slong) * P);

    /* save space using shorts to store the bit sizes temporarily;
       the block algorithm will not be used at extremely high precision */
    A_bits = flint_malloc(sizeof(short) * M * N);
    B_bits = flint_malloc(sizeof(short) * N * P);

    A_exact = B_exact = 1;
    A_max_bits = B_max_bits = 0;
    A_density = B_density = 0;

    /* Build table of bottom exponents (WORD_MIN signifies a zero),
       and also collect some statistics. */
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            t = arb_mat_entry(A, i, j);
            if (arf_is_zero(arb_midref(t)))
            {
                A_bot[i * N + j] = WORD_MIN;
                A_bits[i * N + j] = 0;
            }
            else
            {
                b = arf_bits(arb_midref(t));
                A_bot[i * N + j] = ARF_EXP(arb_midref(t)) - b;
                A_bits[i * N + j] = b;
                A_max_bits = FLINT_MAX(A_max_bits, b);
                A_density++;
            }
            A_exact = A_exact && mag_is_zero(arb_radref(t));
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < P; j++)
        {
            t = arb_mat_entry(B, i, j);
            if (arf_is_zero(arb_midref(t)))
            {
                B_bot[i * P + j] = WORD_MIN;
                B_bits[i * P + j] = 0;
            }
            else
            {
                b = arf_bits(arb_midref(t));
                B_bot[i * P + j] = ARF_EXP(arb_midref(t)) - b;
                B_bits[i * P + j] = b;
                B_max_bits = FLINT_MAX(B_max_bits, b);
                B_density++;
            }
            B_exact = B_exact && mag_is_zero(arb_radref(t));
        }
    }

    A_density = A_density / (M * N);
    B_density = B_density / (N * P);

    /* Don't shift too far when creating integer block matrices. */
    max_height = 1.25 * FLINT_MIN(prec, FLINT_MAX(A_max_bits, B_max_bits)) + 192;

    /* Avoid block algorithm for extremely high-precision matrices? */
    /* Warning: these cutoffs are completely bogus... */
    if (A_max_bits > 8000 || B_max_bits > 8000 ||
        (A_density < 0.1 && B_density < 0.1 && max_height > 1024))
    {
        flint_free(A_bot);
        flint_free(A_max);
        flint_free(A_min);
        flint_free(B_bot);
        flint_free(B_max);
        flint_free(B_min);
        flint_free(A_bits);
        flint_free(B_bits);
        arb_mat_mul_classical(C, A, B, prec);
        return;
    }

    if (arb_mat_mul_block_min_block_size != 0)
        min_block_size = arb_mat_mul_block_min_block_size;
    else
        min_block_size = 30;

    block_start = 0;
    while (block_start < N)
    {
        /* Find a run of columns of A and rows of B such that the
           bottom exponents differ by at most max_height. */

        block_end = block_start + 1;  /* index is exclusive block_end */

        /* begin with this column of A and row of B */
        for (i = 0; i < M; i++)
        {
            A_max[i] = A_min[i] = A_bot[i * N + block_start];
            A_max[i] += (slong) A_bits[i * N + block_start];
        }

        for (i = 0; i < P; i++)
        {
            B_max[i] = B_min[i] = B_bot[block_start * P + i];
            B_max[i] += (slong) B_bits[block_start * P + i];
        }

        while (block_end < N)
        {
            double size;

            /* End block if memory would be excessive. */
            /* Necessary? */
            /* Should also do initial check above, if C alone is too large. */
            size = (block_end - block_start) * M * (double) A_max_bits;
            size += (block_end - block_start) * P * (double) B_max_bits;
            size += (M * P) * (double) (A_max_bits + B_max_bits);
            size /= 8.0;
            if (size > 2e9)
                goto blocks_built;

            /* check if we can extend with column [block_end] of A */
            for (i = 0; i < M; i++)
            {
                bot = A_bot[i * N + block_end];
                /* zeros are irrelevant */
                if (bot == WORD_MIN || A_max[i] == WORD_MIN)
                    continue;
                top = bot + (slong) A_bits[i * N + block_end];
                /* jump will be too big */
                if (top > A_min[i] + max_height || bot < A_max[i] - max_height)
                    goto blocks_built;
            }

            /* check if we can extend with row [block_end] of B */
            for (i = 0; i < P; i++)
            {
                bot = B_bot[block_end * P + i];
                if (bot == WORD_MIN || B_max[i] == WORD_MIN)
                    continue;
                top = bot + (slong) B_bits[block_end * P + i];
                if (top > B_min[i] + max_height || bot < B_max[i] - max_height)
                    goto blocks_built;
            }

            /* second pass to update the extreme values */
            for (i = 0; i < M; i++)
            {
                bot = A_bot[i * N + block_end];
                top = bot + (slong) A_bits[i * N + block_end];
                if (A_max[i] == WORD_MIN)
                {
                    A_max[i] = top;
                    A_min[i] = bot;
                }
                else if (bot != WORD_MIN)
                {
                    if (bot < A_min[i]) A_min[i] = bot;
                    if (top > A_max[i]) A_max[i] = top;
                }
            }

            for (i = 0; i < P; i++)
            {
                bot = B_bot[block_end * P + i];
                top = bot + (slong) B_bits[block_end * P + i];
                if (B_max[i] == WORD_MIN)
                {
                    B_max[i] = top;
                    B_min[i] = bot;
                }
                else if (bot != WORD_MIN)
                {
                    if (bot < B_min[i]) B_min[i] = bot;
                    if (top > B_max[i]) B_max[i] = top;
                }
            }

            block_end++;
        }

    blocks_built:
        if (block_end - block_start < min_block_size)
        {
            block_end = FLINT_MIN(N, block_start + min_block_size);

            arb_mat_mid_addmul_block_fallback(C, A, B,
                block_start, block_end, prec);
        }
        else
        {
            arb_mat_mid_addmul_block_prescaled(C, A, B,
                block_start, block_end, A_min, B_min, prec);
        }

        block_start = block_end;
    }

    flint_free(A_bot);
    flint_free(A_max);
    flint_free(A_min);
    flint_free(B_bot);
    flint_free(B_max);
    flint_free(B_min);
    flint_free(A_bits);
    flint_free(B_bits);

    /* Radius multiplications */
    if (!A_exact || !B_exact)
    {
        mag_ptr AA, BB;

        /* Shallow (since exponents are small!) mag_struct matrices
           represented by linear arrays; B is transposed to improve locality. */
        AA = flint_malloc(M * N * sizeof(mag_struct));
        BB = flint_malloc(P * N * sizeof(mag_struct));

        if (!A_exact && !B_exact)
        {
            /* (A+ar)(B+br) = AB + (A+ar)br + ar B
                            = AB + A br + ar (B + br) */

            /* A + ar */
            for (i = 0; i < M; i++)
                for (j = 0; j < N; j++)
                {
                    mag_fast_init_set_arf(AA + i * N + j,
                        arb_midref(arb_mat_entry(A, i, j)));
                    mag_add(AA + i * N + j, AA + i * N + j,
                        arb_radref(arb_mat_entry(A, i, j)));
                }

            /* br */
            for (i = 0; i < N; i++)
                for (j = 0; j < P; j++)
                    BB[j * N + i] = *arb_radref(arb_mat_entry(B, i, j));

            _arb_mat_addmul_rad_mag_fast(C, AA, BB, M, N, P);

            /* ar */
            for (i = 0; i < M; i++)
                for (j = 0; j < N; j++)
                    AA[i * N + j] = *arb_radref(arb_mat_entry(A, i, j));

            /* B */
            for (i = 0; i < N; i++)
                for (j = 0; j < P; j++)
                    mag_fast_init_set_arf(BB + j * N + i,
                        arb_midref(arb_mat_entry(B, i, j)));

            _arb_mat_addmul_rad_mag_fast(C, AA, BB, M, N, P);
        }
        else if (A_exact)
        {
            /* A(B+br) = AB + A br */

            for (i = 0; i < M; i++)
                for (j = 0; j < N; j++)
                    mag_fast_init_set_arf(AA + i * N + j,
                        arb_midref(arb_mat_entry(A, i, j)));

            for (i = 0; i < N; i++)
                for (j = 0; j < P; j++)
                    BB[j * N + i] = *arb_radref(arb_mat_entry(B, i, j));

            _arb_mat_addmul_rad_mag_fast(C, AA, BB, M, N, P);
        }
        else
        {
            /* (A+ar)B = AB + ar B */

            for (i = 0; i < M; i++)
                for (j = 0; j < N; j++)
                    AA[i * N + j] = *arb_radref(arb_mat_entry(A, i, j));

            for (i = 0; i < N; i++)
                for (j = 0; j < P; j++)
                    mag_fast_init_set_arf(BB + j * N + i,
                        arb_midref(arb_mat_entry(B, i, j)));

            _arb_mat_addmul_rad_mag_fast(C, AA, BB, M, N, P);
        }

        flint_free(AA);
        flint_free(BB);
    }
}

