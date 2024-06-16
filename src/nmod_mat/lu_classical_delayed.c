/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

static ulong
nmod_set_uiuiui(ulong s2, ulong s1, ulong s0, nmod_t mod)
{
    NMOD_RED(s2, s2, mod);
    NMOD_RED3(s0, s2, s1, s0, mod);
    return s0;
}

slong
nmod_mat_lu_classical_delayed_1(slong * P, nmod_mat_t A, int rank_check)
{
    ulong d, e, f, **a;
    nmod_t mod;
    slong i, j, nrows, ncols, rank, row, col, pivot_row, tmp_index;
    nn_ptr tnn_ptr;

    nrows = A->r;
    ncols = A->c;
    a = A->rows;
    mod = A->mod;

    rank = row = col = 0;

    for (i = 0; i < nrows; i++)
        P[i] = i;

    while (row < nrows && col < ncols)
    {
        /* reduce current column */
        /* can be skipped on the first iteration */
        if (col != 0)
            for (j = row; j < nrows; j++)
                NMOD_RED(a[j][col], a[j][col], mod);

        pivot_row = -1;
        for (i = row; i < nrows; i++)
        {
            if (a[i][col] != 0)
            {
                pivot_row = i;
                break;
            }
        }

        if (pivot_row == -1)
        {
            if (rank_check)
            {
                rank = 0;
                break;
            }

            col++;
            continue;
        }

        /* swap rows */
        if (pivot_row != row)
        {
            tnn_ptr = a[pivot_row];
            a[pivot_row] = a[row];
            a[row] = tnn_ptr;

            tmp_index = P[pivot_row];
            P[pivot_row] = P[row];
            P[row] = tmp_index;
        }

        /* reduce current pivot row */
        if (col != 0)
            for (j = col + 1; j < ncols; j++)
                NMOD_RED(a[row][j], a[row][j], mod);

        rank++;

        /* eliminate remaining submatrix */
        d = nmod_inv(a[row][col], mod);
        for (i = row + 1; i < nrows; i++)
        {
            e = nmod_mul(a[i][col], d, mod);
            f = nmod_neg(e, mod);

            for (j = col + 1; j + 4 < ncols; j += 4)
            {
                ulong x0, x1, x2, x3;
                x0 = a[row][j + 0];
                x1 = a[row][j + 1];
                x2 = a[row][j + 2];
                x3 = a[row][j + 3];
                a[i][j + 0] += x0 * f;
                a[i][j + 1] += x1 * f;
                a[i][j + 2] += x2 * f;
                a[i][j + 3] += x3 * f;
            }

            for ( ; j < ncols; j++)
                a[i][j] += a[row][j] * f;

            a[i][col] = 0;
            a[i][rank - 1] = e;
        }
        row++;
        col++;
    }

    return rank;
}

slong
nmod_mat_lu_classical_delayed_2(slong * P, nmod_mat_t A, int rank_check)
{
    ulong d, e, f, **a;
    nmod_t mod;
    slong i, j, nrows, ncols, rank, row, col, pivot_row, tmp_index;
    nn_ptr tnn_ptr;
    nn_ptr b;
    TMP_INIT;

    nrows = A->r;
    ncols = A->c;
    a = A->rows;
    mod = A->mod;

    rank = row = col = 0;

    for (i = 0; i < nrows; i++)
        P[i] = i;

    TMP_START;
    b = TMP_ALLOC(2 * sizeof(ulong) * nrows * ncols);

#define UNREDUCED_LO(ii, jj) b[2 * ((ii) * ncols + jj)]
#define UNREDUCED_HI(ii, jj) b[2 * ((ii) * ncols + jj) + 1]

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            UNREDUCED_LO(i, j) = a[i][j];
            UNREDUCED_HI(i, j) = 0;
        }
    }

    while (row < nrows && col < ncols)
    {
        /* reduce current column */
        /* can be skipped on the first iteration */
        if (col != 0)
            for (j = row; j < nrows; j++)
                NMOD2_RED2(a[j][col], UNREDUCED_HI(j, col), UNREDUCED_LO(j, col), mod);

        pivot_row = -1;
        for (i = row; i < nrows; i++)
        {
            if (a[i][col] != 0)
            {
                pivot_row = i;
                break;
            }
        }

        if (pivot_row == -1)
        {
            if (rank_check)
            {
                rank = 0;
                break;
            }

            col++;
            continue;
        }

        /* swap rows */
        if (pivot_row != row)
        {
            tnn_ptr = a[pivot_row];
            a[pivot_row] = a[row];
            a[row] = tnn_ptr;

            tmp_index = P[pivot_row];
            P[pivot_row] = P[row];
            P[row] = tmp_index;

            /* swap rows in unreduced submatrix, and reduce new pivot row */
            for (j = col + 1; j < ncols; j++)
            {
                ulong hi, lo;
                lo = UNREDUCED_LO(row, j);
                hi = UNREDUCED_HI(row, j);
                NMOD2_RED2(a[row][j], UNREDUCED_HI(pivot_row, j), UNREDUCED_LO(pivot_row, j), mod);
                UNREDUCED_LO(pivot_row, j) = lo;
                UNREDUCED_HI(pivot_row, j) = hi;
            }
        }
        else if (row != 0)
        {
            /* reduce current pivot row */
            for (j = col + 1; j < ncols; j++)
                NMOD2_RED2(a[row][j], UNREDUCED_HI(row, j), UNREDUCED_LO(row, j), mod);
        }

        rank++;

        /* eliminate remaining submatrix */
        d = nmod_inv(a[row][col], mod);
        for (i = row + 1; i < nrows; i++)
        {
            e = nmod_mul(a[i][col], d, mod);
            f = nmod_neg(e, mod);

            if (mod.n <= UWORD(1) << (FLINT_BITS / 2))
            {
                for (j = col + 1; j + 4 < ncols; j += 4)
                {
                    ulong x0, x1, x2, x3;
                    x0 = a[row][j + 0] * f;
                    x1 = a[row][j + 1] * f;
                    x2 = a[row][j + 2] * f;
                    x3 = a[row][j + 3] * f;
                    add_ssaaaa(UNREDUCED_HI(i, j + 0), UNREDUCED_LO(i, j + 0),
                               UNREDUCED_HI(i, j + 0), UNREDUCED_LO(i, j + 0), 0, x0);
                    add_ssaaaa(UNREDUCED_HI(i, j + 1), UNREDUCED_LO(i, j + 1),
                               UNREDUCED_HI(i, j + 1), UNREDUCED_LO(i, j + 1), 0, x1);
                    add_ssaaaa(UNREDUCED_HI(i, j + 2), UNREDUCED_LO(i, j + 2),
                               UNREDUCED_HI(i, j + 2), UNREDUCED_LO(i, j + 2), 0, x2);
                    add_ssaaaa(UNREDUCED_HI(i, j + 3), UNREDUCED_LO(i, j + 3),
                               UNREDUCED_HI(i, j + 3), UNREDUCED_LO(i, j + 3), 0, x3);
                }

                for ( ; j < ncols; j++)
                {
                    ulong hi, lo;
                    hi = 0;
                    lo = a[row][j] * f;
                    add_ssaaaa(UNREDUCED_HI(i, j), UNREDUCED_LO(i, j),
                               UNREDUCED_HI(i, j), UNREDUCED_LO(i, j), hi, lo);
                }
            }
            else
            {
                for (j = col + 1; j + 4 < ncols; j += 4)
                {
                    ulong x0, x1, x2, x3, h0, h1, h2, h3;
                    umul_ppmm(h0, x0, a[row][j + 0], f);
                    umul_ppmm(h1, x1, a[row][j + 1], f);
                    umul_ppmm(h2, x2, a[row][j + 2], f);
                    umul_ppmm(h3, x3, a[row][j + 3], f);
                    add_ssaaaa(UNREDUCED_HI(i, j + 0), UNREDUCED_LO(i, j + 0),
                               UNREDUCED_HI(i, j + 0), UNREDUCED_LO(i, j + 0), h0, x0);
                    add_ssaaaa(UNREDUCED_HI(i, j + 1), UNREDUCED_LO(i, j + 1),
                               UNREDUCED_HI(i, j + 1), UNREDUCED_LO(i, j + 1), h1, x1);
                    add_ssaaaa(UNREDUCED_HI(i, j + 2), UNREDUCED_LO(i, j + 2),
                               UNREDUCED_HI(i, j + 2), UNREDUCED_LO(i, j + 2), h2, x2);
                    add_ssaaaa(UNREDUCED_HI(i, j + 3), UNREDUCED_LO(i, j + 3),
                               UNREDUCED_HI(i, j + 3), UNREDUCED_LO(i, j + 3), h3, x3);
                }

                for ( ; j < ncols; j++)
                {
                    ulong hi, lo;
                    umul_ppmm(hi, lo, a[row][j], f);
                    add_ssaaaa(UNREDUCED_HI(i, j), UNREDUCED_LO(i, j),
                               UNREDUCED_HI(i, j), UNREDUCED_LO(i, j), hi, lo);
                }
            }

            a[i][col] = 0;
            a[i][rank - 1] = e;
        }
        row++;
        col++;
    }

    TMP_END;
    return rank;
}

slong
nmod_mat_lu_classical_delayed_3(slong * P, nmod_mat_t A, int rank_check)
{
    ulong d, e, f, **a;
    nmod_t mod;
    slong i, j, nrows, ncols, rank, row, col, pivot_row, tmp_index;
    nn_ptr tnn_ptr;
    nn_ptr b;
    TMP_INIT;

    nrows = A->r;
    ncols = A->c;
    a = A->rows;
    mod = A->mod;

    rank = row = col = 0;

    for (i = 0; i < nrows; i++)
        P[i] = i;

    TMP_START;
    b = TMP_ALLOC(3 * sizeof(ulong) * nrows * ncols);

#define UNREDUCED3_L0(ii, jj) b[3 * ((ii) * ncols + jj)]
#define UNREDUCED3_L1(ii, jj) b[3 * ((ii) * ncols + jj) + 1]
#define UNREDUCED3_L2(ii, jj) b[3 * ((ii) * ncols + jj) + 2]

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            UNREDUCED3_L0(i, j) = a[i][j];
            UNREDUCED3_L1(i, j) = 0;
            UNREDUCED3_L2(i, j) = 0;
        }
    }

    while (row < nrows && col < ncols)
    {
        /* reduce current column */
        /* can be skipped on the first iteration */
        if (col != 0)
            for (j = row; j < nrows; j++)
                a[j][col] = nmod_set_uiuiui(UNREDUCED3_L2(j, col),
                    UNREDUCED3_L1(j, col),
                    UNREDUCED3_L0(j, col), mod);

        pivot_row = -1;
        for (i = row; i < nrows; i++)
        {
            if (a[i][col] != 0)
            {
                pivot_row = i;
                break;
            }
        }

        if (pivot_row == -1)
        {
            if (rank_check)
            {
                rank = 0;
                break;
            }

            col++;
            continue;
        }

        /* swap rows */
        if (pivot_row != row)
        {
            tnn_ptr = a[pivot_row];
            a[pivot_row] = a[row];
            a[row] = tnn_ptr;

            tmp_index = P[pivot_row];
            P[pivot_row] = P[row];
            P[row] = tmp_index;

            /* swap rows in unreduced submatrix, and reduce new pivot row */
            for (j = col + 1; j < ncols; j++)
            {
                ulong t2, t1, t0;
                t0 = UNREDUCED3_L0(row, j);
                t1 = UNREDUCED3_L1(row, j);
                t2 = UNREDUCED3_L2(row, j);

                a[row][j] = nmod_set_uiuiui(UNREDUCED3_L2(pivot_row, j),
                            UNREDUCED3_L1(pivot_row, j),
                            UNREDUCED3_L0(pivot_row, j), mod);

                UNREDUCED3_L0(pivot_row, j) = t0;
                UNREDUCED3_L1(pivot_row, j) = t1;
                UNREDUCED3_L2(pivot_row, j) = t2;
            }
        }
        else if (row != 0)
        {
            /* reduce current pivot row */
            for (j = col + 1; j < ncols; j++)
                a[row][j] = nmod_set_uiuiui(UNREDUCED3_L2(row, j),
                            UNREDUCED3_L1(row, j),
                            UNREDUCED3_L0(row, j), mod);
        }

        rank++;

        /* eliminate remaining submatrix */
        d = nmod_inv(a[row][col], mod);
        for (i = row + 1; i < nrows; i++)
        {
            e = nmod_mul(a[i][col], d, mod);
            f = nmod_neg(e, mod);

            for (j = col + 1; j < ncols; j++)
            {
                ulong hi, lo;
                umul_ppmm(hi, lo, a[row][j], f);
                add_sssaaaaaa(UNREDUCED3_L2(i, j), UNREDUCED3_L1(i, j), UNREDUCED3_L0(i, j),
                              UNREDUCED3_L2(i, j), UNREDUCED3_L1(i, j), UNREDUCED3_L0(i, j),
                              0, hi, lo);
            }

            a[i][col] = 0;
            a[i][rank - 1] = e;
        }
        row++;
        col++;
    }

    TMP_END;
    return rank;
}

slong
nmod_mat_lu_classical_delayed(slong * P, nmod_mat_t A, int rank_check)
{
    slong nrows, ncols;

    nrows = A->r;
    ncols = A->c;
    const dot_params_t params = _nmod_vec_dot_params(FLINT_MIN(nrows, ncols), A->mod);

    // TODO cases to re-examine after dot product changes?
    if (params.method <= _DOT1)
        return nmod_mat_lu_classical_delayed_1(P, A, rank_check);
    else if (params.method <= _DOT2)
        return nmod_mat_lu_classical_delayed_2(P, A, rank_check);
    else
        return nmod_mat_lu_classical_delayed_3(P, A, rank_check);
}
