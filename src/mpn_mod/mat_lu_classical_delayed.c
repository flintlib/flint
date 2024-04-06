/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

/* for wide add_ssss.... macros. todo; these ought to be provided
   everywhere */
#if FLINT_BITS == 64 && defined(__AVX2__)
#include "crt_helpers.h"
#endif

/* todo: optimize for when 2n rather than 2n+1 limbs suffice */
int
mpn_mod_mat_lu_classical_delayed(slong * res_rank, slong * P, gr_mat_t A, const gr_mat_t A_in, int rank_check, gr_ctx_t ctx)
{
    mp_limb_t d[MPN_MOD_MAX_LIMBS];
    mp_limb_t e[MPN_MOD_MAX_LIMBS];
    mp_limb_t f[MPN_MOD_MAX_LIMBS];
    mp_ptr * a;
    mp_ptr tmprow;
    slong n = MPN_MOD_CTX_NLIMBS(ctx);
    slong i, j, nrows, ncols, rank, row, col, pivot_row, tmp_index;
    int status = GR_SUCCESS;
    mp_ptr tmp_ptr, b;
    TMP_INIT;

    nrows = A->r;
    ncols = A->c;

    if (nrows == 0 || ncols == 0)
    {
        *res_rank = 0;
        return GR_SUCCESS;
    }

    a = (mp_ptr *) A->rows;

    if (A != A_in)
    {
        for (i = 0; i < nrows; i++)
            flint_mpn_copyi(a[i], A_in->rows[i], n * ncols);
    }

    rank = row = col = 0;

    for (i = 0; i < nrows; i++)
        P[i] = i;

    TMP_START;
    b = TMP_ALLOC((2 * n + 1) * sizeof(mp_limb_t) * (nrows + 1) * ncols);
    tmprow = b + (2 * n + 1) * (nrows * ncols);

#define UNREDUCED(ii, jj) (b + (2 * n + 1) * ((ii) * ncols + (jj)))
#define REDUCED(ii, jj) (a[ii] + ((jj) * n))
#define TMPROW(ii) (tmprow + (2 * n + 1) * (ii))

    flint_mpn_zero(tmprow, (2 * n + 1) * ncols);

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            flint_mpn_copyi(UNREDUCED(i, j), REDUCED(i, j), n);
            flint_mpn_zero(UNREDUCED(i, j) + n, n + 1);
        }
    }

    while (row < nrows && col < ncols)
    {
        /* reduce current column */
        /* can be skipped on the first iteration */
        if (col != 0)
            for (j = row; j < nrows; j++)
                mpn_mod_set_mpn(REDUCED(j, col), UNREDUCED(j, col), 2 * n + 1, ctx);

        pivot_row = -1;
        for (i = row; i < nrows; i++)
        {
            if (!flint_mpn_zero_p(REDUCED(i, col), n) != 0)
            {
                pivot_row = i;
                break;
            }
        }

        /* There is certainly no nonzero pivot element. */
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
            tmp_ptr = a[pivot_row];
            a[pivot_row] = a[row];
            a[row] = tmp_ptr;

            tmp_index = P[pivot_row];
            P[pivot_row] = P[row];
            P[row] = tmp_index;

            /* swap rows in unreduced submatrix, and reduce new pivot row */
            for (j = col + 1; j < ncols; j++)
            {
                mpn_mod_set_mpn(REDUCED(row, j), UNREDUCED(pivot_row, j), 2 * n + 1, ctx);
                flint_mpn_copyi(UNREDUCED(pivot_row, j), UNREDUCED(row, j), 2 * n + 1);
            }
        }
        else if (row != 0)
        {
            /* Reduce current pivot row. */
            for (j = col + 1; j < ncols; j++)
                mpn_mod_set_mpn(REDUCED(row, j), UNREDUCED(row, j), 2 * n + 1, ctx);
        }

        rank++;

        /* Eliminate remaining submatrix. */

        /* Must be able to invert pivot element. */
        status = mpn_mod_inv(d, REDUCED(row, col), ctx);
        if (status != GR_SUCCESS)
            break;

        for (i = row + 1; i < nrows; i++)
        {
            mpn_mod_mul(e, REDUCED(i, col), d, ctx);
            mpn_mod_neg(f, e, ctx);

#if defined(add_sssssaaaaaaaaaa)
            if (n == 2)
            {
                for (j = col + 1; j < ncols; j++)
                {
                    mp_limb_t t[4];
                    FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], REDUCED(row, j)[1], REDUCED(row, j)[0], f[1], f[0]);
                    add_sssssaaaaaaaaaa(UNREDUCED(i, j)[4], UNREDUCED(i, j)[3], UNREDUCED(i, j)[2], UNREDUCED(i, j)[1], UNREDUCED(i, j)[0],
                                        UNREDUCED(i, j)[4], UNREDUCED(i, j)[3], UNREDUCED(i, j)[2], UNREDUCED(i, j)[1], UNREDUCED(i, j)[0],
                                        0, t[3], t[2], t[1], t[0]);
                }

                REDUCED(i, col)[0] = 0;
                REDUCED(i, col)[1] = 0;
                REDUCED(i, rank - 1)[0] = e[0];
                REDUCED(i, rank - 1)[1] = e[1];
            }
            else
#endif
            {
                if (col + 1 < ncols)
                {
#if 1
                    for (j = col + 1; j < ncols; j++)
                        flint_mpn_mul_n(TMPROW(j), REDUCED(row, j), f, n);

                    mpn_add_n(UNREDUCED(i, col + 1), UNREDUCED(i, col + 1), TMPROW(col + 1), (2 * n + 1) * (ncols - col - 1));
#else
                    for (j = col + 1; j < ncols; j++)
                    {
                        flint_mpn_mul_n(TMPROW(0), REDUCED(row, j), f, n);
                        mpn_add_n(UNREDUCED(i, j), UNREDUCED(i, j), TMPROW(0), 2 * n + 1);
                    }
#endif
                }

                flint_mpn_zero(REDUCED(i, col), n);
                flint_mpn_copyi(REDUCED(i, rank - 1), e, n);
            }
        }
        row++;
        col++;
    }

    *res_rank = rank;

    TMP_END;
    return status;
}
