/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "fmpz.h"


long _nmod_mat_rowreduce_2_extended(nmod_mat_t mat, int options)
{
    mp_limb_t d, e, **a;
    nmod_t mod;
    long j, k, m, n, rank, length, pivot_row, pivot_col;
    int det_sign, sign = 1;

    m = mat->r;
    n = mat->c / 2;
    a = mat->rows;
    mod = mat->mod;

    rank = pivot_row = pivot_col = 0;

    while (pivot_row < m && pivot_col < n)
    {
        /* Reduce pivot column */
        for (j = pivot_row; j < m; j++)
        {
            NMOD2_RED2(a[j][2*pivot_col], a[j][2*pivot_col+1],
                a[j][2*pivot_col], mod);
            a[j][2*pivot_col+1] = 0;
        }

        det_sign = _nmod_mat_pivot(a, m, pivot_row, 2*pivot_col);

        /* Reduce pivot row */
        for (k = pivot_col + 1; k < n; k++)
        {
            NMOD2_RED2(a[pivot_row][2*k], a[pivot_row][2*k+1],
                a[pivot_row][2*k], mod);
            a[pivot_row][2*k+1] = 0;
        }

        if (!det_sign)
        {
            if (options & ROWREDUCE_FAST_ABORT)
            {
                rank = 0;
                break;
            }
            pivot_col++;
            continue;
        }

        sign *= det_sign;
        rank++;

        for (j = pivot_row + 1; j < m; j++)
        {
            if (j == pivot_row)
                continue;

            d = a[pivot_row][2*pivot_col];
            e = a[j        ][2*pivot_col];

            d = n_invmod(d, mod.n);
            d = n_mulmod2_preinv(e, d, mod.n, mod.ninv);
            e = nmod_neg(d, mod);

            if (1)  /* for loop version */
            {
                for (k = pivot_col + 1; k < n; k++)
                {
                    mp_limb_t u, dh, f, fh;
                    u = a[j][2*k];
                    dh = a[j][2*k+1];
                    f = a[pivot_row][2*k];
                    fh = a[pivot_row][2*k+1];
                    umul_ppmm(fh, f, f, e);
                    add_ssaaaa(dh, u, dh, u, fh, f);
                    a[j][2*k] = u;
                    a[j][2*k+1] = dh;
                }
            }
            else   /* mpn_addmul_1 version */
            {
                length = n - pivot_col - 1;
                if (length > 0)
                {
                    long offset = pivot_col + 1;
                    mpn_addmul_1(a[j] + 2*offset, \
                        a[pivot_row] + 2*offset, 2*length, e);
                }
            }

            if (options & ROWREDUCE_CLEAR_LOWER)
            {
                a[j][2*pivot_col] = 0UL;
                a[j][2*pivot_col+1] = 0UL;
            }
            else
            {
                a[j][2*pivot_col] = d;
                a[j][2*pivot_col+1] = 0UL;
            }
        }

        pivot_row++;
        pivot_col++;
    }

    return rank * sign;
}

long _nmod_mat_rowreduce_2(nmod_mat_t mat, int options)
{
    nmod_mat_t mat2;
    long m = mat->r;
    long n = mat->c;
    long j, k, rank;

    if (FLINT_MIN(m,n) < 1)
        return 0;

    nmod_mat_init(mat2, m, 2*n, mat->mod.n);
    for (j = 0; j < m; j++)
    {
        for (k = 0; k < n; k++)
        {
            mat2->rows[j][2*k] = mat->rows[j][k];
            mat2->rows[j][2*k+1] = 0;
        }
    }

    rank = _nmod_mat_rowreduce_2_extended(mat2, options);

    for (j = 0; j < m; j++)
        mat->rows[j] = mat->entries + 
            (((mat2->rows[j] - mat2->entries) / (2*n)) * mat->c);

    for (j = 0; j < m; j++)
        for (k = 0; k < n; k++)
            mat->rows[j][k] = mat2->rows[j][2*k];

    nmod_mat_clear(mat2);
    return rank;
}
