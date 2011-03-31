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
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


long _nmod_mat_rowreduce_1(nmod_mat_t mat, int options)
{
    mp_limb_t ** a, d, e;
    long j, k, m, n, rank, pivot_row, pivot_col;
    int det_sign, sign = 1;
    nmod_t mod;

    a = mat->rows;
    m = mat->r;
    n = mat->c;
    mod = mat->mod;

    rank = pivot_row = pivot_col = 0;

    while (pivot_row < m && pivot_col < n)
    {
        /* Reduce pivot column */
        for (j = pivot_row; j < m; j++)
            NMOD_RED(a[j][pivot_col], a[j][pivot_col], mod);

        det_sign = _nmod_mat_pivot(a, m, pivot_row, pivot_col);

        if (!det_sign)
        {
            if (options & ROWREDUCE_FAST_ABORT)
            {
                rank = 0L;
                break;
            }
            pivot_col++;
            continue;
        }

        sign *= det_sign;
        rank++;

        /* Reduce pivot column */
        for (k = pivot_col; k < n; k++)
        {
            d = a[pivot_row][k];
            NMOD_RED(d, d, mod);
            a[pivot_row][k] = d;
        }

        for (j = pivot_row + 1; j < m; j++)
        {
            if (j == pivot_row)
                continue;

            d = a[pivot_row][pivot_col];
            e = a[j][pivot_col];
            d = n_invmod(d, mod.n);
            d = n_mulmod2_preinv(e, d, mod.n, mod.ninv);
            e = nmod_neg(d, mod);

            if (1)  /* for loop version */
            {
                for (k = pivot_col + 1; k < n; k++)
                    a[j][k] += a[pivot_row][k] * e;
            }
            else    /* mpn_addmul_1 version */
            {
                mpn_addmul_1(a[j]+pivot_col+1,
                    a[pivot_row]+pivot_col+1, n-pivot_col-1, e);
            }

            if (options & ROWREDUCE_CLEAR_LOWER)
                a[j][pivot_col] = 0L;
            else
                a[j][pivot_col] = d;
        }

        pivot_row++;
        pivot_col++;
    }

    return rank * sign;
}
