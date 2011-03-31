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


long _nmod_mat_rowreduce_r(nmod_mat_t mat, int options)
{
    mp_limb_t d, e, **a;
    nmod_t mod;
    long j, m, n, rank, length, pivot_row, pivot_col;
    int det_sign, sign = 1;

    m = mat->r;
    n = mat->c;
    a = mat->rows;
    mod = mat->mod;

    rank = pivot_row = pivot_col = 0;

    while (pivot_row < m && pivot_col < n)
    {
        det_sign = _nmod_mat_pivot(a, m, pivot_row, pivot_col);

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

        for (j = (options & ROWREDUCE_FULL ? 0 : pivot_row + 1); j < m; j++)
        {
            if (j == pivot_row)
                continue;

            d = a[pivot_row][pivot_col];
            e = a[j][pivot_col];

            d = n_invmod(d, mod.n);
            d = n_mulmod2_preinv(e, d, mod.n, mod.ninv);
            e = nmod_neg(d, mod);

            /* Length argument may not be zero */
            length = n - pivot_col - 1;
            if (length > 0)
                _nmod_vec_scalar_addmul_nmod(a[j] + pivot_col + 1,
                    a[pivot_row] + pivot_col + 1, length, e, mod);

            if (options & ROWREDUCE_CLEAR_LOWER)
                a[j][pivot_col] = 0UL;
            else
                a[j][pivot_col] = d; /* LU decomposition */
        }

        pivot_row++;
        pivot_col++;
    }

    return rank * sign;
}
