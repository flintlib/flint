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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

static int _pivot(mp_limb_t ** rows, long n, long from_row, long in_column)
{
    mp_limb_t * tmp;
    long j;

    if (rows[from_row][in_column] != 0)
        return 1;

    for (j = from_row + 1; j < n; j++)
    {
        if (rows[j][in_column] != 0L)
        {
            tmp = rows[j];
            rows[j] = rows[from_row];
            rows[from_row] = tmp;
            return -1;
        }
    }
    return 0;
}

long _nmod_mat_rowreduce(nmod_mat_t mat, int options)
{
    mp_limb_t ** a;
    long m, n;
    nmod_t mod;

    long j, rank;
    int sign = 1;
    int det_sign;

    long pivot_row;
    long pivot_col;
    long length;

    mp_limb_t d, e;

    m = mat->r;
    n = mat->c;
    a = mat->rows;
    mod = mat->mod;

    rank = 0L;
    pivot_row = 0L;
    pivot_col = 0L;

    while (pivot_row < m && pivot_col < n)
    {
        det_sign = _pivot(a, m, pivot_row, pivot_col);

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
            {
                _nmod_vec_scalar_addmul(a[j] + pivot_col + 1,
                                        a[pivot_row] + pivot_col + 1,
                                        length, mod, e);
            }

            if (options & ROWREDUCE_CLEAR_LOWER)
                a[j][pivot_col] = 0L;
            else
                a[j][pivot_col] = d; /* LU decomposition */

        }

        pivot_row++;
        pivot_col++;
    }

    return rank * sign;
}
