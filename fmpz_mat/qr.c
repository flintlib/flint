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

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_mat.h"

void
fmpz_mat_qr(d_mat_t Q, d_mat_t R, const fmpz_mat_t A)
{
    slong i, k, j, r = A->r, c = A->c;

    if (Q->r != A->r || Q->c != A->c || R->r != A->r || R->c != A->r)
    {
        flint_printf("Exception (fmpz_mat_qr). Incompatible dimensions.\n");
        abort();
    }

    if (A->r == 0)
    {
        return;
    }

    fmpz_mat_get_d_mat(Q, A);

    for (k = 0; k < r; k++)
    {
        d_mat_entry(R, k, k) = _d_vec_norm(Q->rows[k], c);
        d_mat_entry(R, k, k) = sqrt(d_mat_entry(R, k, k));
        for (i = 0; i < c; i++)
            if (d_mat_entry(R, k, k) != 0)
                d_mat_entry(Q, k, i) =
                    d_mat_entry(Q, k, i) / d_mat_entry(R, k, k);
        for (j = k + 1; j < r; j++)
        {
            d_mat_entry(R, j, k) = _d_vec_dot(Q->rows[k], Q->rows[j], c);
            for (i = 0; i < c; i++)
            {
                d_mat_entry(Q, j, i) -=
                    d_mat_entry(R, j, k) * d_mat_entry(Q, k, i);
            }
        }
    }
}
