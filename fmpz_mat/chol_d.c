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
fmpz_mat_chol_d(d_mat_t R, const fmpz_mat_t A)
{
    slong i, k, j, r = A->r;

    if (A->r != A->c || R->r != A->r || R->c != A->c)
    {
        flint_printf
            ("Exception (fmpz_mat_chol_d). Incompatible dimensions.\n");
        abort();
    }

    if (A->r == 0)
    {
        return;
    }

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < i + 1; j++)
        {
            double s = 0;
            for (k = 0; k < j; k++)
            {
                s += d_mat_entry(R, i, k) * d_mat_entry(R, j, k);
            }
            if (i == j)
                d_mat_entry(R, i, j) =
                    sqrt(fmpz_get_d(fmpz_mat_entry(A, i, i)) - s);
            else
                d_mat_entry(R, i, j) =
                    (fmpz_get_d(fmpz_mat_entry(A, i, j)) - s) / d_mat_entry(R,
                                                                            j,
                                                                            j);
        }
    }
}
