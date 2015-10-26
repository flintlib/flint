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

    Copyright (C) 2015 Tommy Hofmann

******************************************************************************/

#include "fmpz_mat.h"

slong
fmpz_mat_howell_form_mod(fmpz_mat_t A, const fmpz_t mod)
{
    slong i, j, n;
    slong k;

    n = A->r;
    k = n;
    fmpz_mat_strong_echelon_form_mod(A, mod);

    for (i = 0; i < n; i++)
    {
        if (fmpz_mat_is_zero_row(A, i))
        {
            k--;
            for (j = i + 1; j < n; j++)
            {
                if (!fmpz_mat_is_zero_row(A, j))
                {
                    fmpz_mat_swap_rows(A, NULL, i, j);
                    j = n;
                    k++;
                }
            }
        }
    }
    return k;
}

