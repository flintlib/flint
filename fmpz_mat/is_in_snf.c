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

    Copyright (C) 2014 Alex J. Best

******************************************************************************/

#include "fmpz_mat.h"

int fmpz_mat_is_in_snf(const fmpz_mat_t A)
{
    slong i, j;
    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            if (i == j)
            {
                if (fmpz_sgn(fmpz_mat_entry(A, i, i)) < 0)
                    return 0;
                if (i > 0)
                {
                    if (!fmpz_is_zero(fmpz_mat_entry(A, i, i)) &&
                            fmpz_is_zero(fmpz_mat_entry(A, i - 1, i - 1)))
                        return 0;
                    if (!fmpz_divisible(fmpz_mat_entry(A, i, i),
                            fmpz_mat_entry(A, i - 1, i - 1)))
                        return 0;
                }
            }
            else if (!fmpz_is_zero(fmpz_mat_entry(A, i, j)))
                return 0;
        }
    }

    return 1;
}
