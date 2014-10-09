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

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_mat.h"

int
fmpz_mat_is_one(const fmpz_mat_t mat)
{
    slong i, j;

    if (mat->r == 0 || mat->c == 0)
        return 1;

    for (i = 0; i < mat->r; i++)
    {
        for (j = 0; j < mat->c; j++)
        {
            if (fmpz_cmp_ui(fmpz_mat_entry(mat, i, j), i == j) != 0)
            {
                return 0;
            }
        }
    }

    return 1;
}
