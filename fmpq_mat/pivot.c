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

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

int
fmpq_mat_pivot(len_t * perm, fmpq_mat_t mat, len_t r, len_t c)
{
    len_t t, j;
    fmpq * u;

    if (!fmpq_is_zero(fmpq_mat_entry(mat, r, c)))
        return 1;

    for (j = r + 1; j < mat->r; j++)
    {
        if (!fmpq_is_zero(fmpq_mat_entry(mat, j, c)))
        {
            if (perm)
            {
                t = perm[j];
                perm[j] = perm[r];
                perm[r] = t;
            }

            u = mat->rows[j];
            mat->rows[j] = mat->rows[r];
            mat->rows[r] = u; 
            return -1;
        }
    }
    return 0;
}
