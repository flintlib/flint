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
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

void
fmpq_mat_get_fmpz_mat_matwise(fmpz_mat_t num, fmpz_t den, const fmpq_mat_t mat)
{
    fmpz_t t, lcm;
    long i, j;

    if (fmpq_mat_is_empty(mat))
        return;

    fmpz_init(t);
    fmpz_init(lcm);
    fmpz_one(lcm);

    /* Compute common denominator of matrix */
    for (i = 0; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            fmpz_lcm(lcm, lcm, fmpq_mat_entry_den(mat, i, j));

    fmpz_set(den, lcm);

    for (i = 0; i < mat->r; i++)
    {
        for (j = 0; j < mat->c; j++)
        {
            /* Rescale numerators */
            if (fmpz_is_one(lcm))
            {
                fmpz_set(fmpz_mat_entry(num, i, j),
                             fmpq_mat_entry_num(mat, i, j));
            }
            else
            {
                fmpz_divexact(t, lcm, fmpq_mat_entry_den(mat, i, j));
                fmpz_mul(fmpz_mat_entry(num, i, j),
                         fmpq_mat_entry_num(mat, i, j), t);
            }
        }
    }

    fmpz_clear(t);
    fmpz_clear(lcm);
}
