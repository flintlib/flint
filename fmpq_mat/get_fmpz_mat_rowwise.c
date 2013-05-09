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
_fmpq_mat_get_fmpz_mat_rowwise(fmpz_mat_struct ** num, fmpz * den,
                        const fmpq_mat_struct ** mat, len_t n)
{
    fmpz_t t, lcm;
    len_t i, j, k;

    if (fmpq_mat_is_empty(mat[0]))
        return;

    fmpz_init(t);
    fmpz_init(lcm);

    for (i = 0; i < mat[0]->r; i++)
    {
        /* Compute common denominator of row */
        fmpz_set(lcm, fmpq_mat_entry_den(mat[0], i, 0));

        for (k = 0; k < n; k++)
            for (j = (k == 0); j < mat[k]->c; j++)
                fmpz_lcm(lcm, lcm, fmpq_mat_entry_den(mat[k], i, j));

        if (den != NULL)
            fmpz_set(den + i, lcm);

        for (k = 0; k < n; k++)
        {
            /* Rescale numerators in row */
            if (fmpz_is_one(lcm))
            {
                for (j = 0; j < mat[k]->c; j++)
                    fmpz_set(fmpz_mat_entry(num[k], i, j),
                             fmpq_mat_entry_num(mat[k], i, j));
            }
            else
            {
                for (j = 0; j < mat[k]->c; j++)
                {
                    fmpz_divexact(t, lcm, fmpq_mat_entry_den(mat[k], i, j));
                    fmpz_mul(fmpz_mat_entry(num[k], i, j),
                             fmpq_mat_entry_num(mat[k], i, j), t);
                }
            }
        }
    }

    fmpz_clear(t);
    fmpz_clear(lcm);
}

void
fmpq_mat_get_fmpz_mat_rowwise(fmpz_mat_t num, fmpz * den, const fmpq_mat_t mat)
{
    _fmpq_mat_get_fmpz_mat_rowwise(&num, den, &mat, 1);
}

void
fmpq_mat_get_fmpz_mat_rowwise_2(fmpz_mat_t num, fmpz_mat_t num2, fmpz * den,
                                   const fmpq_mat_t mat, const fmpq_mat_t mat2)
{
    fmpz_mat_struct * nums[2];
    fmpq_mat_struct * mats[2];

    nums[0] = num;
    nums[1] = num2;
    mats[0] = (fmpq_mat_struct *) mat;
    mats[1] = (fmpq_mat_struct *) mat2;

    _fmpq_mat_get_fmpz_mat_rowwise(nums, den, (const fmpq_mat_struct **) mats, 2);
}
