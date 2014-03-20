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

void fmpz_mat_hnf_classical(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong i, i2, j, k, l;
    fmpz_mat_set(H, A);

    for (i = A->r - 1, k = A->c - 1, l = (i - k)*(i > k); i + 1 != l; i--, k--)
    {
        int row_finished = 1;
        for (j = 0; (j < k) && row_finished; j++)
        {
            row_finished = fmpz_is_zero(fmpz_mat_entry(H, i, j));
        }
        if (row_finished)
        {
            if (fmpz_sgn(fmpz_mat_entry(H, i, k)) < 0)
            {
                for (i2 = 0; i2 < A->r; i2++)
                {
                    fmpz_neg(fmpz_mat_entry(H, i2, k),
                            fmpz_mat_entry(H, i2, k));
                }
            }
            if (fmpz_is_zero(fmpz_mat_entry(H, i, k)))
            {
                k++;
            }
            else
            {
                /* reduce trailing entries of row i with column k */
                for (j = k + 1; j < A->c; j++)
                {
                    fmpz_t q;
                    fmpz_init(q);
                    fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                            fmpz_mat_entry(H, i, k));
                    for (i2 = 0; i2 < A->r; i2++)
                    {
                        fmpz_submul(fmpz_mat_entry(H, i2, j), q,
                                fmpz_mat_entry(H, i2, k));
                    }
                    fmpz_clear(q);
                }
            }
        }
        else
        {
            slong j0 = 0;
            fmpz_t min;
            fmpz_init(min);
            /* locate non-zero entry on row i with lowest absolute value */
            for (j = 0; j <= k; j++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(H, i, j)))
                    continue;
                if (fmpz_is_zero(min) || 
                        fmpz_cmpabs(fmpz_mat_entry(H, i, j), min) < 0)
                {
                    j0 = j;
                    fmpz_abs(min, fmpz_mat_entry(H, i, j));
                }
            }
            fmpz_clear(min);
            /* move the column found to column k */
            if (j0 < k)
            {
                for (i2 = 0; i2 < A->r; i2++)
                {
                    fmpz_swap(fmpz_mat_entry(H, i2, j0),
                            fmpz_mat_entry(H, i2, k));
                }
            }
            if (fmpz_sgn(fmpz_mat_entry(H, i, k)) < 0)
            {
                for (i2 = 0; i2 < A->r; i2++)
                {
                    fmpz_neg(fmpz_mat_entry(H, i2, k),
                            fmpz_mat_entry(H, i2, k));
                }
            }
            /* reduce leading entries of row i with column k */
            for (j = 0; j <= k - 1; j++)
            {
                fmpz_t q;
                fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                        fmpz_mat_entry(H, i, k));
                for (i2 = 0; i2 < A->r; i2++)
                {
                    fmpz_submul(fmpz_mat_entry(H, i2, j), q,
                            fmpz_mat_entry(H, i2, k));
                }
            }
            /* don't move to the next row yet */
            i++;
            k++;
        }
    }
}
