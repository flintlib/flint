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
    slong j, j2, i, k, l;

    fmpz_mat_set(H, A);

    for (j = 0, k = 0, l = (A->c - A->r)*(A->c > A->r); A->c - j != l; j++, k++)
    {
        int col_finished = 1;
        for (i = k + 1; (i < A->r) && col_finished; i++)
        {
            col_finished = fmpz_is_zero(fmpz_mat_entry(H, i, j));
        }
        if (col_finished)
        {
            if (fmpz_sgn(fmpz_mat_entry(H, k, j)) < 0)
            {
                for (j2 = 0; j2 < A->c; j2++)
                {
                    fmpz_neg(fmpz_mat_entry(H, k, j2),
                            fmpz_mat_entry(H, k, j2));
                }
            }
            if (fmpz_is_zero(fmpz_mat_entry(H, k, j)))
            {
                k--;
                if (l > 0)
                    l--;
            }
            else
            {
                /* reduce first entries of column j with row k */
                for (i = 0; i < k; i++)
                {
                    fmpz_t q;
                    fmpz_init(q);
                    fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                            fmpz_mat_entry(H, k, j));
                    for (j2 = 0; j2 < A->c; j2++)
                    {
                        fmpz_submul(fmpz_mat_entry(H, i, j2), q,
                                fmpz_mat_entry(H, k, j2));
                    }
                    fmpz_clear(q);
                }
            }
        }
        else
        {
            slong i0 = 0;
            fmpz_t min;

            fmpz_init(min);
            /* locate non-zero entry in column j below k with lowest absolute value */
            for (i = k + 1; i < A->r; i++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(H, i, j)))
                    continue;
                if (fmpz_is_zero(min) ||
                        fmpz_cmpabs(fmpz_mat_entry(H, i, j), min) < 0)
                {
                    i0 = i;
                    fmpz_abs(min, fmpz_mat_entry(H, i, j));
                }
            }
            fmpz_clear(min);
            /* move the row found to row k */
            if (i0 > k)
            {
                for (j2 = 0; j2 < A->c; j2++)
                {
                    fmpz_swap(fmpz_mat_entry(H, i0, j2),
                            fmpz_mat_entry(H, k, j2));
                }
            }
            if (fmpz_sgn(fmpz_mat_entry(H, k, j)) < 0)
            {
                for (j2 = 0; j2 < A->c; j2++)
                {
                    fmpz_neg(fmpz_mat_entry(H, k, j2),
                            fmpz_mat_entry(H, k, j2));
                }
            }
            /* reduce lower entries of column j with row k */
            for (i = k + 1; i < A->r; i++)
            {
                fmpz_t q;
                fmpz_init(q);
                fmpz_fdiv_q(q, fmpz_mat_entry(H, i, j),
                        fmpz_mat_entry(H, k, j));
                for (j2 = 0; j2 < A->c; j2++)
                {
                    fmpz_submul(fmpz_mat_entry(H, i, j2), q,
                            fmpz_mat_entry(H, k, j2));
                }
                fmpz_clear(q);
            }
            /* don't move to the next column yet */
            j--;
            k--;
        }
    }
}
