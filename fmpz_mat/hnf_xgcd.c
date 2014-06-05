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

void fmpz_mat_hnf_xgcd(fmpz_mat_t H, const fmpz_mat_t A)
{
    slong j, j2, i, k, l;

    fmpz_mat_set(H, A);
    for (j = 0, k = 0, l = (A->c - A->r)*(A->c > A->r); A->c - j != l; j++, k++)
    {
        for (i = k + 1; i != A->r; i++)
        {
            fmpz_t id, kd, b, u, v, d;
            if (fmpz_is_zero(fmpz_mat_entry(H, i, j)))
                continue;
            fmpz_init(id);
            fmpz_init(kd);
            fmpz_init(b);
            fmpz_init(u);
            fmpz_init(v);
            fmpz_init(d);
            fmpz_xgcd(d, u, v, fmpz_mat_entry(H, k, j), fmpz_mat_entry(H, i, j));
            fmpz_divexact(id, fmpz_mat_entry(H, i, j), d);
            fmpz_divexact(kd, fmpz_mat_entry(H, k, j), d);
            for (j2 = 0; j2 < A->c; j2++)
            {
                fmpz_mul(b, u, fmpz_mat_entry(H, k, j2));
                fmpz_addmul(b, v, fmpz_mat_entry(H, i, j2));
                fmpz_mul(fmpz_mat_entry(H, i, j2), kd, fmpz_mat_entry(H, i, j2));
                fmpz_submul(fmpz_mat_entry(H, i, j2), id, fmpz_mat_entry(H, k, j2));
                fmpz_set(fmpz_mat_entry(H, k, j2), b);
            }
            fmpz_clear(id);
            fmpz_clear(kd);
            fmpz_clear(b);
            fmpz_clear(u);
            fmpz_clear(v);
            fmpz_clear(d);
        }
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
            /* reduce higher entries of column j with row k */
            for (i = k - 1; i >= 0; i--)
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
}
