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

#include "fmpq_vec.h"
#include "fmpq_mat.h"

void
fmpq_mat_gso(fmpq_mat_t B, fmpq_mat_t mu, const fmpq_mat_t A)
{
    slong i, j, k, d = A->r, n = A->c;
    fmpq_t num, den;

    if (B->r != A->r || B->c != A->c || mu->r != A->r || mu->c != A->r)
    {
        flint_printf("Exception (fmpq_mat_gso). Incompatible dimensions.\n");
        abort();
    }

    if (B == A)
    {
        fmpq_mat_t t;
        fmpq_mat_init(t, d, n);
        fmpq_mat_gso(t, mu, A);
        fmpq_mat_swap(B, t);
        fmpq_mat_clear(t);
        return;
    }

    if (n == 0)
    {
        fmpq_mat_one(mu);
        return;
    }

    fmpq_init(num);
    fmpq_init(den);

    fmpq_mat_one(mu);
    for (i = 0; i < d; i++)
    {
        for (j = 0; j < n; j++)
        {
            fmpq_set(fmpq_mat_entry(B, i, j), fmpq_mat_entry(A, i, j));
        }

        for (j = 0; j < i; j++)
        {
            _fmpq_vec_dot(num, A->rows[i], B->rows[j], n);
            _fmpq_vec_dot(den, B->rows[j], B->rows[j], n);

            if (fmpq_is_zero(den) == 0)
            {
                fmpq_div(fmpq_mat_entry(mu, i, j), num, den);

                for (k = 0; k < A->c; k++)
                {
                    fmpq_submul(fmpq_mat_entry(B, i, k),
                                fmpq_mat_entry(mu, i, j), fmpq_mat_entry(B, j,
                                                                         k));
                }
            }
        }
    }

    fmpq_clear(num);
    fmpq_clear(den);
}
