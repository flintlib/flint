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

#include "fmpz_mat.h"
#include "fmpq_vec.h"
#include "fmpq_mat.h"

int
fmpz_mat_is_reduced(const fmpz_mat_t A, double delta, double eta)
{
    slong i, j, d = A->r, n = A->c;
    fmpz_mat_t B;
    fmpq_mat_t Bq, Aq;
    mpq_t cx;
    fmpq_t c, p, q, r;

    if (d == 1)
        return 1;

    fmpz_mat_init(B, n, d);
    fmpq_mat_init(Bq, n, d);
    fmpq_mat_init(Aq, d, n);
    mpq_init(cx);
    fmpq_init(c);
    fmpq_init(p);
    fmpq_init(q);
    fmpq_init(r);

    fmpz_mat_transpose(B, A);
    fmpq_mat_set_fmpz_mat(Bq, B);
    fmpq_mat_gso(Bq, Bq);
    fmpq_mat_transpose(Aq, Bq);

    mpq_set_d(cx, delta - eta * eta);
    fmpq_set_mpq(c, cx);

    for (i = 0; i < d; i++)
    {
        _fmpq_vec_dot(p, Aq->rows[i], Aq->rows[i], n);
        for (j = 0; j < i; j++)
        {
            _fmpq_vec_dot(q, Aq->rows[j], Aq->rows[j], n);
            fmpq_pow_si(r, c, j - i);
            fmpq_mul(r, r, p);
            if (fmpq_cmp(q, r) > 0)
            {
                fmpz_mat_clear(B);
                fmpq_mat_clear(Bq);
                fmpq_mat_clear(Aq);
                mpq_clear(cx);
                fmpq_clear(c);
                fmpq_clear(p);
                fmpq_clear(q);
                fmpq_clear(r);
                return 0;
            }
        }
    }

    fmpz_mat_clear(B);
    fmpq_mat_clear(Bq);
    fmpq_mat_clear(Aq);
    mpq_clear(cx);
    fmpq_clear(c);
    fmpq_clear(p);
    fmpq_clear(q);
    fmpq_clear(r);

    return 1;
}
