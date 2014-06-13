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

int
fmpz_mat_is_reduced_gram(const fmpz_mat_t A, double delta, double eta)
{
    int result = 1;
    slong i, d = A->r;
    mpq_t cx;
    fmpz_t det;
    fmpq_t c, q, r;

    if (d == 1)
        return 1;

    mpq_init(cx);
    fmpz_init(det);
    fmpq_init(c);
    fmpq_init(q);
    fmpq_init(r);

    fmpz_mat_det(det, A);
    fmpz_pow_ui(det, det, 2);

    mpq_set_d(cx, delta - eta * eta);
    fmpq_set_mpq(c, cx);

    fmpq_pow_si(r, c, d * (1 - d));
    fmpq_mul_fmpz(r, r, det);

    fmpq_one(q);
    for (i = 0; i < d; i++)
    {
        fmpq_mul_fmpz(q, q, fmpz_mat_entry(A, i, i));
    }
    fmpq_pow_si(q, q, 2);

    if (fmpq_cmp(q, r) > 0)
        result = 0;

    mpq_clear(cx);
    fmpz_clear(det);
    fmpq_clear(c);
    fmpq_clear(q);
    fmpq_clear(r);

    return result;
}
