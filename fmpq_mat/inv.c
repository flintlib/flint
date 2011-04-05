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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

void fmpq_mat_inv(fmpq_mat_t B, const fmpq_mat_t A)
{
    long n, i, j;

    fmpz_mat_t Aclear;
    fmpz_mat_t Bclear;
    fmpz_mat_t I;

    fmpz * den;

    n = A->r;

    if (n == 0)
    {
        return;
    }
    else if (n == 1)
    {
        fmpq_inv(fmpq_mat_entry(B, 0, 0), fmpq_mat_entry(A, 0, 0));
        return;
    }

    fmpz_mat_init(Aclear, n, n);
    fmpz_mat_init(Bclear, n, n);
    fmpz_mat_init(I, n, n);

    den = _fmpz_vec_init(n);

    fmpq_mat_get_fmpz_mat_rowwise(Aclear, den, A);

    for (i = 0; i < n; i++)
        fmpz_set(fmpz_mat_entry(I, i, i), den + i);

    fmpz_mat_solve_mat(Bclear, den, Aclear, I);

    for (i = 0; i < B->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            fmpz_set(fmpq_mat_entry_num(B, i, j), fmpz_mat_entry(Bclear, i, j));
            fmpz_set(fmpq_mat_entry_den(B, i, j), den);
            fmpq_canonicalise(fmpq_mat_entry(B, i, j));
        }
    }

    fmpz_mat_clear(Aclear);
    fmpz_mat_clear(Bclear);
    fmpz_mat_clear(I);

    _fmpz_vec_clear(den, A->r);
}
