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

int fmpq_mat_inv(fmpq_mat_t B, const fmpq_mat_t A)
{
    long n = A->r;

    if (n == 0)
    {
        return 1;
    }
    else if (n == 1)
    {
        if (fmpq_is_zero(fmpq_mat_entry(A, 0, 0)))
            return 0;
        fmpq_inv(fmpq_mat_entry(B, 0, 0), fmpq_mat_entry(A, 0, 0));
        return 1;
    }
    else if (n == 2)
    {
        fmpq_t d;
        int success;

        fmpq_init(d);

        fmpq_mul(d, fmpq_mat_entry(A, 0, 0), fmpq_mat_entry(A, 1, 1));
        fmpq_submul(d, fmpq_mat_entry(A, 0, 1), fmpq_mat_entry(A, 1, 0));
        success = !fmpq_is_zero(d);

        if (success)
        {
            fmpq_t t00, t01, t10, t11;
            fmpq_inv(d, d);

            fmpq_init(t00);
            fmpq_init(t01);
            fmpq_init(t10);
            fmpq_init(t11);

            fmpq_mul(t00, fmpq_mat_entry(A, 1, 1), d);
            fmpq_mul(t01, fmpq_mat_entry(A, 0, 1), d);
            fmpq_mul(t10, fmpq_mat_entry(A, 1, 0), d);
            fmpq_mul(t11, fmpq_mat_entry(A, 0, 0), d);

            fmpq_set(fmpq_mat_entry(B, 0, 0), t00);
            fmpq_neg(fmpq_mat_entry(B, 0, 1), t01);
            fmpq_neg(fmpq_mat_entry(B, 1, 0), t10);
            fmpq_set(fmpq_mat_entry(B, 1, 1), t11);

            fmpq_clear(t00);
            fmpq_clear(t01);
            fmpq_clear(t10);
            fmpq_clear(t11);
        }

        fmpq_clear(d);
        return success;
    }
    else
    {
        fmpz_mat_t Aclear, Bclear, I;
        fmpz * den;
        long i;
        int success;

        fmpz_mat_init(Aclear, n, n);
        fmpz_mat_init(Bclear, n, n);
        fmpz_mat_init(I, n, n);
        den = _fmpz_vec_init(n);

        fmpq_mat_get_fmpz_mat_rowwise(Aclear, den, A);
        for (i = 0; i < n; i++)
            fmpz_set(fmpz_mat_entry(I, i, i), den + i);

        success = fmpz_mat_solve(Bclear, den, Aclear, I);
        if (success)
            fmpq_mat_set_fmpz_mat_div_fmpz(B, Bclear, den);

        fmpz_mat_clear(Aclear);
        fmpz_mat_clear(Bclear);
        fmpz_mat_clear(I);
        _fmpz_vec_clear(den, A->r);

        return success;
    }
}
