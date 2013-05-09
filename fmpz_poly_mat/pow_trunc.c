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

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_pow_trunc(fmpz_poly_mat_t B, const fmpz_poly_mat_t A, ulong exp,
                            len_t len)
{
    len_t d = fmpz_poly_mat_nrows(A);

    if (len < 1)
    {
        fmpz_poly_mat_zero(B);
    }
    else if (exp == 0 || d == 0)
    {
        fmpz_poly_mat_one(B);
    }
    else if (exp == 1)
    {
        fmpz_poly_mat_set(B, A);
        fmpz_poly_mat_truncate(B, len);
    }
    else if (exp == 2)
    {
        fmpz_poly_mat_sqrlow(B, A, len);
    }
    else if (d == 1)
    {
        fmpz_poly_pow_trunc(fmpz_poly_mat_entry(B, 0, 0),
                        fmpz_poly_mat_entry(A, 0, 0), exp, len);
    }
    else
    {
        fmpz_poly_mat_t T, U;
        len_t i;

        fmpz_poly_mat_init_set(T, A);
        fmpz_poly_mat_truncate(T, len);
        fmpz_poly_mat_init(U, d, d);

        for (i = ((len_t) FLINT_BIT_COUNT(exp)) - 2; i >= 0; i--)
        {
            fmpz_poly_mat_sqrlow(U, T, len);

            if (exp & (1L << i))
                fmpz_poly_mat_mullow(T, U, A, len);
            else
                fmpz_poly_mat_swap(T, U);
        }

        fmpz_poly_mat_swap(B, T);
        fmpz_poly_mat_clear(T);
        fmpz_poly_mat_clear(U);
    }
}
