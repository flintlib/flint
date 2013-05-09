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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"
#include "perm.h"

#define E fmpz_poly_mat_entry

int
fmpz_poly_mat_inv(fmpz_poly_mat_t Ainv, fmpz_poly_t den,
                    const fmpz_poly_mat_t A)
{
    len_t n = fmpz_poly_mat_nrows(A);

    if (n == 0)
    {
        fmpz_poly_one(den);
        return 1;
    }
    else if (n == 1)
    {
        fmpz_poly_set(den, E(A, 0, 0));
        fmpz_poly_one(E(Ainv, 0, 0));
        return !fmpz_poly_is_zero(den);
    }
    else if (n == 2)
    {
        fmpz_poly_mat_det(den, A);
        if (fmpz_poly_is_zero(den))
        {
            return 0;
        }
        else if (Ainv == A)
        {
            fmpz_poly_swap(E(A, 0, 0), E(A, 1, 1));
            fmpz_poly_neg(E(A, 0, 1), E(A, 0, 1));
            fmpz_poly_neg(E(A, 1, 0), E(A, 1, 0));
            return 1;
        }
        else
        {
            fmpz_poly_set(E(Ainv, 0, 0), E(A, 1, 1));
            fmpz_poly_set(E(Ainv, 1, 1), E(A, 0, 0));
            fmpz_poly_neg(E(Ainv, 0, 1), E(A, 0, 1));
            fmpz_poly_neg(E(Ainv, 1, 0), E(A, 1, 0));
            return 1;
        }
    }
    else
    {
        fmpz_poly_mat_t LU, I;
        len_t * perm;
        int result;

        perm = _perm_init(n);
        fmpz_poly_mat_init_set(LU, A);
        result = (fmpz_poly_mat_fflu(LU, den, perm, LU, 1) == n);

        if (result)
        {
            fmpz_poly_mat_init(I, n, n);
            fmpz_poly_mat_one(I);
            fmpz_poly_mat_solve_fflu_precomp(Ainv, perm, LU, I);
            fmpz_poly_mat_clear(I);
        }
        else
            fmpz_poly_zero(den);

        if (_perm_parity(perm, n))
        {
            fmpz_poly_mat_neg(Ainv, Ainv);
            fmpz_poly_neg(den, den);
        }

        _perm_clear(perm);
        fmpz_poly_mat_clear(LU);
        return result;
    }
}
