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

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "perm.h"

#define E nmod_poly_mat_entry

int
nmod_poly_mat_inv(nmod_poly_mat_t Ainv, nmod_poly_t den,
                    const nmod_poly_mat_t A)
{
    len_t n = nmod_poly_mat_nrows(A);

    if (n == 0)
    {
        nmod_poly_one(den);
        return 1;
    }
    else if (n == 1)
    {
        nmod_poly_set(den, E(A, 0, 0));
        nmod_poly_one(E(Ainv, 0, 0));
        return !nmod_poly_is_zero(den);
    }
    else if (n == 2)
    {
        nmod_poly_mat_det(den, A);

        if (nmod_poly_is_zero(den))
        {
            return 0;
        }
        else if (Ainv == A)
        {
            nmod_poly_swap(E(A, 0, 0), E(A, 1, 1));
            nmod_poly_neg(E(A, 0, 1), E(A, 0, 1));
            nmod_poly_neg(E(A, 1, 0), E(A, 1, 0));
            return 1;
        }
        else
        {
            nmod_poly_set(E(Ainv, 0, 0), E(A, 1, 1));
            nmod_poly_set(E(Ainv, 1, 1), E(A, 0, 0));
            nmod_poly_neg(E(Ainv, 0, 1), E(A, 0, 1));
            nmod_poly_neg(E(Ainv, 1, 0), E(A, 1, 0));
            return 1;
        }
    }
    else
    {
        nmod_poly_mat_t LU, I;
        len_t * perm;
        int result;

        perm = _perm_init(n);
        nmod_poly_mat_init_set(LU, A);
        result = (nmod_poly_mat_fflu(LU, den, perm, LU, 1) == n);

        if (result)
        {
            nmod_poly_mat_init(I, n, n, nmod_poly_mat_modulus(A));
            nmod_poly_mat_one(I);
            nmod_poly_mat_solve_fflu_precomp(Ainv, perm, LU, I);
            nmod_poly_mat_clear(I);
        }
        else
            nmod_poly_zero(den);

        if (_perm_parity(perm, n))
        {
            nmod_poly_mat_neg(Ainv, Ainv);
            nmod_poly_neg(den, den);
        }

        _perm_clear(perm);
        nmod_poly_mat_clear(LU);
        return result;
    }
}
