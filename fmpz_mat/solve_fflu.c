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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "perm.h"

int
fmpz_mat_solve_fflu(fmpz_mat_t X, fmpz_t den,
                    const fmpz_mat_t A, const fmpz_mat_t B)
{
    fmpz_mat_t LU;
    len_t dim, *perm;
    int result;

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
    {
        fmpz_one(den);
        return 1;
    }

    dim = fmpz_mat_nrows(A);
    perm = _perm_init(dim);
    fmpz_mat_init_set(LU, A);
    result = (fmpz_mat_fflu(LU, den, perm, LU, 1) == dim);

    if (result)
        fmpz_mat_solve_fflu_precomp(X, perm, LU, B);
    else
        fmpz_zero(den);

    _perm_clear(perm);
    fmpz_mat_clear(LU);
    return result;
}
