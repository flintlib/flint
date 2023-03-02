/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"
#include "perm.h"

int
fmpz_poly_mat_solve_fflu(fmpz_poly_mat_t X, fmpz_poly_t den,
                    const fmpz_poly_mat_t A, const fmpz_poly_mat_t B)
{
    fmpz_poly_mat_t LU;
    slong dim, *perm;
    int result;

    if (fmpz_poly_mat_is_empty(B))
    {
        fmpz_poly_one(den);
        return 1;
    }

    dim = fmpz_poly_mat_nrows(A);
    perm = _perm_init(dim);
    fmpz_poly_mat_init_set(LU, A);
    result = (fmpz_poly_mat_fflu(LU, den, perm, LU, 1) == dim);

    if (result)
    {
        fmpz_poly_mat_solve_fflu_precomp(X, perm, LU, B);

        if (_perm_parity(perm, dim))
	{
            fmpz_poly_neg(den, den);

	    fmpz_poly_mat_neg(X, X);
        } 
    } else
        fmpz_poly_zero(den);

    _perm_clear(perm);
    fmpz_poly_mat_clear(LU);
    return result;
}
