/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "perm.h"

int
nmod_poly_mat_solve_fflu(nmod_poly_mat_t X, nmod_poly_t den,
                    const nmod_poly_mat_t A, const nmod_poly_mat_t B)
{
    nmod_poly_mat_t LU;
    slong dim, *perm;
    int result;

    if (nmod_poly_mat_is_empty(B))
    {
        nmod_poly_one(den);
        return 1;
    }

    dim = nmod_poly_mat_nrows(A);
    perm = _perm_init(dim);
    nmod_poly_mat_init_set(LU, A);
    result = (nmod_poly_mat_fflu(LU, den, perm, LU, 1) == dim);

    if (result)
    {
        nmod_poly_mat_solve_fflu_precomp(X, perm, LU, B);

        if (_perm_parity(perm, dim))
        {
            nmod_poly_neg(den, den);

	    nmod_poly_mat_neg(X, X);
        }
    } else
        nmod_poly_zero(den);

    _perm_clear(perm);
    nmod_poly_mat_clear(LU);
    return result;
}
