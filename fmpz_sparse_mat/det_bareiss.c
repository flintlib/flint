/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"
#include "perm.h"

void
fmpz_sparse_mat_det_bareiss(fmpz_t det, const fmpz_sparse_mat_t M)
{
    slong i, rk, *P, *Q;
    fmpz *D;
    fmpz_sparse_mat_t L, U;
    if (M->r != M->c) {fmpz_zero(det); return;}
    fmpz_one(det);
    if (M->r == UWORD(0)) return; 
    
    P = flint_malloc(M->r*sizeof(*P));
    Q = flint_malloc(M->c*sizeof(*P));
    D = _fmpz_vec_init(M->r);
    fmpz_sparse_mat_init(L, M->r, M->c);
    fmpz_sparse_mat_init(U, M->r, M->c);
    rk = fmpz_sparse_mat_fflu(D, P, Q, L, U, M);
    if (rk != M->r) fmpz_zero(det);
    else
    {
        for (i = 0; i < M->r; ++i)
            fmpz_mul(det, det, U->rows[i].entries[0].val);
        for (i = 0; i < M->r; ++i)
            fmpz_divexact(det, det, &D[i]);
        if (_perm_parity(P, M->r) ^ _perm_parity(Q, M->c))
            fmpz_neg(det, det);
    }
    
    flint_free(P);
    flint_free(Q);
    _fmpz_vec_clear(D, M->r);
    fmpz_sparse_mat_clear(L);
    fmpz_sparse_mat_clear(U);
}
