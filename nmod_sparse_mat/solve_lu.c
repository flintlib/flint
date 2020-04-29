/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_vec.h"
#include "nmod_sparse_mat.h"

/* PAQ = LU, Ax = b => set b' = Pb, solve Ly = b', solve Ux' = y, set x = Qx' */
int nmod_sparse_mat_solve_lu(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b)
{
    int good = 1;
    slong rk, *P, *Q, i;
    nmod_sparse_mat_t L, U;
    mp_ptr bp, y, xp;
    if (_nmod_vec_is_zero(b, M->c))
    {
        _nmod_vec_zero(x, M->c);
        return 1;
    }

    P = flint_malloc(M->r * sizeof(*P));
    Q = flint_malloc(M->c * sizeof(*Q));
    nmod_sparse_mat_init(L, M->r, M->c, M->mod);
    nmod_sparse_mat_init(U, M->r, M->c, M->mod);
    rk = nmod_sparse_mat_lu(P, Q, L, U, M);
    L->c = U->r = rk;

    /* Solve Ly = b' = Pb */
    bp = flint_malloc(M->r * sizeof(*bp));
    y = flint_calloc(rk, sizeof(*y));
    for (i = 0; i < M->r; ++i) bp[P[i]] = b[i];

    flint_free(P);
    for (i = 0; i < rk; ++i)
        y[i] = nmod_sub(bp[i], nmod_sparse_vec_dot_dense(&L->rows[i], y, M->mod), M->mod);
    for (i = rk; i < M->r; ++i)
        if (bp[i] != nmod_sparse_vec_dot_dense(&L->rows[i], y, M->mod)) {good = 0; break;}
    nmod_sparse_mat_mul_vec(bp, L, y);

    flint_free(bp);
    nmod_sparse_mat_clear(L);

    if (good) 
    {
        /* Find a solution for Ux' = y */
        xp = flint_calloc(M->c, sizeof(*xp));
        for (i = rk-1; i >= 0; --i) 
            xp[i] = nmod_div(nmod_sub(y[i], nmod_sparse_vec_dot_dense(&U->rows[i], xp, M->mod), M->mod), U->rows[i].entries[0].val, M->mod);
        nmod_sparse_mat_mul_vec(y, U, xp);

        for (i = 0; i < M->c; ++i) x[i] = xp[Q[i]];
        flint_free(xp);
    }
    flint_free(Q);
    flint_free(y);
    U->r = M->r;
    nmod_sparse_mat_clear(U);
    return good;
}
