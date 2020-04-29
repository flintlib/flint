/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include <string.h>
#include "templates.h"

/* PAQ = LU, Ax = b => set b' = Pb, solve Ly = b', solve Ux' = y, set x = Qx' */
int TEMPLATE(T, sparse_mat_solve_lu) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *b, const TEMPLATE(T, ctx_t) ctx)
{
    int good = 1;
    slong rk, *P, *Q, i;
    TEMPLATE(T, t) cc;
    TEMPLATE(T, struct) *bp, *y, *xp;
    TEMPLATE(T, sparse_mat_t) L, U;
    if (_TEMPLATE(T, vec_is_zero) (b, M->r, ctx))
    {
        _TEMPLATE(T, vec_zero) (x, M->c, ctx);
        return 1;
    }

    P = flint_malloc(M->r * sizeof(*P));
    Q = flint_malloc(M->c * sizeof(*Q));
    TEMPLATE(T, init) (cc, ctx);
    bp = _TEMPLATE(T, vec_init) (M->r, ctx);
    xp = _TEMPLATE(T, vec_init) (M->c, ctx);
    TEMPLATE(T, sparse_mat_init) (L, M->r, M->c, ctx);
    TEMPLATE(T, sparse_mat_init) (U, M->r, M->c, ctx);

    rk = TEMPLATE(T, sparse_mat_lu) (P, Q, L, U, M, ctx);
    y = _TEMPLATE(T, vec_init) (rk, ctx);

    /* Solve Ly = b' = Pb */
    for (i = 0; i < M->r; ++i) TEMPLATE(T, set) (&bp[P[i]], &b[i], ctx);

    for (i = 0; i < rk; ++i)
    {
        TEMPLATE(T, sparse_vec_dot_dense) (cc, &L->rows[i], y, ctx);
        TEMPLATE(T, sub) (&y[i], &bp[i], cc, ctx);
    }
    for (i = rk; i < M->r; ++i)
    {
        TEMPLATE(T, sparse_vec_dot_dense) (cc, &L->rows[i], y, ctx);
        if (!TEMPLATE(T, equal) (&bp[i], cc, ctx)) {good = 0; break;}
    }

    if (good) 
    {
        /* Find a solution for Ux' = y */
        for (i = rk-1; i >= 0; --i) 
        {
            TEMPLATE(T, sparse_vec_dot_dense) (&xp[i], &U->rows[i], xp, ctx);
            TEMPLATE(T, sub) (&xp[i], &y[i], &xp[i], ctx);
            TEMPLATE(T, div) (&xp[i], &xp[i], U->rows[i].entries[0].val, ctx);
        }
        for (i = 0; i < M->c; ++i) TEMPLATE(T, set) (&x[i], &xp[Q[i]], ctx);
    }
    TEMPLATE(T, sparse_mat_clear) (L, ctx);
    TEMPLATE(T, sparse_mat_clear) (U, ctx);
    _TEMPLATE(T, vec_clear) (xp, M->c, ctx);
    _TEMPLATE(T, vec_clear) (y, rk, ctx);
    _TEMPLATE(T, vec_clear) (bp, M->r, ctx);
    return good;
}

#endif
