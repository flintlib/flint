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

int TEMPLATE(T, sparse_mat_solve_lanczos) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, struct) *b, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx)
{
    slong j, ret;

    /* We assume that M is not symmetric, and work with A = M^t M */
    TEMPLATE(T, t) cinv, cc, AvtAv, vMtb, alpha, beta;
    TEMPLATE(T, struct) *v[2], *Mv, *Av, *Mtb, vtAv[2];
    TEMPLATE(T, sparse_mat_t) Mt;

    if (_TEMPLATE(T, vec_is_zero) (b, M->r, ctx))
    {
        _TEMPLATE(T, vec_zero) (x, M->c, ctx);
        return 1;
    }

    TEMPLATE(T, init) (cinv, ctx);
    TEMPLATE(T, init) (cc, ctx);
    TEMPLATE(T, init) (AvtAv, ctx);
    TEMPLATE(T, init) (vMtb, ctx);
    TEMPLATE(T, init) (alpha, ctx);
    TEMPLATE(T, init) (beta, ctx);
    TEMPLATE(T, init) (&vtAv[0], ctx);
    TEMPLATE(T, init) (&vtAv[1], ctx);
    v[0] = _TEMPLATE(T, vec_init) (M->c, ctx);
    v[1] = _TEMPLATE(T, vec_init) (M->c, ctx);
    Mv = _TEMPLATE(T, vec_init) (M->r, ctx);
    Av = _TEMPLATE(T, vec_init) (M->c, ctx);
    Mtb = _TEMPLATE(T, vec_init) (M->c, ctx);
    TEMPLATE(T, sparse_mat_init) (Mt, M->c, M->r, ctx);

    /* Construct transpose */
    TEMPLATE(T, sparse_mat_transpose) (Mt, M, ctx);


    /* Set starting data */
    _TEMPLATE(T, vec_randtest) (v[0], state, M->c, ctx);
    _TEMPLATE(T, vec_zero) (v[1], M->c, ctx); 
    TEMPLATE(T, one) (&vtAv[1], ctx);
    _TEMPLATE(T, vec_zero) (x, M->c, ctx);
    TEMPLATE(T, sparse_mat_mul_vec) (Mtb, Mt, b, ctx);
    for (j = 0; ; j = 1-j)
    {
        /* Compute M^T M v_j and check if it is orthogonal to v_j */
        TEMPLATE(T, sparse_mat_mul_vec) (Mv, M, v[j], ctx);
        TEMPLATE(T, sparse_mat_mul_vec) (Av, Mt, Mv, ctx);
        _TEMPLATE(T, vec_dot) (&vtAv[j], v[j], Av, M->c, ctx);
        if (TEMPLATE(T, is_zero) (&vtAv[j], ctx)) break; /* Can't make any more progress */

        /* Update putative solution by <v_j, M^T b>/delta_j * v_j */
        _TEMPLATE(T, vec_dot) (cc, v[j], Mtb, M->c, ctx);
        TEMPLATE(T, div) (cc, cc, &vtAv[j], ctx);
        _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T)) (x, v[j], M->c, cc, ctx);

        /* v_{j+1} = AtAv - alpha*v_j - beta*v_{j-1}, where */
        /*    alpha = <Av_j, Av_j>/delta_j, and */
        /*    beta = delta_j/delta_{j-1} */
        TEMPLATE(T, div) (cc, &vtAv[j], &vtAv[1-j], ctx);
        TEMPLATE(T, neg) (cc, cc, ctx);
        _TEMPLATE(T, TEMPLATE(vec_scalar_mul, T))(v[1-j], v[1-j], M->c, cc, ctx);
        _TEMPLATE(T, vec_dot) (cc, Av, Av, M->c, ctx);
        TEMPLATE(T, div) (cc, cc, &vtAv[j], ctx);
        TEMPLATE(T, neg) (cc, cc, ctx);
        _TEMPLATE(T, TEMPLATE(vec_scalar_addmul, T)) (v[1-j], v[j], M->c, cc, ctx);
        _TEMPLATE(T, vec_add) (v[1-j], v[1-j], Av, M->c, ctx);
    }
    /* Check result */
    TEMPLATE(T, sparse_mat_mul_vec) (Mv, M, x, ctx);
    TEMPLATE(T, sparse_mat_mul_vec) (Av, Mt, Mv, ctx);
    ret = _TEMPLATE(T, vec_equal) (Av, Mtb, M->c, ctx);

    /* Clear auxiliary vectors and transpose */
    TEMPLATE(T, clear) (cc, ctx);
    TEMPLATE(T, clear) (cinv, ctx);
    TEMPLATE(T, clear) (AvtAv, ctx);
    TEMPLATE(T, clear) (vMtb, ctx);
    TEMPLATE(T, clear) (alpha, ctx);
    TEMPLATE(T, clear) (beta, ctx);
    TEMPLATE(T, clear) (&vtAv[0], ctx);
    TEMPLATE(T, clear) (&vtAv[1], ctx);
    _TEMPLATE(T, vec_clear) (v[0], M->c, ctx);
    _TEMPLATE(T, vec_clear) (v[1], M->c, ctx);
    _TEMPLATE(T, vec_clear) (Mv, M->r, ctx);
    _TEMPLATE(T, vec_clear) (Av, M->c, ctx);
    _TEMPLATE(T, vec_clear) (Mtb, M->c, ctx);
    TEMPLATE(T, sparse_mat_clear) (Mt, ctx);
    return ret;
}

int TEMPLATE(T, sparse_mat_nullvector_lanczos) (TEMPLATE(T, struct) *x, const TEMPLATE(T, sparse_mat_t) M, flint_rand_t state, const TEMPLATE(T, ctx_t) ctx) 
{
    int ret = 1;
    TEMPLATE(T, struct) *x2, *b;
    x2 = _TEMPLATE(T, vec_init) (M->c, ctx);
    b = _TEMPLATE(T, vec_init) (M->r, ctx);

    _TEMPLATE(T, vec_randtest) (x, state, M->c, ctx);
    TEMPLATE(T, sparse_mat_mul_vec) (b, M, x, ctx);
    if (TEMPLATE(T, sparse_mat_solve_lanczos) (x2, M, b, state, ctx) == 0) ret = 0; /* Lanczos failed */
    if (ret)
    {
        _TEMPLATE(T, vec_sub) (x, x, x2, M->c, ctx);
        TEMPLATE(T, sparse_mat_mul_vec) (b, M, x, ctx);
        ret = !_TEMPLATE(T, vec_is_zero) (x, M->c, ctx) && _TEMPLATE(T, vec_is_zero) (b, M->r, ctx);
    }
    _TEMPLATE(T, vec_clear) (x2, M->c, ctx);
    _TEMPLATE(T, vec_clear) (b, M->r, ctx);
    return ret;
}

#endif
