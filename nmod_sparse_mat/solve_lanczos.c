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
#include "nmod_sparse_mat.h"


int nmod_sparse_mat_solve_lanczos(mp_ptr x, const nmod_sparse_mat_t A, const mp_ptr b, flint_rand_t state) {
    _nmod_vec_zero(x, A->c);

    /* Construct transpose */
    nmod_sparse_mat_t At;
    nmod_sparse_mat_init(At, A->c, A->r, A->mod);
    nmod_sparse_mat_transpose(At, A);

    /* Construct auxiliary vectors */
    /* Rather than storing the whole sequence of values w_j, we alternate between two vectors */
    slong j, iter;
    const slong nlimbs = _nmod_vec_dot_bound_limbs(A->c, A->mod);
    mp_ptr w[2], Aw, AtAw, Atb;
    mp_limb_t delta[2];
    w[0] = _nmod_vec_init(A->c);
    w[1] = _nmod_vec_init(A->c);
    Aw = _nmod_vec_init(A->r);
    AtAw = _nmod_vec_init(A->c);
    Atb = _nmod_vec_init(A->c);
    nmod_sparse_mat_mul_vec(Atb, At, b);

    /* Make 0th vector random (and -1st vector trivial) */
    _nmod_vec_randtest(w[0], state, A->c, A->mod);
    _nmod_vec_zero(w[1], A->c); delta[1] = 1;  
    for(j=0; ; j=1-j) {
        /* Compute A^T A w_j and check if it is orthogonal to w_j */
        nmod_sparse_mat_mul_vec(Aw, A, w[j]);
        nmod_sparse_mat_mul_vec(AtAw, At, Aw);
        delta[j] = _nmod_vec_dot(w[j], AtAw, A->c, A->mod, nlimbs);
        if (delta[j]==UWORD(0)) break; // Can't make any more progress

        /* Update putative solution by <w_j, A^T b>/delta_j * w_j */
        const mp_limb_t wAtb = nmod_div(_nmod_vec_dot(w[j], Atb, A->c, A->mod, nlimbs), delta[j], A->mod);
        _nmod_vec_scalar_addmul_nmod(x, w[j], A->c, wAtb, A->mod);

        /* w_{j+1} = AtAw - alpha*w_j - beta*w_{j-1}, where */
        /*    alpha = <Aw_j, Aw_j>/delta_j, and */
        /*    beta = delta_j/delta_{j-1} */
        const mp_limb_t alpha = nmod_div(_nmod_vec_dot(AtAw, AtAw, A->c, A->mod, nlimbs), delta[j], A->mod);
        const mp_limb_t beta = nmod_div(delta[j], delta[1-j], A->mod);
        _nmod_vec_scalar_mul_nmod(w[1-j], w[1-j], A->c, nmod_neg(beta, A->mod), A->mod);
        _nmod_vec_scalar_addmul_nmod(w[1-j], w[j], A->c, nmod_neg(alpha, A->mod), A->mod);
        _nmod_vec_add(w[1-j], w[1-j], AtAw, A->c, A->mod);
    }
    /* Check result */
    nmod_sparse_mat_mul_vec(Aw, A, x);
    nmod_sparse_mat_mul_vec(AtAw, At, Aw);
    int ret = _nmod_vec_equal(AtAw, Atb, A->c);

    /* Clear auxiliary vectors and transpose */
    _nmod_vec_clear(w[0]);
    _nmod_vec_clear(w[1]);
    _nmod_vec_clear(Aw);
    _nmod_vec_clear(AtAw);
    _nmod_vec_clear(Atb);
    nmod_sparse_mat_clear(At);
    return ret;
}
