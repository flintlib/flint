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


int nmod_sparse_mat_solve_lanczos(mp_ptr x, const nmod_sparse_mat_t M, mp_srcptr b, flint_rand_t state)
{
    slong j, ret;
    const slong nlimbs = _nmod_vec_dot_bound_limbs(M->c, M->mod);

    /* We assume that M is not symmetric, and work with A = M^t M */
    nmod_sparse_mat_t Mt;
    mp_ptr v[2], Mv, Av, Mtb;
    mp_limb_t vtAv[2], AvtAv, vMtb;

    _nmod_vec_zero(x, M->c);
    if (_nmod_vec_is_zero(b, M->c)) return 1;

    /* Construct transpose */
    nmod_sparse_mat_init(Mt, M->c, M->r, M->mod);
    nmod_sparse_mat_transpose(Mt, M);

    /* Construct auxiliary vectors */
    /* Rather than storing the whole sequence of values v_j, we alternate between two vectors */
    v[0] = _nmod_vec_init(M->c);
    v[1] = _nmod_vec_init(M->c);
    Mv = _nmod_vec_init(M->r);
    Av = _nmod_vec_init(M->c);
    Mtb = _nmod_vec_init(M->c);
    nmod_sparse_mat_mul_vec(Mtb, Mt, b);

    /* Make 0th vector random (and -1st vector trivial) */
    /*_nmod_vec_set(v[0], Mtb, M->c);
    for (j = 0; j < M->c; ++j) v[0][j] = n_randint(state, M->mod.n); */
    _nmod_vec_randtest(v[0], state, M->c, M->mod);
    _nmod_vec_zero(v[1], M->c); vtAv[1] = 1;  
    for (j = 0; ; j = 1-j)
    {
        /* Compute M^T M v_j and check if it is orthogonal to v_j */
        nmod_sparse_mat_mul_vec(Mv, M, v[j]);
        nmod_sparse_mat_mul_vec(Av, Mt, Mv);
        vtAv[j] = _nmod_vec_dot(v[j], Av, M->c, M->mod, nlimbs);
        if (vtAv[j] == UWORD(0)) break; /* Can't make any more progress */

        /* Update putative solution by <v_j, M^T b>/delta_j * v_j */
        vMtb = nmod_div(_nmod_vec_dot(v[j], Mtb, M->c, M->mod, nlimbs), vtAv[j], M->mod);
        _nmod_vec_scalar_addmul_nmod(x, v[j], M->c, vMtb, M->mod);

        /* v_{j+1} = MtMv - alpha*v_j - beta*v_{j-1}, where */
        /*    alpha = <Mv_j, Mv_j>/delta_j, and */
        /*    beta = delta_j/delta_{j-1} */
        AvtAv = _nmod_vec_dot(Av, Av, M->c, M->mod, nlimbs);
        _nmod_vec_scalar_mul_nmod(v[1-j], v[1-j], M->c, nmod_neg(nmod_div(vtAv[j], vtAv[1-j], M->mod), M->mod), M->mod);
        _nmod_vec_scalar_addmul_nmod(v[1-j], v[j], M->c, nmod_neg(nmod_div(AvtAv, vtAv[j], M->mod), M->mod), M->mod);
        _nmod_vec_add(v[1-j], v[1-j], Av, M->c, M->mod);
    }
    /* Check result */
    nmod_sparse_mat_mul_vec(Mv, M, x);
    nmod_sparse_mat_mul_vec(Av, Mt, Mv);
    ret = _nmod_vec_equal(Av, Mtb, M->c);

    /* Clear auxiliary vectors and transpose */
    _nmod_vec_clear(v[0]);
    _nmod_vec_clear(v[1]);
    _nmod_vec_clear(Mv);
    _nmod_vec_clear(Av);
    _nmod_vec_clear(Mtb);
    nmod_sparse_mat_clear(Mt);
    return ret;
}

int nmod_sparse_mat_nullvector_lanczos(mp_ptr x, const nmod_sparse_mat_t M, flint_rand_t state) 
{
    int ret = 1;
    mp_ptr x2, b;
    x2 = _nmod_vec_init(M->c);
    b = _nmod_vec_init(M->r);

    _nmod_vec_randtest(x, state, M->c, M->mod);
    nmod_sparse_mat_mul_vec(b, M, x);
    if (nmod_sparse_mat_solve_lanczos(x2, M, b, state) == 0) ret = 0; /* Lanczos failed */
    if (ret)
    {
        _nmod_vec_sub(x, x, x2, M->c, M->mod);
        nmod_sparse_mat_mul_vec(b, M, x);
        ret = !_nmod_vec_is_zero(x, M->c) && _nmod_vec_is_zero(b, M->r);
    }
    _nmod_vec_clear(x2);
    _nmod_vec_clear(b);
    return ret;
}
