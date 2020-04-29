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
slong nmod_sparse_mat_inv(nmod_sparse_mat_t Mi, const nmod_sparse_mat_t M)
{
    slong rk;
    nmod_sparse_mat_t I, MI;

    /* Create block matrix [M | I] */
    nmod_sparse_mat_init(I, M->r, M->r, M->mod);
    nmod_sparse_mat_one(I);
    nmod_sparse_mat_init(MI, M->r, M->r + M->c, M->mod);
    nmod_sparse_mat_concat_horizontal(MI, M, I);

    /* Run Gaussian elimination on first half */
    MI->c = M->c;
    rk = nmod_sparse_mat_rref(MI);
    MI->c = M->c+M->r;
    nmod_sparse_mat_split_horizontal(I, Mi, MI, M->c);
    nmod_sparse_mat_clear(I);
    nmod_sparse_mat_clear(MI); 
    return rk;
}
