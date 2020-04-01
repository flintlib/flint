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
void nmod_sparse_mat_inv(nmod_sparse_mat_t Ai, const nmod_sparse_mat_t A)
{
    nmod_sparse_mat_t I, AI, window;
    nmod_sparse_vec_struct *row;
    nmod_sparse_entry_struct *le, *re;

    /* Create block matrix [A | I] */
    nmod_sparse_mat_init(I, A->r, A->r, A->mod);
    nmod_sparse_mat_one(I);
    nmod_sparse_mat_init(AI, A->r, A->c, A->mod);
    nmod_sparse_mat_concat_horizontal(AI, A, I);
    nmod_sparse_mat_clear(I);

    /* Run Gaussian elimination on first half */
    AI->c = A->c;
    nmod_sparse_mat_rref(AI);
    AI->c = A->c+A->r;
    nmod_sparse_mat_window_init(window, AI, 0, A->c, A->r, A->c+A->r);
    nmod_sparse_mat_set(Ai, window);

    /* TODO: ?? not invertible */
    nmod_sparse_mat_window_clear(window);
    nmod_sparse_mat_clear(AI);
}
