/*
<<<<<<< HEAD
    Copyright (C) 2011 Fredrik Johansson
=======
    Copyright (C) 2010 Fredrik Johansson
>>>>>>> Initial code for sparse matrices mod limb size integers, just construction and arithmetic for starters

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"

void nmod_sparse_mat_window_init(nmod_sparse_mat_t window, const nmod_sparse_mat_t mat, slong r1, slong c1, slong r2, slong c2) {
    slong i;
    r2 = FLINT_MIN(r2, mat->r), r1 = FLINT_MIN(r1, r2);
    c2 = FLINT_MIN(c2, mat->c), c1 = FLINT_MIN(c1, c2);
    window->mod = mat->mod;
    window->r = r2-r1;
    window->c = c2-c1;
    window->c_off = c1;
    window->rows = flint_malloc(window->r*sizeof(*window->rows));
    for(i=0; i<window->r; ++i)
        nmod_sparse_vec_window_init(&window->rows[i], &mat->rows[i+r1], c1, c2);
}
