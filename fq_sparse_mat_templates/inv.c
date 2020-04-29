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

slong TEMPLATE(T, sparse_mat_inv) (TEMPLATE(T, sparse_mat_t) Mi, const TEMPLATE(T, sparse_mat_t) M, const TEMPLATE(T, ctx_t) ctx)
{
    slong rk;
    TEMPLATE(T, sparse_mat_t) I, MI;

    /* Create block matrix [M | I] */
    TEMPLATE(T, sparse_mat_init) (I, M->r, M->r, ctx);
    TEMPLATE(T, sparse_mat_one) (I, ctx);
    TEMPLATE(T, sparse_mat_init) (MI, M->r, M->r + M->c, ctx);
    TEMPLATE(T, sparse_mat_concat_horizontal) (MI, M, I, ctx);

    /* Run Gaussian elimination on first half */
    MI->c = M->c;
    rk = TEMPLATE(T, sparse_mat_rref) (MI, ctx);
    MI->c = M->c+M->r;
    TEMPLATE(T, sparse_mat_split_horizontal) (I, Mi, MI, M->c, ctx);
    TEMPLATE(T, sparse_mat_clear) (I, ctx);
    TEMPLATE(T, sparse_mat_clear) (MI, ctx); 
    return rk;
}

#endif
