/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#include "templates.h"
#include "gr.h"
#include "gr_mat.h"

void
TEMPLATE(T, mat_minpoly) (TEMPLATE(T, poly_t) p,
                      const TEMPLATE(T, mat_t) X, const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;

    if (X->r != X->c)
        flint_throw(FLINT_ERROR, "Exception (fq_mat_minpoly).  Non-square matrix.\n");

    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_mat_minpoly_field((gr_poly_struct *) p, (const gr_mat_struct *) X, gr_ctx));
}

#endif
