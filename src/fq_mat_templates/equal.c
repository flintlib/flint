/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "gr.h"
#include "gr_mat.h"
#include "templates.h"

int TEMPLATE(T, mat_equal) (const TEMPLATE(T, mat_t) mat1,
                            const TEMPLATE(T, mat_t) mat2,
                            const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    return gr_mat_equal((const gr_mat_struct *) mat1,
                        (const gr_mat_struct *) mat2, gr_ctx) == T_TRUE;
}

#endif
