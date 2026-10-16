/*
    Copyright (C) 2015 Elena Sergeicheva

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

void
TEMPLATE(T, mat_concat_vertical) (TEMPLATE(T, mat_t) res,
		                            const TEMPLATE(T, mat_t) mat1,
		                            const TEMPLATE(T, mat_t) mat2,
		                            const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_mat_concat_vertical((gr_mat_struct *) res,
        (const gr_mat_struct *) mat1, (const gr_mat_struct *) mat2, gr_ctx));
}

#endif
