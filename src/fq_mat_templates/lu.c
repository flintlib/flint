/*
    Copyright (C) 2011 Fredrik Johansson
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

slong
TEMPLATE(T, mat_lu) (slong * P, TEMPLATE(T, mat_t) A, int rank_check,
                     const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    slong rank;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_mat_lu(&rank, P, (gr_mat_struct *) A, (const gr_mat_struct *) A, rank_check, gr_ctx));
    return rank;
}


#endif
