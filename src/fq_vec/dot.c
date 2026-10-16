/*
    Copyright (C) 2019 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"
#include "fq_vec.h"
#include "gr_vec.h"

void _fq_vec_dot(fq_t res, const fq_struct * vec1,
         const fq_struct * vec2, slong len2, const fq_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fq_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_vec_dot(res, NULL, 0, vec1, vec2, len2, gr_ctx));
}
