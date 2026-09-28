/*
    Copyright (C) 2019 William Hart
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "fq_nmod_vec.h"
#include "gr_vec.h"

void _fq_nmod_vec_dot(fq_nmod_t res, const fq_nmod_struct * vec1,
         const fq_nmod_struct * vec2, slong len2, const fq_nmod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fq_nmod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_vec_dot(res, NULL, 0, vec1, vec2, len2, gr_ctx));
}
