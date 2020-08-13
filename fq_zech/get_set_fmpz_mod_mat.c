/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

void fq_zech_get_nmod_mat(nmod_mat_t col,
                          const fq_zech_t a,
                          const fq_zech_ctx_t ctx)
{
    fq_nmod_t tmp;
    fq_nmod_init(tmp, ctx->fq_nmod_ctx);
    fq_zech_get_fq_nmod(tmp, a, ctx);
    fq_nmod_get_nmod_mat(col, tmp, ctx->fq_nmod_ctx);
    fq_nmod_clear(tmp, ctx->fq_nmod_ctx);
}

void fq_zech_set_nmod_mat(fq_zech_t a,
                          const nmod_mat_t col,
                          const fq_zech_ctx_t ctx)
{
    fq_nmod_t tmp;
    fq_nmod_init(tmp, ctx->fq_nmod_ctx);
    fq_nmod_set_nmod_mat(tmp, col, ctx->fq_nmod_ctx);
    fq_zech_set_fq_nmod(a, tmp, ctx);
    fq_nmod_clear(tmp, ctx->fq_nmod_ctx);
}
