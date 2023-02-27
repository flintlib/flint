/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_embed.h"

void fq_zech_embed_composition_matrix_sub(nmod_mat_t matrix,
                                    const fq_zech_t gen,
                                    const fq_zech_ctx_t ctx,
                                    slong trunc)
{
    fq_nmod_t gen_nmod;
    fq_nmod_ctx_struct *modulus = ctx->fq_nmod_ctx;
    fq_nmod_init(gen_nmod, modulus);
    fq_zech_get_fq_nmod(gen_nmod, gen, ctx);
    fq_nmod_embed_composition_matrix_sub(matrix, gen_nmod, modulus, trunc);
    fq_nmod_clear(gen_nmod, modulus);
}
