/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_embed.h"

void fq_nmod_embed_mul_matrix(nmod_mat_t matrix,
                        const fq_nmod_t gen,
                        const fq_nmod_ctx_t ctx) {
    slong i, j, len = fq_nmod_ctx_degree(ctx);
    const nmod_poly_struct *modulus = ctx->modulus;
    const nmod_t mod = modulus->mod;
    mp_limb_t lead;

    /* This is usually 1, unless the context is non-monic */
    lead = nmod_inv(modulus->coeffs[len], mod);

    for (i = 0; i < gen->length; i++)
        nmod_mat_entry(matrix, i, 0) =  gen->coeffs[i];
    for (i = gen->length; i < len; i++)
        nmod_mat_entry(matrix, i, 0) = 0;

    /* M[i, j] = M[i - 1, j - 1] - M[len - 1, j - 1] * lead * ctx->modulus->coeffs[i] */
    for (j = 1; j < len; j++) {
        nmod_mat_entry(matrix, len - 1, j) =
            nmod_mul(nmod_mat_entry(matrix, len - 1, j - 1), lead, mod);
        for (i = 0; i < len; i++) {
            nmod_mat_entry(matrix, i, j) =
                nmod_mul(nmod_mat_entry(matrix, len - 1, j),
                         modulus->coeffs[i], mod);
            if (i > 0)
                nmod_mat_entry(matrix, i, j) =
                    nmod_sub(nmod_mat_entry(matrix, i, j),
                             nmod_mat_entry(matrix, i - 1, j - 1),
                             mod);
            nmod_mat_entry(matrix, i, j) = 
                nmod_neg(nmod_mat_entry(matrix, i, j), mod);
        }
    }
}
