/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

void
fq_zech_get_nmod_poly(nmod_poly_t rop, const fq_zech_t op,
                                                       const fq_zech_ctx_t ctx)
{
    slong i;
    mp_limb_t q, r;

    rop->mod = ctx->fq_nmod_ctx->modulus->mod;
    
    nmod_poly_fit_length(rop, fq_zech_ctx_degree(ctx));

    q = ctx->eval_table[op->value];
    i = 0;
    while (q >= ctx->p)
    {
        r = n_divrem2_precomp(&q, q, ctx->p, ctx->ppre);
        nmod_poly_set_coeff_ui(rop, i, r);
        i ++;
    }
    nmod_poly_set_coeff_ui(rop, i, q);
}
