/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "flint-impl.h"

void _fq_nmod_dense_reduce(mp_limb_t* R, slong lenR, const fq_nmod_ctx_t ctx)
{
    mp_limb_t  *q, *r;

    if (lenR < ctx->modulus->length)
    {
        _nmod_vec_reduce(R, R, lenR, ctx->mod);
        return;
    }

    q = _nmod_vec_init(lenR - ctx->modulus->length + 1);
    r = _nmod_vec_init(ctx->modulus->length - 1);

    _nmod_poly_divrem_newton_n_preinv(q, r, R, lenR, 
                                      ctx->modulus->coeffs, ctx->modulus->length,
                                      ctx->inv->coeffs, ctx->inv->length,
                                      ctx->mod);

    _NMOD_VEC_SET(R, r, ctx->modulus->length - 1);
    _nmod_vec_clear(q);
    _nmod_vec_clear(r);

}

