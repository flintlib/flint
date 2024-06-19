/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fq_nmod.h"

void _fq_nmod_sparse_reduce(ulong *R, slong lenR, const fq_nmod_ctx_t ctx)
{
    slong i, k;
    const slong d = ctx->j[ctx->len - 1];

    NMOD_VEC_NORM(R, lenR);

    for (i = lenR - 1; i >= d; i--)
    {
        for (k = ctx->len - 2; k >= 0; k--)
        {
            /* TODO clean this mess up */
            R[ctx->j[k] + i - d] = n_submod(R[ctx->j[k] + i - d],
                                            n_mulmod2_preinv(R[i], ctx->a[k], ctx->mod.n, ctx->mod.ninv),
                                            ctx->mod.n);
        }
        R[i] = UWORD(0);
    }
}

void _fq_nmod_dense_reduce(ulong* R, slong lenR, const fq_nmod_ctx_t ctx)
{
    ulong  *q, *r;

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

    _nmod_vec_set(R, r, ctx->modulus->length - 1);
    _nmod_vec_clear(q);
    _nmod_vec_clear(r);

}

void _fq_nmod_reduce(ulong* R, slong lenR, const fq_nmod_ctx_t ctx)
{
    if (ctx->sparse_modulus)
        _fq_nmod_sparse_reduce(R, lenR, ctx);
    else
        _fq_nmod_dense_reduce(R, lenR, ctx);
}

void fq_nmod_reduce(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(rop->length <= 2*(ctx->modulus->length - 1));
    _fq_nmod_reduce(rop->coeffs, rop->length, ctx);
    rop->length = FLINT_MIN(rop->length, ctx->modulus->length - 1);
    _nmod_poly_normalise(rop);
}
