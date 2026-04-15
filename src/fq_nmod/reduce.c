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

/* Todo: merge code with _nmod_poly_divrem_try_sparse */

void _fq_nmod_sparse_reduce(ulong *R, slong lenR, const fq_nmod_ctx_t ctx)
{
    slong i, k;
    ulong c, cbound;
    const slong d = ctx->j[ctx->len - 1];
    nmod_t mod = ctx->mod;

    cbound = 0;
    for (i = 0; i < ctx->len - 1; i++)
        cbound |= ctx->a[i];

    NMOD_VEC_NORM(R, lenR);

    if (ctx->mod.n == 2)
    {
        for (i = lenR - 1; i >= d; i--)
        {
            c = R[i];

            for (k = ctx->len - 2; k >= 0; k--)
            {
                R[ctx->j[k] + i - d] ^= c;
            }
            R[i] = 0;
        }
    }
    else if (cbound == 1)
    {
        for (i = lenR - 1; i >= d; i--)
        {
            c = R[i];

            for (k = ctx->len - 2; k >= 0; k--)
            {
                R[ctx->j[k] + i - d] = nmod_sub(R[ctx->j[k] + i - d], c, mod);
            }
            R[i] = 0;
        }
    }
    else if (lenR > d + 2 && NMOD_BITS(mod) < FLINT_BITS / 2)
    {
        ulong ninv = n_barrett_precomp(mod.n);

        for (i = lenR - 1; i >= d; i--)
        {
            for (k = ctx->len - 2; k >= 0; k--)
            {
                R[ctx->j[k] + i - d] = n_mod_barrett(R[ctx->j[k] + i - d]
                    + R[i] * (mod.n - ctx->a[k]), mod.n, ninv);
            }
            R[i] = 0;
        }
    }
    else
    {
        for (i = lenR - 1; i >= d; i--)
        {
            for (k = ctx->len - 2; k >= 0; k--)
            {
                R[ctx->j[k] + i - d] = nmod_addmul(R[ctx->j[k] + i - d], R[i], mod.n - ctx->a[k], ctx->mod);
            }
            R[i] = 0;
        }
    }
}

void _fq_nmod_dense_reduce(ulong* R, slong lenR, const fq_nmod_ctx_t ctx)
{
    ulong  *q, *r;
    slong lenq, lenr;

    if (lenR < ctx->modulus->length)
        return;

    TMP_INIT;
    TMP_START;

    lenq = (lenR - ctx->modulus->length + 1);
    lenr = ctx->modulus->length - 1;
    q = TMP_ALLOC((lenq + lenr) * sizeof(ulong));
    r = q + lenq;

    _nmod_poly_divrem_newton_n_preinv(q, r, R, lenR,
                                      ctx->modulus->coeffs, ctx->modulus->length,
                                      ctx->inv->coeffs, ctx->inv->length,
                                      ctx->mod);

    _nmod_vec_set(R, r, lenr);

    TMP_END;

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

