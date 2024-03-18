/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"
#include "fq.h"

void _fq_sparse_reduce(fmpz *R, slong lenR, const fq_ctx_t ctx)
{
    const slong d = ctx->j[ctx->len - 1];

    FMPZ_VEC_NORM(R, lenR);

    if (lenR > d)
    {
        slong i, k;

        for (i = lenR - 1; i >= d; i--)
        {
            for (k = ctx->len - 2; k >= 0; k--)
            {
                fmpz_submul(R + ctx->j[k] + i - d, R + i, ctx->a + k);
            }
            fmpz_zero(R + i);
        }
    }

    _fmpz_mod_vec_set_fmpz_vec(R, R, FLINT_MIN(d, lenR), ctx->ctxp);
}

void _fq_dense_reduce(fmpz* R, slong lenR, const fq_ctx_t ctx)
{
    fmpz  *q, *r;

    if (lenR < ctx->modulus->length)
    {
        _fmpz_mod_vec_set_fmpz_vec(R, R, lenR, ctx->ctxp);
        return;
    }

    q = _fmpz_vec_init(lenR - ctx->modulus->length + 1);
    r = _fmpz_vec_init(ctx->modulus->length - 1);

    _fmpz_mod_vec_set_fmpz_vec(R, R, lenR, ctx->ctxp);
    _fmpz_mod_poly_divrem_newton_n_preinv(q, r, R, lenR,
                                        ctx->modulus->coeffs, ctx->modulus->length,
                                        ctx->inv->coeffs, ctx->inv->length,
                                        ctx->ctxp);

    _fmpz_vec_set(R, r, ctx->modulus->length - 1);
    _fmpz_vec_clear(q, lenR - ctx->modulus->length + 1);
    _fmpz_vec_clear(r, ctx->modulus->length - 1);

}

void _fq_reduce(fmpz* R, slong lenR, const fq_ctx_t ctx)
{
    if (ctx->sparse_modulus)
        _fq_sparse_reduce(R, lenR, ctx);
    else
        _fq_dense_reduce(R, lenR, ctx);
}

void fq_reduce(fq_t rop, const fq_ctx_t ctx)
{
    _fq_reduce(rop->coeffs, rop->length, ctx);
    rop->length = FLINT_MIN(rop->length, ctx->modulus->length - 1);
    _fmpz_poly_normalise(rop);
}
