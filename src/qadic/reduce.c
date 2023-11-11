/*
    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "padic.h"
#include "qadic.h"

void _fmpz_poly_reduce(fmpz * R, slong lenR, const fmpz * a, const slong * j, slong len)
{
    const slong d = j[len - 1];
    slong i, k;

    FMPZ_VEC_NORM(R, lenR);

    for (i = lenR - 1; i >= d; i--)
    {
        for (k = len - 2; k >= 0; k--)
        {
            fmpz_submul(R + j[k] + i - d, R + i, a + k);
        }
        fmpz_zero(R + i);
    }
}

void _fmpz_mod_poly_reduce(fmpz * R, slong lenR, const fmpz * a, const slong * j, slong len, const fmpz_t p)
{
    const slong d = j[len - 1];

    if (lenR > d)
    {
        _fmpz_poly_reduce(R, lenR, a, j, len);
        _fmpz_vec_scalar_mod_fmpz(R, R, d, p);
    }
    else
    {
        _fmpz_vec_scalar_mod_fmpz(R, R, lenR, p);
    }
}

void qadic_reduce(qadic_t x, const qadic_ctx_t ctx)
{
    const slong N = qadic_prec(x);
    const slong d = ctx->j[ctx->len - 1];

    if (x->length == 0 || x->val >= N)
    {
        padic_poly_zero(x);
    }
    else
    {
        fmpz_t pow;
        int alloc;

        alloc = _padic_ctx_pow_ui(pow, N - x->val, &ctx->pctx);

        _fmpz_mod_poly_reduce(x->coeffs, x->length, ctx->a, ctx->j, ctx->len, pow);
        _padic_poly_set_length(x, FLINT_MIN(x->length, d));
        _padic_poly_normalise(x);
        padic_poly_canonicalise(x, (&ctx->pctx)->p);

        if (alloc)
            fmpz_clear(pow);
    }
}
