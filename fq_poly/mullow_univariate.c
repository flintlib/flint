/*
    Copyright (C) 2016 Jean-Pierre Flori

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_poly.h"

void
_fq_poly_mullow_univariate (fq_struct * rop,
                        const fq_struct * op1, slong len1,
                        const fq_struct * op2, slong len2,
                        slong n, const fq_ctx_t ctx)
{
    const slong fqlen = ctx->modulus->length - 1;
    const slong pfqlen = 2*fqlen - 1;
    const slong rlen = len1 + len2 - 1;
    const slong m = FLINT_MIN(n, rlen);
    const slong cmlen = pfqlen*m;
    const slong clen1 = pfqlen*len1;
    const slong clen2 = pfqlen*len2;
    slong i;
    slong len;

    fmpz *cop1, *cop2, *crop;

    if (!len1 || !len2)
    {
        _fq_poly_zero(rop, n, ctx);
        return;
    }

    cop1 = _fmpz_vec_init(clen1);
    for (i = 0; i < len1; i++)
    {
        _fmpz_vec_set(cop1 + pfqlen*i, (op1 + i)->coeffs, (op1 + i)->length);
        _fmpz_vec_zero(cop1 + pfqlen*i + (op1 + i)->length, pfqlen - (op1 + i)->length);
    }

    if (op2 != op1)
    {
        cop2 = _fmpz_vec_init(clen2);
        for (i = 0; i < len2; i++)
        {
            _fmpz_vec_set(cop2 + pfqlen*i, (op2 + i)->coeffs,(op2 + i)->length);
            _fmpz_vec_zero(cop2 + pfqlen*i + (op2 + i)->length, pfqlen - (op2 + i)->length);
        }
    }
    else
    {
        cop2 = cop1;
    }

    crop = _fmpz_vec_init(cmlen);
    if (clen1 >= clen2)
        _fmpz_poly_mullow(crop, cop1, clen1, cop2, clen2, cmlen);
    else
        _fmpz_poly_mullow(crop, cop2, clen2, cop1, clen1, cmlen);

    for (i = 0; i < m; i++)
    {
        _fq_reduce(crop + pfqlen*i, pfqlen, ctx);
        len = fqlen;
        while (len && !(*(crop + pfqlen*i + len - 1))) len--;
        fmpz_poly_fit_length(rop + i, len);
        (rop + i)->length = len;
        _fmpz_vec_set((rop + i)->coeffs, crop + pfqlen*i, len);
    }
    for (; i < n; i++)
    {
        (rop + i)->length = 0;
    }

    _fmpz_vec_clear(cop1, clen1);
    if (op2 != op1)
    {
        _fmpz_vec_clear(cop2, clen2);
    }
    _fmpz_vec_clear(crop, cmlen);
}

void
fq_poly_mullow_univariate (fq_poly_t rop,
                       const fq_poly_t op1,
                       const fq_poly_t op2,
                       slong n, const fq_ctx_t ctx)
{
    const slong len1 = op1->length;
    const slong len2 = op2->length;
    const slong rlen = op1->length + op2->length - 1;
    const slong m = FLINT_MIN(n, rlen);

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        fq_poly_zero(rop, ctx);
        return;
    }
    else
    {
        fq_poly_fit_length(rop, m, ctx);
        _fq_poly_mullow_univariate(rop->coeffs, op1->coeffs, len1,
                                   op2->coeffs, len2, m, ctx);
        _fq_poly_set_length(rop, m, ctx);
        _fq_poly_normalise(rop, ctx);
    }
}

