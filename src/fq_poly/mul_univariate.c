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
_fq_poly_mul_univariate (fq_struct * rop,
                        const fq_struct * op1, slong len1,
                        const fq_struct * op2, slong len2,
                        const fq_ctx_t ctx)
{
    const slong fqlen = ctx->modulus->length - 1;
    const slong pfqlen = 2*fqlen - 1;
    const slong rlen = len1 + len2 - 1;
    const slong llen1 = (op1 + (len1 - 1))->length;
    const slong llen2 = (op2 + (len2 - 1))->length;
    const slong clen1 = pfqlen*(len1-1) + llen1;
    const slong clen2 = pfqlen*(len2-1) + llen2;
    const slong crlen = clen1 + clen2 - 1;
    const slong lrlen = llen1 + llen2 - 1;
    slong i;
    slong len;

    fmpz *cop1, *cop2, *crop;

    cop1 = _fmpz_vec_init(clen1);
    for (i = 0; i < len1 - 1; i++)
    {
        _fmpz_vec_set(cop1 + pfqlen*i, (op1 + i)->coeffs, (op1 + i)->length);
        _fmpz_vec_zero(cop1 + pfqlen*i + (op1 + i)->length, pfqlen - (op1 + i)->length);
    }
    {
        _fmpz_vec_set(cop1 + pfqlen*i, (op1 + i)->coeffs, (op1 + i)->length);
    }

    if (op2 != op1)
    {
        cop2 = _fmpz_vec_init(clen2);
        for (i = 0; i < len2 - 1; i++)
        {
            _fmpz_vec_set(cop2 + pfqlen*i, (op2 + i)->coeffs,(op2 + i)->length);
            _fmpz_vec_zero(cop2 + pfqlen*i + (op2 + i)->length, pfqlen - (op2 + i)->length);
        }
        {
            _fmpz_vec_set(cop2 + pfqlen*i, (op2 + i)->coeffs, (op2 + i)->length);
        }
    }
    else
    {
        cop2 = cop1;
    }

    crop = _fmpz_vec_init(crlen);
    if (clen1 >= clen2)
        _fmpz_poly_mul(crop, cop1, clen1, cop2, clen2);
    else
        _fmpz_poly_mul(crop, cop2, clen2, cop1, clen1);

    for (i = 0; i < rlen - 1; i++)
    {
        _fq_reduce(crop + pfqlen*i, pfqlen, ctx);
        len = fqlen;
        while (len && !(*(crop + pfqlen*i + len - 1))) len--;
        fmpz_poly_fit_length(rop + i, len);
        (rop + i)->length = len;
        _fmpz_vec_set((rop + i)->coeffs, crop + pfqlen*i, len);
    }
    {
        _fq_reduce(crop + pfqlen*i, lrlen, ctx);
        len = FLINT_MIN(fqlen, lrlen);
        while (len && !(*(crop + pfqlen*i + len - 1))) len--;
        fmpz_poly_fit_length(rop + i, len);
        (rop + i)->length = len;
        _fmpz_vec_set((rop + i)->coeffs, crop + pfqlen*i, len);
    }

    _fmpz_vec_clear(cop1, clen1);
    if (op2 != op1)
    {
        _fmpz_vec_clear(cop2, clen2);
    }
    _fmpz_vec_clear(crop, crlen);
}

void
fq_poly_mul_univariate (fq_poly_t rop,
                       const fq_poly_t op1,
                       const fq_poly_t op2,
                       const fq_ctx_t ctx)
{
    const slong len1 = op1->length;
    const slong len2 = op2->length;
    const slong rlen = op1->length + op2->length - 1;

    if (len1 == 0 || len2 == 0)
    {
        fq_poly_zero(rop, ctx);
        return;
    }
    else
    {
        fq_poly_fit_length(rop, rlen, ctx);
        _fq_poly_mul_univariate(rop->coeffs, op1->coeffs, len1,
                                   op2->coeffs, len2, ctx);
        _fq_poly_set_length(rop, rlen, ctx);
    }
}

