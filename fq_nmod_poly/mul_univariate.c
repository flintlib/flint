/*
    Copyright (C) 2016 Jean-Pierre Flori

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_poly.h"

void
_fq_nmod_poly_mul_univariate_no_pad (fq_nmod_struct * rop,
                        const fq_nmod_struct * op1, slong len1,
                        const fq_nmod_struct * op2, slong len2,
                        const fq_nmod_ctx_t ctx)
{
    const slong fqlen = ctx->modulus->length - 1;
    const slong pfqlen = 2*fqlen - 1;
    const nmod_t mod = ctx->mod;
    const slong rlen = len1 + len2 - 1;
    const slong llen1 = (op1 + (len1 - 1))->length;
    const slong llen2 = (op2 + (len2 - 1))->length;
    const slong clen1 = pfqlen*(len1-1) + llen1;
    const slong clen2 = pfqlen*(len2-1) + llen2;
    const slong crlen = clen1 + clen2 - 1;
    const slong lrlen = llen1 + llen2 - 1;
    slong i;
    slong len;

    mp_ptr cop1, cop2, crop;

    cop1 = (mp_limb_t *) flint_malloc(clen1*sizeof(mp_limb_t));
    for (i = 0; i < len1 - 1; i++)
    {
        flint_mpn_copyi(cop1 + pfqlen*i, (op1 + i)->coeffs, (op1 + i)->length);
        flint_mpn_zero(cop1 + pfqlen*i + (op1 + i)->length, pfqlen - (op1 + i)->length);
    }
    {
        flint_mpn_copyi(cop1 + pfqlen*i, (op1 + i)->coeffs, (op1 + i)->length);
    }

    if (op2 != op1)
    {
        cop2 = (mp_limb_t *) flint_malloc(clen2*sizeof(mp_limb_t));
        for (i = 0; i < len2 - 1; i++)
        {
            flint_mpn_copyi(cop2 + pfqlen*i, (op2 + i)->coeffs,(op2 + i)->length);
            flint_mpn_zero(cop2 + pfqlen*i + (op2 + i)->length, pfqlen - (op2 + i)->length);
        }
        {
            flint_mpn_copyi(cop2 + pfqlen*i, (op2 + i)->coeffs, (op2 + i)->length);
        }
    }
    else
    {
        cop2 = cop1;
    }

    crop = (mp_limb_t *) flint_malloc(crlen*sizeof(mp_limb_t));

    if (clen1 >= clen2)
        _nmod_poly_mul(crop, cop1, clen1, cop2, clen2, mod);
    else
        _nmod_poly_mul(crop, cop2, clen2, cop1, clen1, mod);

    for (i = 0; i < rlen - 1; i++)
    {
        _fq_nmod_reduce(crop + pfqlen*i, pfqlen, ctx);
        len = fqlen;
        while (len && (*(crop + pfqlen*i + len - 1) == UWORD(0))) len--;
        nmod_poly_fit_length(rop + i, len);
        (rop + i)->length = len;
        flint_mpn_copyi((rop + i)->coeffs, crop + pfqlen*i, len);
    }
    {
        len = lrlen;
        if (len > fqlen)
        {
            _fq_nmod_reduce(crop + pfqlen*i, len, ctx);
            len = fqlen;
            while (len && (*(crop + pfqlen*i + len - 1) == UWORD(0))) len--;
        }
        nmod_poly_fit_length(rop + i, len);
        (rop + i)->length = len;
        flint_mpn_copyi((rop + i)->coeffs, crop + pfqlen*i, len);
    }

    flint_free(cop1);
    if (op2 != op1)
    {
        flint_free(cop2);
    }
    flint_free(crop);
}

void
_fq_nmod_poly_mul_univariate (fq_nmod_struct * rop,                                             const fq_nmod_struct * op1, slong len1,                                        const fq_nmod_struct * op2, slong len2,                                        const fq_nmod_ctx_t ctx)
{
   slong len1n = len1, len2n = len2;

   _fq_nmod_poly_normalise2(op1, &len1n, ctx);
   _fq_nmod_poly_normalise2(op2, &len2n, ctx);

   if (len1n == 0)
      _fq_nmod_vec_zero(rop, len1 + len2 - 1, ctx);
   else if (len2n == 0)
      _fq_nmod_vec_zero(rop, len1 + len2 - 1, ctx);
   else
   {
       _fq_nmod_poly_mul_univariate_no_pad(rop, op1, len1n, op2, len2n, ctx);

       _fq_nmod_vec_zero(rop + len1n + len2n - 1,
		      len1 + len2 - len1n - len2n, ctx);
   }
}

void
fq_nmod_poly_mul_univariate (fq_nmod_poly_t rop,
                       const fq_nmod_poly_t op1,
                       const fq_nmod_poly_t op2,
                       const fq_nmod_ctx_t ctx)
{
    const slong len1 = op1->length;
    const slong len2 = op2->length;
    const slong rlen = op1->length + op2->length - 1;

    if (len1 == 0 || len2 == 0)
    {
        fq_nmod_poly_zero(rop, ctx);
        return;
    }
    else
    {
        fq_nmod_poly_fit_length(rop, rlen, ctx);
        _fq_nmod_poly_mul_univariate(rop->coeffs, op1->coeffs, len1,
                                   op2->coeffs, len2, ctx);
        _fq_nmod_poly_set_length(rop, rlen, ctx);
    }
}

