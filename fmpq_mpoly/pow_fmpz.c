/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_pow_fmpz(fmpq_mpoly_t qpoly1, const fmpq_mpoly_t qpoly2,
                                 const  fmpz_t pow, const fmpq_mpoly_ctx_t qctx)
{
    slong i;
    fmpq_t c;
    fmpz * max_fields2;
    mp_bitcnt_t exp_bits;
    fmpz_mpoly_struct * poly1 = qpoly1->zpoly;
    const fmpz_mpoly_struct * poly2 = qpoly2->zpoly;
    const fmpz_mpoly_ctx_struct * ctx = qctx->zctx;
    TMP_INIT;

    if (fmpz_sgn(pow) < 0)
        flint_throw(FLINT_ERROR, "Negative power in fmpq_mpoly_pow_fmpz");

    if (fmpz_fits_si(pow))
    {
        if (fmpz_is_one(fmpq_denref(qpoly2->content))
              && fmpz_is_pm1(fmpq_numref(qpoly2->content)))
        {
            if (fmpz_is_one(fmpq_numref(qpoly2->content))
                || fmpz_is_even(pow))
            {
                fmpq_set_si(qpoly1->content, +WORD(1), UWORD(1));
            } else {
                fmpq_set_si(qpoly1->content, -WORD(1), UWORD(1));
            }

        } else {
            fmpq_pow_si(qpoly1->content, qpoly2->content, fmpz_get_si(pow));
        }

        fmpz_mpoly_pow_fps(poly1, poly2, fmpz_get_ui(pow), ctx);
        return;
    }

    /*
        we are raising a polynomial to an unreasonable exponent.
        It must either be zero or a monomial with coefficient +-1.
    */

    if (fmpz_mpoly_is_zero(poly2, ctx))
    {
        fmpq_mpoly_zero(qpoly1, qctx);
        return;
    }

    if (poly2->length != WORD(1))
        flint_throw(FLINT_ERROR, "Multinomial in fmpz_mpoly_pow_fmpz");


    fmpq_init(c);
    fmpq_abs(c, qpoly2->content);
    if (!fmpq_is_one(c))
        flint_throw(FLINT_ERROR, "Non-unit coefficient in fmpz_mpoly_pow_fmpz");
    fmpq_clear(c);

    TMP_START;

    max_fields2 = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(max_fields2 + i);
    mpoly_max_fields_fmpz(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);

    _fmpz_vec_scalar_mul_fmpz(max_fields2, max_fields2, ctx->minfo->nfields, pow);

    exp_bits = _fmpz_vec_max_bits(max_fields2, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    fmpz_mpoly_fit_length(poly1, 1, ctx);
    fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
    poly1->bits = exp_bits;

    fmpz_set_si(poly1->coeffs + 0, UWORD(1));
    if (fmpq_is_one(qpoly2->content) || fmpz_is_even(pow))
        fmpq_set_si(qpoly1->content, +WORD(1), UWORD(1));
    else
        fmpq_set_si(qpoly1->content, -WORD(1), UWORD(1));

    mpoly_pack_vec_fmpz(poly1->exps + 0, max_fields2, exp_bits, ctx->minfo->nfields, 1);

    _fmpz_mpoly_set_length(poly1, 1, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_clear(max_fields2 + i);

    TMP_END;
}
