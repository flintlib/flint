/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void _fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t qpoly,
                                  const fmpz * exp, const fmpq_mpoly_ctx_t qctx)
{
    slong N, index, exp_bits;
    ulong * cmpmask, * packed_exp;
    int exists;
    const fmpz_mpoly_struct * poly = qpoly->zpoly;
    const fmpz_mpoly_ctx_struct * ctx = qctx->zctx;
    TMP_INIT;

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);

    if (exp_bits > poly->bits) /* exponent too large to be poly exponent */
    {
        fmpq_zero(c);
        return;
    }

    TMP_START;
   
    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, poly->bits, ctx->minfo);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_set_monomial_ffmpz(packed_exp, exp, poly->bits, ctx->minfo);

    exists = mpoly_monomial_exists(&index, poly->exps,
                                  packed_exp, poly->length, N, cmpmask);

    if (!exists)
    {
        fmpq_zero(c);
    } else
    {
        fmpq_mul_fmpz(c, qpoly->content, poly->coeffs + index);
    }

    TMP_END; 
}

void fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t poly,
                                fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->zctx->minfo->nvars;
    fmpz * newexp;
    TMP_INIT;

    TMP_START;
    newexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
    {
        fmpz_init(newexp + i);
        fmpz_set(newexp + i, exp[i]);
    }

    _fmpq_mpoly_get_coeff_fmpq_fmpz(c, poly, newexp, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;
}
