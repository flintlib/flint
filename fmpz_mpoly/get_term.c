/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void _fmpz_mpoly_get_term_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t poly,
                                  const fmpz * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong N, index, exp_bits;
    ulong * cmpmask, * packed_exp;
    int exists;
    TMP_INIT;

    exp_bits = mpoly_exp_bits_required_fmpz(exp, ctx->minfo);

    if (exp_bits > poly->bits) /* exponent too large to be poly exponent */
    {
        fmpz_zero(c);
        return;
    }

    TMP_START;
   
    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, poly->bits, ctx->minfo);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_set_monomial_fmpz(packed_exp, exp, poly->bits, ctx->minfo);

    exists = mpoly_monomial_exists(&index, poly->exps,
                                  packed_exp, poly->length, N, cmpmask);

    if (!exists)
        fmpz_zero(c);
    else
        fmpz_set(c, poly->coeffs + index);

    TMP_END; 
}

void fmpz_mpoly_get_term_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t poly,
                                 const fmpz ** exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * newexp;
    TMP_INIT;

    TMP_START;
    newexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
    {
        fmpz_init(newexp + i);
        fmpz_set(newexp + i, exp[i]);
    }

    _fmpz_mpoly_get_term_fmpz_fmpz(c, poly, newexp, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;
}

ulong fmpz_mpoly_get_term_ui_fmpz(const fmpz_mpoly_t poly,
                                 const fmpz ** exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t newc;
    ulong ret;

    fmpz_init(newc);
    fmpz_mpoly_get_term_fmpz_fmpz(newc, poly, exp, ctx);

    ret = fmpz_get_ui(newc);
    fmpz_clear(newc);
    return ret;
}

slong fmpz_mpoly_get_term_si_fmpz(const fmpz_mpoly_t poly,
                                 const fmpz ** exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t newc;
    slong ret;

    fmpz_init(newc);
    fmpz_mpoly_get_term_fmpz_fmpz(newc, poly, exp, ctx);

    ret = fmpz_get_si(newc);
    fmpz_clear(newc);
    return ret;
}

void fmpz_mpoly_get_term_fmpz_ui(fmpz_t c, const fmpz_mpoly_t poly,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * newexp;
    TMP_INIT;

    TMP_START;
    newexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
        fmpz_init_set_ui(newexp + i, exp[i]);

    _fmpz_mpoly_get_term_fmpz_fmpz(c, poly, newexp, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;
}

ulong fmpz_mpoly_get_term_ui_ui(const fmpz_mpoly_t poly,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t newc;
    ulong ret;

    fmpz_init(newc);
    fmpz_mpoly_get_term_fmpz_ui(newc, poly, exp, ctx);

    ret = fmpz_get_ui(newc);
    fmpz_clear(newc);
    return ret;
}

slong fmpz_mpoly_get_term_si_ui(const fmpz_mpoly_t poly,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t newc;
    slong ret;

    fmpz_init(newc);
    fmpz_mpoly_get_term_fmpz_ui(newc, poly, exp, ctx);

    ret = fmpz_get_si(newc);
    fmpz_clear(newc);
    return ret;
}

