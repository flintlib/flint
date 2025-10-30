/*
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mpoly.h"
#include "fmpz_mpoly.h"
#include "mpoly.h"

void
fmpz_mod_mpoly_set_fmpz_mpoly(fmpz_mod_mpoly_t res, const fmpz_mpoly_t f, fmpz_mod_mpoly_ctx_t ctxm, fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(f->bits, ctx->minfo);
    slong res_len, i;
    FLINT_ASSERT(ctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(ctx->minfo->ord == ctx->minfo->ord);
    fmpz_mod_mpoly_fit_length_reset_bits(res, f->length, f->bits, ctxm);
    res_len = 0;
    for (i = 0; i < f->length; i++)
    {

        fmpz_mod_set_fmpz(res->coeffs+res_len,f->coeffs + i, ctxm->ffinfo);

        if (res->coeffs[res_len] == 0)
            continue;

        mpoly_monomial_set(res->exps + N*res_len, f->exps + N*i, N);
        res_len++;
    }
    res->length = res_len;
}

void
fmpz_mod_mpoly_get_fmpz_mpoly(fmpz_mpoly_t res, const fmpz_mod_mpoly_t f, fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(f->bits, ctx->minfo);
    slong res_len, i;
    FLINT_ASSERT(ctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(ctx->minfo->ord == ctx->minfo->ord);
    fmpz_mpoly_fit_length_reset_bits(res, f->length, f->bits, ctx);
    res_len = 0;
    for (i = 0; i < f->length; i++)
    {
        fmpz_set(res->coeffs+res_len, f->coeffs + i);
        
        if (res->coeffs[res_len] == 0)
            continue;
            
        mpoly_monomial_set(res->exps + N*res_len, f->exps + N*i, N);
        res_len++;
    }
    res->length = res_len;
}
