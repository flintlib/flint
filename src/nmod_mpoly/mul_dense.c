/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_mpoly.h"

int _nmod_mpoly_mul_dense(nmod_mpoly_t P,
                                 const nmod_mpoly_t A, fmpz * maxAfields,
                                 const nmod_mpoly_t B, fmpz * maxBfields,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i;
    slong nvars = ctx->minfo->nvars;
    nmod_mpolyd_ctx_t dctx;
    nmod_mpolyd_t Ad, Bd, Pd;
    nmod_poly_t Au, Bu, Pu;
    slong * Abounds, * Bbounds, * Pbounds;
    TMP_INIT;

    FLINT_ASSERT(A->length != 0);
    FLINT_ASSERT(B->length != 0);
    FLINT_ASSERT(nvars > 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    TMP_START;

    /* set the ordering of variables for the dense representation */
    nmod_mpolyd_ctx_init(dctx, nvars);

    /*
        for each variable v except for the outermost variable,
        we need to pack to degree deg_v(A) + deg_v(B)
    */
    Abounds = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Bbounds = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Pbounds = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Abounds, maxAfields, ctx->minfo);
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Bbounds, maxBfields, ctx->minfo);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        Abounds[i] = Abounds[i] + 1;
        Bbounds[i] = Bbounds[i] + 1;
        Pbounds[i] = Abounds[i] + Bbounds[i] - 1;
        if ((Abounds[i] | Bbounds[i] | Pbounds[i]) < WORD(0))
        {
            goto failed_stage1;
        }
        if (i != dctx->perm[0])
        {
            /* variable of index i is not the outermost */
            Abounds[i] = Pbounds[i];
            Bbounds[i] = Pbounds[i];
        }
    }

    nmod_mpolyd_init(Ad, nvars);
    nmod_mpolyd_init(Bd, nvars);
    nmod_mpolyd_init(Pd, nvars);

    success = 1;
    success = success && nmod_mpolyd_set_degbounds_perm(Ad, dctx, Abounds);
    success = success && nmod_mpolyd_set_degbounds_perm(Bd, dctx, Bbounds);
    success = success && nmod_mpolyd_set_degbounds_perm(Pd, dctx, Pbounds);
    if (!success)
    {
        goto failed_stage2;
    }

    nmod_mpoly_convert_to_nmod_mpolyd_degbound(Ad, dctx, A, ctx);
    nmod_mpoly_convert_to_nmod_mpolyd_degbound(Bd, dctx, B, ctx);

    /* let Au and Bu borrow Ad and Bd */
    Au->alloc  = Ad->coeff_alloc;
    Au->coeffs = Ad->coeffs;
    Au->length = nmod_mpolyd_length(Ad);
    Au->mod    = ctx->mod;

    Bu->alloc  = Bd->coeff_alloc;
    Bu->coeffs = Bd->coeffs;
    Bu->length = nmod_mpolyd_length(Bd);
    Bu->mod    = ctx->mod;

    /* manually move P to Pu */
    Pu->alloc  = Pd->coeff_alloc;
    Pu->coeffs = Pd->coeffs;
    Pu->length = 0;
    Pu->mod    = ctx->mod;

    nmod_poly_mul(Pu, Au, Bu);

    /* manually move Pu to P */
    Pd->coeff_alloc = Pu->alloc;
    Pd->coeffs      = Pu->coeffs;
    for (i = Pu->length; i < Pd->coeff_alloc; i++)
        Pd->coeffs[i] = UWORD(0);

    nmod_mpolyd_clear(Bd);
    nmod_mpolyd_clear(Ad);
    nmod_mpoly_convert_from_nmod_mpolyd(P, ctx, Pd, dctx);
    nmod_mpolyd_clear(Pd);

    nmod_mpolyd_ctx_clear(dctx);

done:
    TMP_END;
    return success;

failed_stage2:
    nmod_mpolyd_clear(Ad);
    nmod_mpolyd_clear(Bd);
    nmod_mpolyd_clear(Pd);

failed_stage1:
    nmod_mpolyd_ctx_clear(dctx);
    success = 0;
    goto done;
}


int nmod_mpoly_mul_dense(nmod_mpoly_t A, const nmod_mpoly_t B,
                              const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    int success;
    fmpz * maxBfields, * maxCfields;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
        return 1;
    }

    if (B->bits > FLINT_BITS || C->bits > FLINT_BITS ||
        ctx->minfo->nvars < 1)
    {
        return 0;
    }

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxCfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(maxBfields + i);
        fmpz_init(maxCfields + i);
    }
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);

    success = _nmod_mpoly_mul_dense(A, B, maxBfields, C, maxCfields, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
    return success;
}
