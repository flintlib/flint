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

/*
    Convert B to A if the degrees of A are <= expected_deg
    If not, return 0 and set A to 0.
*/
int nmod_mpoly_convert_from_nmod_mpolyd_degbound(
                                  nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx,
                           const nmod_mpolyd_t B, const nmod_mpolyd_ctx_t dctx,
                                                            slong * expect_deg)
{
    int ret;
    slong off, j, k, N;
    slong bits, nvars = ctx->minfo->nvars;
    slong Alen;
    slong * rexpect_deg, * perm = dctx->perm;
    slong perm_nontrivial;
    ulong topmask, outrange;
    ulong * exps, * pcurexp, * pexps, * rangemask;
    TMP_INIT;

    FLINT_ASSERT(nvars == B->nvars);
    FLINT_ASSERT(nvars <= FLINT_BITS);

    TMP_START;
    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    rexpect_deg = (slong *) TMP_ALLOC(nvars*sizeof(slong));

    /* find bits needed for the result */
    off = 1;
    perm_nontrivial = 0;
    for (j = 0; j < nvars; j++)
    {
        rexpect_deg[j] = expect_deg[perm[j]];
        FLINT_ASSERT(rexpect_deg[j] >= 0);
        off *= B->deg_bounds[j];
        exps[perm[j]] = rexpect_deg[j];
        perm_nontrivial |= j ^ perm[j];
    }

    FLINT_ASSERT(off <= B->coeff_alloc);

    bits = mpoly_exp_bits_required_ui(exps, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* we are going to push back terms manually */
    nmod_mpoly_fit_length_reset_bits(A, 0, bits, ctx);
    Alen = 0;

    /* find exponent vector for all variables */
    pexps = (ulong *) TMP_ALLOC(N*nvars*sizeof(ulong));
    for (k = 0; k < nvars; k++)
    {
        for (j = 0; j < nvars; j++)
            exps[perm[j]] = (j == k);
        mpoly_set_monomial_ui(pexps + k*N, exps, bits, ctx->minfo);
    }

    /* get most significant exponent in exps and its vector in ptempexp */
    off--;
    pcurexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    rangemask = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    outrange = 0;
    mpoly_monomial_zero(pcurexp, N);
    k = off;
    for (j = nvars - 1; j >= 0; j--) 
    {
        exps[j] = k % B->deg_bounds[j];
        rangemask[j] = UWORD(1) << j;
        outrange ^= (outrange ^ FLINT_SIGN_EXT(rexpect_deg[j]
                                                    - exps[j])) & rangemask[j];
        k = k / B->deg_bounds[j];
        mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + N*j, N);
    }

    /* scan down through the exponents */
    topmask = 0;
    for (; off >= 0; off--)
    {
        if (B->coeffs[off] != UWORD(0))
        {
            if (outrange)
                goto failed_out_range;

            _nmod_mpoly_fit_length(&A->coeffs, &A->coeffs_alloc,
                                   &A->exps, &A->exps_alloc, N, Alen + 1);
            A->coeffs[Alen] = B->coeffs[off];
            mpoly_monomial_set(A->exps + N*Alen, pcurexp, N);
            topmask |= (A->exps + N*Alen)[N - 1];
            Alen++;
        }

        j = nvars - 1;
        do {
            --exps[j];
            outrange ^= (outrange ^ FLINT_SIGN_EXT(rexpect_deg[j]
                                                    - exps[j])) & rangemask[j];
            if ((slong)(exps[j]) < WORD(0))
            {
                FLINT_ASSERT(off == 0 || j > 0);
                FLINT_ASSERT(exps[j] == -UWORD(1));
                exps[j] = B->deg_bounds[j] - 1;
                outrange ^= (outrange ^ FLINT_SIGN_EXT(rexpect_deg[j]
                                                    - exps[j])) & rangemask[j];
                mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + N*j, N);
            } else
            {
                mpoly_monomial_sub_mp(pcurexp, pcurexp, pexps + N*j, N);
                break;
            }
        } while (--j >= 0);
    }
    _nmod_mpoly_set_length(A, Alen, ctx);

    /* sort the exponents if needed */
    if (ctx->minfo->ord != ORD_LEX || perm_nontrivial != WORD(0))
    {
        slong msb;
        mpoly_get_cmpmask(pcurexp, N, bits, ctx->minfo);
        if (topmask != UWORD(0))
        {
            count_leading_zeros(msb, topmask);
            msb = (FLINT_BITS - 1)^msb;
        } else
        {
            msb = -WORD(1);
        }
        if (N == 1) {
            if (msb >= WORD(0))
            {
                _nmod_mpoly_radix_sort1(A, 0, A->length,
                                                   msb, pcurexp[0], topmask);
            }
        } else {
            _nmod_mpoly_radix_sort(A, 0, A->length,
                                        (N - 1)*FLINT_BITS + msb, N, pcurexp);
        }
    }

    ret = 1;

done:
    TMP_END;
    return ret;

failed_out_range:
    ret = 0;
    _nmod_mpoly_set_length(A, WORD(0), ctx);
    goto done;
}

/*
    return  -1 : function failed
             0 : B does not divide A
             1 : B divides A and quotient is in Q
*/
int nmod_mpoly_divides_dense(nmod_mpoly_t Q,
                        const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int ret, success;
    slong i;
    slong nvars = ctx->minfo->nvars;
    nmod_mpolyd_ctx_t dctx;
    nmod_mpolyd_t Ad, Bd, Qd;
    nmod_poly_t Au, Bu, Qu, Ru;
    slong * Abounds, * Bbounds, * Qbounds, * Edegs;
    TMP_INIT;

    if (B->length == 0)
    {
        if (A->length == 0 || nmod_mpoly_ctx_modulus(ctx) == 1)
        {
            nmod_mpoly_set(Q, A, ctx);
            return 1;
        }
        else
        {
            flint_throw(FLINT_DIVZERO, "nmod_mpoly_divides_dense: divide by zero");
        }
    }

    if (A->length == 0)
    {
        nmod_mpoly_zero(Q, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS ||
        ctx->minfo->nvars > FLINT_BITS ||
        ctx->minfo->nvars < 1)
    {
        return -1;
    }

    TMP_START;

    /* set the ordering of variables for the dense representation */
    nmod_mpolyd_ctx_init(dctx, nvars);

    /*
        for each variable v we need to pack to degree deg_v(A)
        except for the outermost variable
    */
    Abounds = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Bbounds = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Qbounds = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    Edegs   = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));

    mpoly_degrees_si(Abounds, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(Bbounds, B->exps, B->length, B->bits, ctx->minfo);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        /* if divides, expected degrees */
        Edegs[i] = Abounds[i] - Bbounds[i];

        if (Abounds[i] < Bbounds[i])
        {
            ret = 0;
            nmod_mpoly_zero(Q, ctx);
            goto cleanup_stage1;
        }

        if (i != dctx->perm[0])
        {
            /* variable of index i is not the outermost */
            Qbounds[i] = Abounds[i] + 1;
            Bbounds[i] = Abounds[i] + 1;
        } else
        {
            /* variable of index i is the outermost */
            Qbounds[i] = Abounds[i] - Bbounds[i] + 1;
            Bbounds[i] = Bbounds[i] + 1;
        }
        Abounds[i] = Abounds[i] + 1;            

        if (Abounds[i] < WORD(0))
        {
            ret = -1;
            nmod_mpoly_zero(Q, ctx);
            goto cleanup_stage1;
        }
    }

    nmod_mpolyd_init(Ad, nvars);
    nmod_mpolyd_init(Bd, nvars);
    nmod_mpolyd_init(Qd, nvars);

    success = 1;
    success = success && nmod_mpolyd_set_degbounds_perm(Ad, dctx, Abounds);
    success = success && nmod_mpolyd_set_degbounds_perm(Bd, dctx, Bbounds);
    success = success && nmod_mpolyd_set_degbounds_perm(Qd, dctx, Qbounds);
    if (!success)
    {
        ret = -1;
        goto cleanup_stage2;
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

    /* manually move Qd to Qu */
    Qu->alloc  = Qd->coeff_alloc;
    Qu->coeffs = Qd->coeffs;
    Qu->length = 0;
    Qu->mod    = ctx->mod;

    nmod_poly_init_mod(Ru, ctx->mod);
    nmod_poly_divrem(Qu, Ru, Au, Bu);
    if (!nmod_poly_is_zero(Ru))
    {
        ret = 0;
        goto cleanup_stage3;
    }
    nmod_poly_clear(Ru);

    /* manually move Qu to Qd */
    Qd->coeff_alloc = Qu->alloc;
    Qd->coeffs      = Qu->coeffs;
    for (i = Qu->length; i < Qd->coeff_alloc; i++)
        Qd->coeffs[i] = UWORD(0);

    nmod_mpolyd_clear(Bd);
    nmod_mpolyd_clear(Ad);
    ret = nmod_mpoly_convert_from_nmod_mpolyd_degbound(Q, ctx, Qd, dctx, Edegs);
    nmod_mpolyd_clear(Qd);
    nmod_mpolyd_ctx_clear(dctx);

done:
    TMP_END;
    return ret;

cleanup_stage3:
    nmod_poly_clear(Ru);    

cleanup_stage2:
    nmod_mpolyd_clear(Ad);
    nmod_mpolyd_clear(Bd);
    nmod_mpolyd_clear(Qd);

cleanup_stage1:
    nmod_mpolyd_ctx_clear(dctx);
    goto done;
}
