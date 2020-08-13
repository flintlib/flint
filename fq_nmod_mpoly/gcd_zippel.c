/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


int fq_nmod_mpolyu_gcdm_zippel_bivar(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t Abar,
    fq_nmod_mpolyu_t Bbar,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    slong var = 0;
    slong Alastdeg, Blastdeg;
    slong bound;
    slong lastdeg;
    int success = 0, changed, have_enough;
    fq_nmod_mpolyun_t An, Bn, H, Ht;
    fq_nmod_poly_t modulus, gamma, hc;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx;
    fq_nmod_mpolyu_t Aeval, Beval, Geval;
    fq_nmod_t t, geval;

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(Abar->bits == A->bits);
    FLINT_ASSERT(Bbar->bits == A->bits);
    FLINT_ASSERT(B->bits == A->bits);

    fq_nmod_mpolyun_init(An, A->bits, ctx);
    fq_nmod_mpolyun_init(Bn, A->bits, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    FLINT_ASSERT(An->bits == B->bits);
    FLINT_ASSERT(An->bits == G->bits);
    FLINT_ASSERT(An->length > 0);
    FLINT_ASSERT(Bn->length > 0);
    FLINT_ASSERT(An->exps[A->length - 1] == 0);
    FLINT_ASSERT(Bn->exps[B->length - 1] == 0);

    fq_nmod_poly_init(gamma, ctx->fqctx);
    fq_nmod_poly_gcd(gamma, fq_nmod_mpolyun_leadcoeff_poly(An, ctx),
                            fq_nmod_mpolyun_leadcoeff_poly(Bn, ctx), ctx->fqctx);

    /* bound on the number of images */
    Alastdeg = fq_nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = fq_nmod_mpolyun_lastdeg(Bn, ctx);
    bound = 1 + fq_nmod_poly_degree(gamma, ctx->fqctx)
              + FLINT_MIN(Alastdeg, Blastdeg);

    fq_nmod_poly_init(hc, ctx->fqctx);
    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_mpolyun_init(H, A->bits, ctx);
    fq_nmod_mpolyun_init(Ht, A->bits, ctx);

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, ectx, ctx, randstate);

    /*
        Once Aeval, Beval, ..., t are inited in ectx->fqctx, they do not need
        to be cleared and reinited when ectx->fqctx changes.
    */
    fq_nmod_mpolyu_init(Aeval, A->bits, ectx);
    fq_nmod_mpolyu_init(Beval, A->bits, ectx);
    fq_nmod_mpolyu_init(Geval, A->bits, ectx);
    fq_nmod_init(geval, ectx->fqctx);
    fq_nmod_init(t, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime_outer:

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, randstate);
    if (cur_emb == NULL)
    {
        /* ran out of primes */
        success = 0;
        goto finished;
    }

have_prime:

    FLINT_ASSERT(cur_emb != NULL);

    /* make sure reduction does not kill both lc(A) and lc(B) */
    bad_fq_nmod_embed_sm_to_lg(geval, gamma, cur_emb);
    if (fq_nmod_is_zero(geval, ectx->fqctx))
        goto choose_prime_outer;

    /* make sure reduction does not kill either A or B */
    fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Aeval, An, ectx, ctx, cur_emb);
    fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Beval, Bn, ectx, ctx, cur_emb);
    if (Aeval->length == 0 || Beval->length == 0)
        goto choose_prime_outer;

    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Aeval, ectx));
    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Beval, ectx));

    fq_nmod_mpolyu_gcdp_zippel_univar_no_cofactors(Geval, Aeval, Beval, ectx);

    if (fq_nmod_mpolyu_is_one(Geval, ectx))
    {
        fq_nmod_mpolyu_one(G, ctx);
        fq_nmod_mpolyu_swap(Abar, A);
        fq_nmod_mpolyu_swap(Bbar, B);
        success = 1;
        goto finished;
    }

    FLINT_ASSERT(Geval->length > 0);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        if (Geval->exps[0] > H->exps[0])
        {
            goto choose_prime_outer;
        }
        else if (Geval->exps[0] < H->exps[0])
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
        }
    }

    fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff(Geval, ectx), ectx->fqctx);
    fq_nmod_mul(t, t, geval, ectx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, t, ectx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        changed = fq_nmod_mpolyun_interp_crt_lg_mpolyu(&lastdeg,
                                    H, Ht, modulus, ctx, Geval, ectx, cur_emb);
        fq_nmod_poly_mul(modulus, modulus, cur_emb->h, ctx->fqctx);

        have_enough = fq_nmod_poly_degree(modulus, ctx->fqctx) >= bound;

        if (changed && !have_enough)
        {
            goto choose_prime_outer;
        }

        if (!changed || have_enough)
        {
            fq_nmod_mpolyun_content_poly(hc, H, ctx);
            fq_nmod_mpolyun_divexact_poly(Ht, H, hc, ctx);
            fq_nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
            if (    fq_nmod_mpolyuu_divides(Abar, A, G, 1, ctx)
                 && fq_nmod_mpolyuu_divides(Bbar, B, G, 1, ctx))
            {
                success = 1;
                goto finished;
            }
        }

        if (have_enough)
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
            goto choose_prime_outer;
        }
    }
    else
    {
        fq_nmod_mpolyun_interp_lift_lg_mpolyu(H, ctx, Geval, ectx, cur_emb);
        fq_nmod_poly_set(modulus, cur_emb->h, ctx->fqctx);
    }

    goto choose_prime_outer;

finished:

    fq_nmod_poly_clear(gamma, ctx->fqctx);
    fq_nmod_poly_clear(hc, ctx->fqctx);
    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_mpolyun_clear(An, ctx);
    fq_nmod_mpolyun_clear(Bn, ctx);
    fq_nmod_mpolyun_clear(H, ctx);
    fq_nmod_mpolyun_clear(Ht, ctx);

    fq_nmod_mpolyu_clear(Aeval, ectx);
    fq_nmod_mpolyu_clear(Beval, ectx);
    fq_nmod_mpolyu_clear(Geval, ectx);
    fq_nmod_clear(geval, ectx->fqctx);
    fq_nmod_clear(t, ectx->fqctx);

    bad_fq_nmod_mpoly_embed_chooser_clear(embc, ectx, ctx, randstate);

    return success;
}


int fq_nmod_mpolyu_gcdm_zippel(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t Abar,
    fq_nmod_mpolyu_t Bbar,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    slong degbound;
    slong bound;
    slong Alastdeg, Blastdeg;
    slong lastdeg;
    int success, changed, have_enough;
    fq_nmod_mpolyun_t An, Bn, Hn, Ht;
    fq_nmod_poly_t modulus, gamma, hc;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx;
    fq_nmod_mpolyu_t Aeval, Beval, Geval, Abareval, Bbareval, Gform;
    fq_nmod_t t, gammaeval;

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(B->bits == A->bits);

    success = fq_nmod_mpolyu_gcdp_zippel(G, Abar, Bbar, A, B,
                                 ctx->minfo->nvars - 1, ctx, zinfo, randstate);
    if (success)
    {
        return 1;
    }

    /* bivariate more comfortable separated */
    if (ctx->minfo->nvars == 1)
    {
        return fq_nmod_mpolyu_gcdm_zippel_bivar(G, Abar, Bbar, A, B,
                                                        ctx, zinfo, randstate);
    }

    FLINT_ASSERT(ctx->minfo->nvars > 1);

    fq_nmod_poly_init(hc, ctx->fqctx);
    fq_nmod_poly_init(modulus, ctx->fqctx);

    fq_nmod_mpolyun_init(An, A->bits, ctx);
    fq_nmod_mpolyun_init(Bn, A->bits, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(An, A, ctx->minfo->nvars - 1, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(Bn, B, ctx->minfo->nvars - 1, ctx);

    FLINT_ASSERT(An->bits == A->bits);
    FLINT_ASSERT(Bn->bits == A->bits);
    FLINT_ASSERT(An->length > 0);
    FLINT_ASSERT(Bn->length > 0);
    FLINT_ASSERT(An->exps[A->length - 1] == 0);
    FLINT_ASSERT(Bn->exps[B->length - 1] == 0);

    fq_nmod_poly_init(gamma, ctx->fqctx);
    fq_nmod_poly_gcd(gamma, fq_nmod_mpolyun_leadcoeff_poly(An, ctx),
                            fq_nmod_mpolyun_leadcoeff_poly(Bn, ctx), ctx->fqctx);

    /* bound on the number of images */
    Alastdeg = fq_nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = fq_nmod_mpolyun_lastdeg(Bn, ctx);
    bound = 1 + fq_nmod_poly_degree(gamma, ctx->fqctx)
              + FLINT_MIN(Alastdeg, Blastdeg);

    /* degree bound on the gcd */
    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    fq_nmod_poly_one(modulus, ctx->fqctx);

    fq_nmod_mpolyun_init(Hn, A->bits, ctx);
    fq_nmod_mpolyun_init(Ht, A->bits, ctx);

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, ectx, ctx, randstate);

    /*
        Once Aeval, Beval, ..., t are inited in ectx->fqctx, they do not need
        to be cleared and reinited when ectx->fqctx changes.
    */
    fq_nmod_mpolyu_init(Aeval, A->bits, ectx);
    fq_nmod_mpolyu_init(Beval, A->bits, ectx);
    fq_nmod_mpolyu_init(Geval, A->bits, ectx);
    fq_nmod_mpolyu_init(Abareval, A->bits, ectx);
    fq_nmod_mpolyu_init(Bbareval, A->bits, ectx);
    fq_nmod_mpolyu_init(Gform, A->bits, ectx);
    fq_nmod_init(gammaeval, ectx->fqctx);
    fq_nmod_init(t, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime_outer:

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, randstate);
    if (cur_emb == NULL)
    {
        /* ran out of primes */
        success = 0;
        goto finished;
    }

have_prime:

    FLINT_ASSERT(cur_emb != NULL);

    /* make sure reduction does not kill both lc(A) and lc(B) */
    bad_fq_nmod_embed_sm_to_lg(gammaeval, gamma, cur_emb);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
    {
        goto choose_prime_outer;
    }

    /* make sure reduction does not kill either A or B */
    fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Aeval, An, ectx, ctx, cur_emb);
    fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Beval, Bn, ectx, ctx, cur_emb);
    if (Aeval->length == 0 || Beval->length == 0)
        goto choose_prime_outer;

    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Aeval, ectx));
    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Beval, ectx));

    success = fq_nmod_mpolyu_gcdp_zippel(Geval, Abareval, Bbareval, Aeval, Beval,
                                ctx->minfo->nvars - 2, ectx, zinfo, randstate);

    FLINT_ASSERT(!success || Geval->length > 0);
    if (!success || Geval->exps[0] > degbound)
        goto choose_prime_outer;

    degbound = Geval->exps[0];

    if (Geval->length == 1 && Geval->exps[0] == 0)
    {
        FLINT_ASSERT(fq_nmod_mpoly_is_one(Geval->coeffs + 0, ectx));
        FLINT_ASSERT(!fq_nmod_mpoly_is_zero(Geval->coeffs + 0, ectx));
        fq_nmod_mpolyu_one(G, ctx);
        fq_nmod_mpolyu_swap(Abar, A);
        fq_nmod_mpolyu_swap(Bbar, B);
        success = 1;
        goto finished;
    }

    fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff(Geval, ectx), ectx->fqctx);
    fq_nmod_mul(t, t, gammaeval, ectx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, t, ectx);

    fq_nmod_mpolyu_setform(Gform, Geval, ectx);
    fq_nmod_mpolyun_interp_lift_lg_mpolyu(Hn, ctx, Geval, ectx, cur_emb);

    fq_nmod_poly_set(modulus, cur_emb->h, ctx->fqctx);

choose_prime_inner:

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, randstate);
    if (cur_emb == NULL)
    {
        success = 0;
        goto finished;
    }

    /* make sure reduction does not kill both lc(A) and lc(B) */
    bad_fq_nmod_embed_sm_to_lg(gammaeval, gamma, cur_emb);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
        goto choose_prime_inner;

    /* make sure reduction does not kill either A or B */
    fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Aeval, An, ectx, ctx, cur_emb);
    fq_nmod_mpolyun_interp_reduce_lg_mpolyu(Beval, Bn, ectx, ctx, cur_emb);
    if (Aeval->length == 0 || Beval->length == 0)
        goto choose_prime_inner;

    switch (fq_nmod_mpolyu_gcds_zippel(Geval, Aeval, Beval, Gform,
                           ctx->minfo->nvars - 1, ectx, randstate, &degbound))
    {
        default:
            FLINT_ASSERT(0);
        case nmod_gcds_form_main_degree_too_high:
        case nmod_gcds_form_wrong:
        case nmod_gcds_no_solution:
            goto choose_prime_outer;
        case nmod_gcds_scales_not_found:
        case nmod_gcds_eval_point_not_found:
        case nmod_gcds_eval_gcd_deg_too_high:
            goto choose_prime_inner;
        case nmod_gcds_success:
            break;
    }

    if (fq_nmod_is_zero(fq_nmod_mpolyu_leadcoeff(Geval, ectx), ectx->fqctx))
        goto choose_prime_inner;

    fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff(Geval, ectx), ectx->fqctx);
    fq_nmod_mul(t, t, gammaeval, ectx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, t, ectx);

    changed = fq_nmod_mpolyun_interp_mcrt_lg_mpolyu(&lastdeg, Hn, ctx, modulus, Geval, ectx, cur_emb);
    fq_nmod_poly_mul(modulus, modulus, cur_emb->h, ctx->fqctx);

    have_enough = fq_nmod_poly_degree(modulus, ctx->fqctx) >= bound;

    if (changed && !have_enough)
    {
        goto choose_prime_inner;
    }

    if (!changed || have_enough)
    {
        fq_nmod_mpolyun_content_poly(hc, Hn, ctx);
        fq_nmod_mpolyun_divexact_poly(Ht, Hn, hc, ctx);
        fq_nmod_mpolyu_cvtfrom_mpolyun(G, Ht, ctx->minfo->nvars - 1, ctx);
        if (    fq_nmod_mpolyuu_divides(Abar, A, G, 1, ctx)
             && fq_nmod_mpolyuu_divides(Bbar, B, G, 1, ctx))
        {
            success = 1;
            goto finished;
        }
    }

    if (have_enough)
    {
        fq_nmod_poly_one(modulus, ctx->fqctx);
        goto choose_prime_outer;
    }

    goto choose_prime_inner;

finished:

    fq_nmod_poly_clear(gamma, ctx->fqctx);

    fq_nmod_poly_clear(hc, ctx->fqctx);
    fq_nmod_poly_clear(modulus, ctx->fqctx);

    fq_nmod_mpolyun_clear(An, ctx);
    fq_nmod_mpolyun_clear(Bn, ctx);
    fq_nmod_mpolyun_clear(Hn, ctx);
    fq_nmod_mpolyun_clear(Ht, ctx);

    fq_nmod_mpolyu_clear(Aeval, ectx);
    fq_nmod_mpolyu_clear(Beval, ectx);
    fq_nmod_mpolyu_clear(Geval, ectx);
    fq_nmod_mpolyu_clear(Abareval, ectx);
    fq_nmod_mpolyu_clear(Bbareval, ectx);
    fq_nmod_mpolyu_clear(Gform, ectx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(t, ectx->fqctx);

    bad_fq_nmod_mpoly_embed_chooser_clear(embc, ectx, ctx, randstate);

    return success;
}


int fq_nmod_mpoly_gcd_zippel(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                        const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t wbits;
    flint_rand_t randstate;
    int success = 0;
    mpoly_zipinfo_t zinfo;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    fq_nmod_mpoly_t Ac, Bc, Gc;
    ulong * shift, * stride;
    ulong amin, bmin;

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        if (fq_nmod_mpoly_is_zero(B, ctx))
        {
            fq_nmod_mpoly_zero(G, ctx);
        }
        else
        {
            fq_nmod_mpoly_make_monic(G, B, ctx);
        }
        return 1;
    }

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
        return 0;

    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    if (ctx->minfo->nvars == 1)
    {
        fq_nmod_poly_t a, b, g;
        fq_nmod_poly_init(a, ctx->fqctx);
        fq_nmod_poly_init(b, ctx->fqctx);
        fq_nmod_poly_init(g, ctx->fqctx);
        _fq_nmod_mpoly_to_fq_nmod_poly_deflate(a, A, 0, shift, stride, ctx);
        _fq_nmod_mpoly_to_fq_nmod_poly_deflate(b, B, 0, shift, stride, ctx);
        fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
        _fq_nmod_mpoly_from_fq_nmod_poly_inflate(G, A->bits, g, 0, shift, stride, ctx);
        fq_nmod_poly_clear(a, ctx->fqctx);
        fq_nmod_poly_clear(b, ctx->fqctx);
        fq_nmod_poly_clear(g, ctx->fqctx);

        success = 1;
        goto cleanup1;
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(!fq_nmod_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fq_nmod_mpoly_is_zero(B, ctx));
    FLINT_ASSERT(ctx->minfo->nvars > WORD(1));

    flint_randinit(randstate);

    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    fq_nmod_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    fq_nmod_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
        zinfo->perm[i] = i;

    wbits = FLINT_MAX(A->bits, B->bits);

    fq_nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->fqctx);
    fq_nmod_mpolyu_init(Au, wbits, uctx);
    fq_nmod_mpolyu_init(Bu, wbits, uctx);
    fq_nmod_mpolyu_init(Gu, wbits, uctx);
    fq_nmod_mpolyu_init(Abaru, wbits, uctx);
    fq_nmod_mpolyu_init(Bbaru, wbits, uctx);
    fq_nmod_mpoly_init3(Ac, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Bc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Gc, 0, wbits, uctx);

    fq_nmod_mpoly_to_mpolyu_perm_deflate(Au, uctx, A, ctx, zinfo->perm,
                                                                shift, stride);
    fq_nmod_mpoly_to_mpolyu_perm_deflate(Bu, uctx, B, ctx, zinfo->perm,
                                                                shift, stride);

    FLINT_ASSERT(Au->length > 0);
    FLINT_ASSERT(Bu->length > 0);
    amin = Au->exps[Au->length - 1];
    bmin = Bu->exps[Bu->length - 1];
    fq_nmod_mpolyu_shift_right(Au, amin);
    fq_nmod_mpolyu_shift_right(Bu, bmin);

    success = fq_nmod_mpolyu_content_mpoly(Ac, Au, uctx);
    success = success && fq_nmod_mpolyu_content_mpoly(Bc, Bu, uctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpolyu_divexact_mpoly_inplace(Au, Ac, uctx);
    fq_nmod_mpolyu_divexact_mpoly_inplace(Bu, Bc, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */
    success = fq_nmod_mpolyu_gcdm_zippel(Gu, Abaru, Bbaru, Au, Bu,
                                                       uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    success = _fq_nmod_mpoly_gcd(Gc, wbits, Ac, Bc, uctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);
    fq_nmod_mpolyu_shift_left(Gu, FLINT_MIN(amin, bmin));

    fq_nmod_mpoly_from_mpolyu_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                         Gu, uctx, zinfo->perm, shift, stride);
    success = 1;

    fq_nmod_mpoly_make_monic(G, G, ctx);

cleanup:

    fq_nmod_mpolyu_clear(Au, uctx);
    fq_nmod_mpolyu_clear(Bu, uctx);
    fq_nmod_mpolyu_clear(Gu, uctx);
    fq_nmod_mpolyu_clear(Abaru, uctx);
    fq_nmod_mpolyu_clear(Bbaru, uctx);
    fq_nmod_mpoly_clear(Ac, uctx);
    fq_nmod_mpoly_clear(Bc, uctx);
    fq_nmod_mpoly_clear(Gc, uctx);

    fq_nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

cleanup1:

    flint_free(shift);
    flint_free(stride);

    return success;
}
