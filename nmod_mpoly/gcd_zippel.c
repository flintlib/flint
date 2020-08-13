/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"

int nmod_mpolyu_gcdm_zippel_bivar(
    nmod_mpolyu_t G,
    nmod_mpolyu_t Abar,
    nmod_mpolyu_t Bbar,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    slong var = 0;
    slong Alastdeg, Blastdeg;
    slong bound;
    slong lastdeg;
    int success = 0, changed, have_enough;
    nmod_mpolyun_t An, Bn, H, Ht;
    slong deg;
    fq_nmod_mpoly_ctx_t ffctx;
    fq_nmod_mpolyu_t Aeval, Beval, Geval;
    nmod_poly_t modulus, gamma, hc;
    fq_nmod_t t, geval;

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(Abar->bits == A->bits);
    FLINT_ASSERT(Bbar->bits == A->bits);
    FLINT_ASSERT(B->bits == A->bits);

    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);
    nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    FLINT_ASSERT(An->bits == B->bits);
    FLINT_ASSERT(An->bits == G->bits);
    FLINT_ASSERT(An->length > 0);
    FLINT_ASSERT(Bn->length > 0);
    FLINT_ASSERT(An->exps[A->length - 1] == 0);
    FLINT_ASSERT(Bn->exps[B->length - 1] == 0);

    nmod_poly_init_mod(gamma, ctx->ffinfo->mod);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_poly(An, ctx),
                         nmod_mpolyun_leadcoeff_poly(Bn, ctx));

    /* bound on the number of images */
    Alastdeg = nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = nmod_mpolyun_lastdeg(Bn, ctx);
    bound = 1 + nmod_poly_degree(gamma) + FLINT_MIN(Alastdeg, Blastdeg);

    nmod_poly_init_mod(hc, ctx->ffinfo->mod);
    nmod_poly_init_mod(modulus, ctx->ffinfo->mod);
    nmod_poly_one(modulus);
    nmod_mpolyun_init(H, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    deg = WORD(20)/(FLINT_BIT_COUNT(ctx->ffinfo->mod.n));
    deg = FLINT_MAX(WORD(2), deg);

    fq_nmod_mpoly_ctx_init_deg(ffctx, ctx->minfo->nvars, ORD_LEX,
                                                      ctx->ffinfo->mod.n, deg);

    fq_nmod_mpolyu_init(Aeval, A->bits, ffctx);
    fq_nmod_mpolyu_init(Beval, A->bits, ffctx);
    fq_nmod_mpolyu_init(Geval, A->bits, ffctx);
    fq_nmod_init(geval, ffctx->fqctx);
    fq_nmod_init(t, ffctx->fqctx);

    while (1)
    {
        /* TODO:
            instead of simply increasing the degree everytime, try to find more
            irreducibles of the same degree
        */
        deg++;
        if (deg > 10000)
        {
            /* ran out of primes */
            success = 0;
            goto finished;
        }

        fq_nmod_mpolyu_clear(Aeval, ffctx);
        fq_nmod_mpolyu_clear(Beval, ffctx);
        fq_nmod_mpolyu_clear(Geval, ffctx);
        fq_nmod_clear(geval, ffctx->fqctx);
        fq_nmod_clear(t, ffctx->fqctx);

        fq_nmod_mpoly_ctx_change_modulus(ffctx, deg);

        fq_nmod_mpolyu_init(Aeval, A->bits, ffctx);
        fq_nmod_mpolyu_init(Beval, A->bits, ffctx);
        fq_nmod_mpolyu_init(Geval, A->bits, ffctx);
        fq_nmod_init(geval, ffctx->fqctx);
        fq_nmod_init(t, ffctx->fqctx);

        /* make sure reduction does not kill both lc(A) and lc(B) */
        nmod_poly_rem(geval, gamma, ffctx->fqctx->modulus);
        if (fq_nmod_is_zero(geval, ffctx->fqctx))
            goto outer_continue;

        /* make sure reduction does not kill either A or B */
        nmod_mpolyun_interp_reduce_lg_mpolyu(Aeval, An, ffctx, ctx);
        nmod_mpolyun_interp_reduce_lg_mpolyu(Beval, Bn, ffctx, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
            goto outer_continue;

        FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Aeval, ffctx));
        FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Beval, ffctx));

        fq_nmod_mpolyu_gcdp_zippel_univar_no_cofactors(Geval, Aeval, Beval, ffctx);

        if (fq_nmod_mpolyu_is_one(Geval, ffctx))
        {
            nmod_mpolyu_one(G, ctx);
            nmod_mpolyu_swap(Abar, A, ctx);
            nmod_mpolyu_swap(Bbar, B, ctx);
            success = 1;
            goto finished;
        }

        FLINT_ASSERT(Geval->length > 0);

        if (nmod_poly_degree(modulus) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;
            }
            else if (Geval->exps[0] < H->exps[0])
            {
                nmod_poly_one(modulus);                
            }
        }

        fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff(Geval, ffctx), ffctx->fqctx);
        fq_nmod_mul(t, t, geval, ffctx->fqctx);
        fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, t, ffctx);

        if (nmod_poly_degree(modulus) > 0)
        {
            changed = nmod_mpolyun_interp_crt_lg_mpolyu(&lastdeg, H, Ht,
                                                   modulus, ctx, Geval, ffctx);
            nmod_poly_mul(modulus, modulus, ffctx->fqctx->modulus);

            have_enough = nmod_poly_degree(modulus) >= bound;

            if (changed && !have_enough)
            {
                goto outer_continue;
            }

            if (!changed || have_enough)
            {
                nmod_mpolyun_content_last(hc, H, ctx);
                nmod_mpolyun_set(Ht, H, ctx);
                nmod_mpolyun_divexact_last(Ht, hc, ctx);
                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (   nmod_mpolyuu_divides(Abar, A, G, 1, ctx)
                    && nmod_mpolyuu_divides(Bbar, B, G, 1, ctx))
                {
                    success = 1;
                    goto finished;
                }
            }

            if (have_enough)
            {
                nmod_poly_one(modulus);
                goto outer_continue;
            }
        }
        else
        {
            nmod_mpolyun_interp_lift_lg_mpolyu(H, ctx, Geval, ffctx);
            nmod_poly_set(modulus, ffctx->fqctx->modulus);
        }

outer_continue:;
    }

finished:

    nmod_poly_clear(gamma);
    nmod_poly_clear(hc);
    nmod_poly_clear(modulus);
    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(H, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    fq_nmod_clear(geval, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);
    fq_nmod_mpolyu_clear(Aeval, ffctx);
    fq_nmod_mpolyu_clear(Beval, ffctx);
    fq_nmod_mpolyu_clear(Geval, ffctx);
    fq_nmod_mpoly_ctx_clear(ffctx);

    return success;
}


int nmod_mpolyu_gcdm_zippel(
    nmod_mpolyu_t G,
    nmod_mpolyu_t Abar,
    nmod_mpolyu_t Bbar,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    slong degbound;
    slong bound;
    slong Alastdeg, Blastdeg;
    slong lastdeg;
    int success, changed, have_enough;
    nmod_mpolyun_t An, Bn, Hn, Ht;
    slong deg;
    fq_nmod_mpoly_ctx_t ffctx;
    fq_nmod_mpolyu_t Aff, Bff, Gff, Abarff, Bbarff, Gform;
    nmod_poly_t modulus, gamma, hc;
    fq_nmod_t t, gammaff;

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(Abar->bits == A->bits);
    FLINT_ASSERT(Bbar->bits == A->bits);
    FLINT_ASSERT(B->bits == A->bits);

    success = nmod_mpolyu_gcdp_zippel(G, Abar, Bbar, A, B,
                                 ctx->minfo->nvars - 1, ctx, zinfo, randstate);
    if (success)
    {
        return 1;
    }

    /* bivariate more comfortable separated */
    if (ctx->minfo->nvars == 1)
    {
        return nmod_mpolyu_gcdm_zippel_bivar(G, Abar, Bbar, A, B,
                                                        ctx, zinfo, randstate);
    }

    FLINT_ASSERT(ctx->minfo->nvars > 1);

    nmod_poly_init(hc, ctx->ffinfo->mod.n);
    nmod_poly_init(modulus, ctx->ffinfo->mod.n);

    nmod_mpolyun_init(An, A->bits, ctx);
    nmod_mpolyun_init(Bn, A->bits, ctx);
    nmod_mpolyu_cvtto_mpolyun(An, A, ctx->minfo->nvars - 1, ctx);
    nmod_mpolyu_cvtto_mpolyun(Bn, B, ctx->minfo->nvars - 1, ctx);

    FLINT_ASSERT(An->bits == B->bits);
    FLINT_ASSERT(An->bits == G->bits);
    FLINT_ASSERT(An->length > 0);
    FLINT_ASSERT(Bn->length > 0);
    FLINT_ASSERT(An->exps[A->length - 1] == 0);
    FLINT_ASSERT(Bn->exps[B->length - 1] == 0);

    nmod_poly_init_mod(gamma, ctx->ffinfo->mod);
    nmod_poly_gcd(gamma, nmod_mpolyun_leadcoeff_poly(An, ctx),
                         nmod_mpolyun_leadcoeff_poly(Bn, ctx));

    /* bound on the number of images */
    Alastdeg = nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = nmod_mpolyun_lastdeg(Bn, ctx);
    bound = 1 + nmod_poly_degree(gamma) + FLINT_MIN(Alastdeg, Blastdeg);

    /* degree bound on the gcd */
    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    nmod_poly_one(modulus);

    nmod_mpolyun_init(Hn, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    deg = WORD(20)/(FLINT_BIT_COUNT(ctx->ffinfo->mod.n));
    deg = FLINT_MAX(WORD(2), deg);

    fq_nmod_mpoly_ctx_init_deg(ffctx, ctx->minfo->nvars, ORD_LEX,
                                                      ctx->ffinfo->mod.n, deg);

    fq_nmod_mpolyu_init(Aff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Abarff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bbarff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gform, A->bits, ffctx);
    fq_nmod_init(gammaff, ffctx->fqctx);
    fq_nmod_init(t, ffctx->fqctx);

choose_prime_outer:

    /* same TODO */
    deg++;
    if (deg > 10000)
    {
        /* ran out of primes */
        success = 0;
        goto finished;
    }

    fq_nmod_mpolyu_clear(Aff, ffctx);
    fq_nmod_mpolyu_clear(Bff, ffctx);
    fq_nmod_mpolyu_clear(Gff, ffctx);
    fq_nmod_mpolyu_clear(Abarff, ffctx);
    fq_nmod_mpolyu_clear(Bbarff, ffctx);
    fq_nmod_mpolyu_clear(Gform, ffctx);
    fq_nmod_clear(gammaff, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);

    fq_nmod_mpoly_ctx_change_modulus(ffctx, deg);

    fq_nmod_mpolyu_init(Aff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Abarff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bbarff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gform, A->bits, ffctx);
    fq_nmod_init(gammaff, ffctx->fqctx);
    fq_nmod_init(t, ffctx->fqctx);

    /* make sure reduction does not kill both lc(A) and lc(B) */
    nmod_poly_rem(gammaff, gamma, ffctx->fqctx->modulus);
    if (fq_nmod_is_zero(gammaff, ffctx->fqctx))
        goto choose_prime_outer;

    /* make sure reduction does not kill either A or B */
    nmod_mpolyun_interp_reduce_lg_mpolyu(Aff, An, ffctx, ctx);
    nmod_mpolyun_interp_reduce_lg_mpolyu(Bff, Bn, ffctx, ctx);
    if (Aff->length == 0 || Bff->length == 0)
        goto choose_prime_outer;

    success = fq_nmod_mpolyu_gcdp_zippel(Gff, Abarff, Bbarff, Aff, Bff,
                               ctx->minfo->nvars - 2, ffctx, zinfo, randstate);
    if (!success || Gff->exps[0] > degbound)
        goto choose_prime_outer;
    degbound = Gff->exps[0];

    if (Gff->length == 1 && Gff->exps[0] == 0)
    {
        FLINT_ASSERT(fq_nmod_mpoly_is_one(Gff->coeffs + 0, ffctx));
        FLINT_ASSERT(!fq_nmod_mpoly_is_zero(Gff->coeffs + 0, ffctx));
        nmod_mpolyu_one(G, ctx);
        nmod_mpolyu_swap(Abar, A, ctx);
        nmod_mpolyu_swap(Bbar, B, ctx);
        success = 1;
        goto finished;
    }

    fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff(Gff, ffctx), ffctx->fqctx);
    fq_nmod_mul(t, t, gammaff, ffctx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Gff, t, ffctx);

    fq_nmod_mpolyu_setform(Gform, Gff, ffctx);
    nmod_mpolyun_interp_lift_lg_mpolyu(Hn, ctx, Gff, ffctx);

    nmod_poly_set(modulus, ffctx->fqctx->modulus);

choose_prime_inner:

    deg++;
    if (deg > 1000)
    {
        /* ran out of primes */
        success = 0;
        goto finished;
    }

    fq_nmod_mpolyu_clear(Aff, ffctx);
    fq_nmod_mpolyu_clear(Bff, ffctx);
    fq_nmod_mpolyu_clear(Gff, ffctx);
    fq_nmod_mpolyu_clear(Abarff, ffctx);
    fq_nmod_mpolyu_clear(Bbarff, ffctx);
    fq_nmod_clear(gammaff, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);

    fq_nmod_mpoly_ctx_change_modulus(ffctx, deg);

    fq_nmod_mpolyu_init(Aff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Abarff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bbarff, A->bits, ffctx);
    fq_nmod_init(gammaff, ffctx->fqctx);
    fq_nmod_init(t, ffctx->fqctx);

    /* make sure reduction does not kill both lc(A) and lc(B) */
    nmod_poly_rem(gammaff, gamma, ffctx->fqctx->modulus);
    if (fq_nmod_is_zero(gammaff, ffctx->fqctx))
        goto choose_prime_inner;

    /* make sure reduction does not kill either A or B */
    nmod_mpolyun_interp_reduce_lg_mpolyu(Aff, An, ffctx, ctx);
    nmod_mpolyun_interp_reduce_lg_mpolyu(Bff, Bn, ffctx, ctx);
    if (Aff->length == 0 || Bff->length == 0)
        goto choose_prime_inner;

    switch (fq_nmod_mpolyu_gcds_zippel(Gff, Aff, Bff, Gform,
                           ctx->minfo->nvars - 1, ffctx, randstate, &degbound))
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

    if (fq_nmod_is_zero(fq_nmod_mpolyu_leadcoeff(Gff, ffctx), ffctx->fqctx))
        goto choose_prime_inner;

    fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff(Gff, ffctx), ffctx->fqctx);
    fq_nmod_mul(t, t, gammaff, ffctx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Gff, t, ffctx);

    changed = nmod_mpolyun_interp_mcrt_lg_mpolyu(&lastdeg, Hn, ctx, modulus, Gff, ffctx);
    nmod_poly_mul(modulus, modulus, ffctx->fqctx->modulus);

    have_enough = nmod_poly_degree(modulus) >= bound;

    if (changed && !have_enough)
    {
        goto choose_prime_inner;
    }

    if (!changed || have_enough)
    {
        nmod_mpolyun_content_last(hc, Hn, ctx);
        nmod_mpolyun_set(Ht, Hn, ctx);
        nmod_mpolyun_divexact_last(Ht, hc, ctx);
        nmod_mpolyu_cvtfrom_mpolyun(G, Ht, ctx->minfo->nvars - 1, ctx);
        if (   nmod_mpolyuu_divides(Abar, A, G, 1, ctx)
            && nmod_mpolyuu_divides(Bbar, B, G, 1, ctx))
        {
            success = 1;
            goto finished;
        }
    }

    if (have_enough)
    {
        nmod_poly_one(modulus);
        goto choose_prime_outer;
    }

    goto choose_prime_inner;

finished:

    nmod_poly_clear(gamma);

    nmod_poly_clear(hc);
    nmod_poly_clear(modulus);

    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(Hn, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    fq_nmod_mpolyu_clear(Aff, ffctx);
    fq_nmod_mpolyu_clear(Bff, ffctx);
    fq_nmod_mpolyu_clear(Gff, ffctx);
    fq_nmod_mpolyu_clear(Abarff, ffctx);
    fq_nmod_mpolyu_clear(Bbarff, ffctx);
    fq_nmod_mpolyu_clear(Gform, ffctx);
    fq_nmod_clear(gammaff, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);

    fq_nmod_mpoly_ctx_clear(ffctx);

    return success;
}


int nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t wbits;
    flint_rand_t randstate;
    int success = 0;
    mpoly_zipinfo_t zinfo;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    nmod_mpoly_t Ac, Bc, Gc;
    ulong * shift, * stride;
    ulong amin, bmin;

    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
        {
            nmod_mpoly_zero(G, ctx);
        }
        else
        {
            nmod_mpoly_make_monic(G, B, ctx);
        }
        return 1;
    }

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_make_monic(G, A, ctx);
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
        nmod_poly_t a, b, g;
        nmod_poly_init_mod(a, ctx->ffinfo->mod);
        nmod_poly_init_mod(b, ctx->ffinfo->mod);
        nmod_poly_init_mod(g, ctx->ffinfo->mod);
        _nmod_mpoly_to_nmod_poly_deflate(a, A, 0, shift, stride, ctx);
        _nmod_mpoly_to_nmod_poly_deflate(b, B, 0, shift, stride, ctx);
        nmod_poly_gcd(g, a, b);
        _nmod_mpoly_from_nmod_poly_inflate(G, A->bits, g, 0, shift, stride, ctx);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);

        success = 1;
        goto cleanup1;
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(!nmod_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!nmod_mpoly_is_zero(B, ctx));
    FLINT_ASSERT(ctx->minfo->nvars > WORD(1));

    flint_randinit(randstate);

    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    nmod_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    nmod_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
        zinfo->perm[i] = i;

    wbits = FLINT_MAX(A->bits, B->bits);

    nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->ffinfo->mod.n);

    nmod_mpolyu_init(Au, wbits, uctx);
    nmod_mpolyu_init(Bu, wbits, uctx);
    nmod_mpolyu_init(Gu, wbits, uctx);
    nmod_mpolyu_init(Abaru, wbits, uctx);
    nmod_mpolyu_init(Bbaru, wbits, uctx);
    nmod_mpoly_init3(Ac, 0, wbits, uctx);
    nmod_mpoly_init3(Bc, 0, wbits, uctx);
    nmod_mpoly_init3(Gc, 0, wbits, uctx);

    nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(Au, uctx, A, ctx, zinfo->perm,
                                                       shift, stride, NULL, 0);
    nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(Bu, uctx, B, ctx, zinfo->perm,
                                                       shift, stride, NULL, 0);

    FLINT_ASSERT(Au->length > 0);
    FLINT_ASSERT(Bu->length > 0);
    amin = Au->exps[Au->length - 1];
    bmin = Bu->exps[Bu->length - 1];
    nmod_mpolyu_shift_right(Au, amin);
    nmod_mpolyu_shift_right(Bu, bmin);

    success = nmod_mpolyu_content_mpoly_threaded_pool(Ac, Au, uctx, NULL, 0);
    success = success &&
              nmod_mpolyu_content_mpoly_threaded_pool(Bc, Bu, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    nmod_mpolyu_divexact_mpoly_inplace(Au, Ac, uctx);
    nmod_mpolyu_divexact_mpoly_inplace(Bu, Bc, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */
    success = nmod_mpolyu_gcdm_zippel(Gu, Abaru, Bbaru, Au, Bu,
                                                       uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    success = _nmod_mpoly_gcd_threaded_pool(Gc, wbits, Ac, Bc, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    nmod_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);
    nmod_mpolyu_shift_left(Gu, FLINT_MIN(amin, bmin));

    nmod_mpoly_from_mpolyu_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                         Gu, uctx, zinfo->perm, shift, stride);
    success = 1;

    nmod_mpoly_make_monic(G, G, ctx);

cleanup:

    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Gu, uctx);
    nmod_mpolyu_clear(Abaru, uctx);
    nmod_mpolyu_clear(Bbaru, uctx);
    nmod_mpoly_clear(Ac, uctx);
    nmod_mpoly_clear(Bc, uctx);
    nmod_mpoly_clear(Gc, uctx);

    nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

cleanup1:

    flint_free(shift);
    flint_free(stride);

    return success;
}
