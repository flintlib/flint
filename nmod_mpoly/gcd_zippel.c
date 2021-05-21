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
    n_poly_t modulus, gamma, hc;
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

    n_poly_init(gamma);
    n_poly_mod_gcd(gamma, nmod_mpolyun_leadcoeff_poly(An, ctx),
                          nmod_mpolyun_leadcoeff_poly(Bn, ctx), ctx->mod);

    /* bound on the number of images */
    Alastdeg = nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = nmod_mpolyun_lastdeg(Bn, ctx);
    bound = 1 + n_poly_degree(gamma) + FLINT_MIN(Alastdeg, Blastdeg);

    n_poly_init(hc);
    n_poly_init(modulus);
    n_poly_one(modulus);
    nmod_mpolyun_init(H, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    deg = WORD(20)/(FLINT_BIT_COUNT(ctx->mod.n));
    deg = FLINT_MAX(WORD(2), deg);

    fq_nmod_mpoly_ctx_init_deg(ffctx, ctx->minfo->nvars, ORD_LEX, ctx->mod.n, deg);

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
        n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(geval), gamma,
                   evil_const_cast_nmod_poly_to_n_poly(ffctx->fqctx->modulus),
                                                                     ctx->mod);
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

        if (n_poly_degree(modulus) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;
            }
            else if (Geval->exps[0] < H->exps[0])
            {
                n_poly_one(modulus);                
            }
        }

        n_fq_get_fq_nmod(t, fq_nmod_mpolyu_leadcoeff(Geval, ffctx), ffctx->fqctx);
        fq_nmod_inv(t, t, ffctx->fqctx);
        fq_nmod_mul(t, t, geval, ffctx->fqctx);
        fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, t, ffctx);

        if (n_poly_degree(modulus) > 0)
        {
            changed = nmod_mpolyun_interp_crt_lg_mpolyu(&lastdeg, H, Ht,
                                                   modulus, ctx, Geval, ffctx);
            n_poly_mod_mul(modulus, modulus,
                 evil_const_cast_nmod_poly_to_n_poly(ffctx->fqctx->modulus),
                                                                     ctx->mod);

            have_enough = n_poly_degree(modulus) >= bound;

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
                n_poly_one(modulus);
                goto outer_continue;
            }
        }
        else
        {
            nmod_mpolyun_interp_lift_lg_mpolyu(H, ctx, Geval, ffctx);
            n_poly_set_nmod_poly(modulus, ffctx->fqctx->modulus);
        }

outer_continue:;
    }

finished:

    n_poly_clear(gamma);
    n_poly_clear(hc);
    n_poly_clear(modulus);
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
    n_poly_t modulus, gamma, hc;
    fq_nmod_t t, gammaff;

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(Abar->bits == A->bits);
    FLINT_ASSERT(Bbar->bits == A->bits);
    FLINT_ASSERT(B->bits == A->bits);

    success = nmod_mpolyu_gcdp_zippel(G, Abar, Bbar, A, B,
                                        ctx->minfo->nvars - 1, ctx, randstate);
    if (success)
    {
        return 1;
    }

    /* bivariate more comfortable separated */
    if (ctx->minfo->nvars == 1)
    {
        return nmod_mpolyu_gcdm_zippel_bivar(G, Abar, Bbar, A, B, ctx, randstate);
    }

    FLINT_ASSERT(ctx->minfo->nvars > 1);

    n_poly_init(hc);
    n_poly_init(modulus);

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

    n_poly_init(gamma);
    n_poly_mod_gcd(gamma, nmod_mpolyun_leadcoeff_poly(An, ctx),
                          nmod_mpolyun_leadcoeff_poly(Bn, ctx), ctx->mod);

    /* bound on the number of images */
    Alastdeg = nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = nmod_mpolyun_lastdeg(Bn, ctx);
    bound = 1 + n_poly_degree(gamma) + FLINT_MIN(Alastdeg, Blastdeg);

    /* degree bound on the gcd */
    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    n_poly_one(modulus);

    nmod_mpolyun_init(Hn, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    deg = WORD(20)/(FLINT_BIT_COUNT(ctx->mod.n));
    deg = FLINT_MAX(WORD(2), deg);

    fq_nmod_mpoly_ctx_init_deg(ffctx, ctx->minfo->nvars, ORD_LEX, ctx->mod.n, deg);

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
    n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(gammaff), gamma,
         evil_const_cast_nmod_poly_to_n_poly(ffctx->fqctx->modulus), ctx->mod);
    if (fq_nmod_is_zero(gammaff, ffctx->fqctx))
        goto choose_prime_outer;

    /* make sure reduction does not kill either A or B */
    nmod_mpolyun_interp_reduce_lg_mpolyu(Aff, An, ffctx, ctx);
    nmod_mpolyun_interp_reduce_lg_mpolyu(Bff, Bn, ffctx, ctx);
    if (Aff->length == 0 || Bff->length == 0)
        goto choose_prime_outer;

    success = fq_nmod_mpolyu_gcdp_zippel(Gff, Abarff, Bbarff, Aff, Bff,
                                      ctx->minfo->nvars - 2, ffctx, randstate);
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

    n_fq_get_fq_nmod(t, fq_nmod_mpolyu_leadcoeff(Gff, ffctx), ffctx->fqctx);
    fq_nmod_inv(t, t, ffctx->fqctx);
    fq_nmod_mul(t, t, gammaff, ffctx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Gff, t, ffctx);

    fq_nmod_mpolyu_setform(Gform, Gff, ffctx);
    nmod_mpolyun_interp_lift_lg_mpolyu(Hn, ctx, Gff, ffctx);

    n_poly_set_nmod_poly(modulus, ffctx->fqctx->modulus);

choose_prime_inner:

    deg++;
    if (deg > 1000)
    {
        /* ran out of primes */
        success = 0;
        goto finished;
    }

    /* Gform need not be reinited because only its exponents are used */

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
    n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(gammaff), gamma,
         evil_const_cast_nmod_poly_to_n_poly(ffctx->fqctx->modulus), ctx->mod);
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

    n_fq_get_fq_nmod(t, fq_nmod_mpolyu_leadcoeff(Gff, ffctx), ffctx->fqctx);

    if (fq_nmod_is_zero(t, ffctx->fqctx))
        goto choose_prime_inner;

    fq_nmod_inv(t, t, ffctx->fqctx);
    fq_nmod_mul(t, t, gammaff, ffctx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Gff, t, ffctx);

    changed = nmod_mpolyun_interp_mcrt_lg_mpolyu(&lastdeg, Hn, ctx, modulus, Gff, ffctx);
    n_poly_mod_mul(modulus, modulus,
         evil_const_cast_nmod_poly_to_n_poly(ffctx->fqctx->modulus), ctx->mod);

    have_enough = n_poly_degree(modulus) >= bound;

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
        n_poly_one(modulus);
        goto choose_prime_outer;
    }

    goto choose_prime_inner;

finished:

    n_poly_clear(gamma);
    n_poly_clear(hc);
    n_poly_clear(modulus);

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


/* should find its way back here in interesting cases */
int nmod_mpoly_gcd_zippel(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    if (nmod_mpoly_is_zero(A, ctx) || nmod_mpoly_is_zero(B, ctx))
        return nmod_mpoly_gcd(G, A, B, ctx);

    return _nmod_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_ZIPPEL);
}

