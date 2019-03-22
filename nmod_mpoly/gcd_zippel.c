/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"

int nmod_mpolyu_gcdm_zippel_bivar(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    slong var = 0;
    slong Alastdeg;
    slong Blastdeg;
    slong bound;
    slong lastdeg;
    int success = 0, changed, have_enough;
    nmod_mpolyun_t An, Bn, H, Ht;
    slong deg;
    fq_nmod_mpoly_ctx_t ffctx;
    fq_nmod_mpolyu_t Aeval, Beval, Geval;
    nmod_polydr_t modulus, a,b,c,g;
    fq_nmod_t t, geval;

    FLINT_ASSERT(G->bits == A->bits);
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

    nmod_polydr_init(a, ctx->ffinfo);
    nmod_polydr_init(b, ctx->ffinfo);
    nmod_polydr_init(c, ctx->ffinfo);
    nmod_polydr_init(g, ctx->ffinfo);
    nmod_mpolyun_content_last(a, An, ctx);
    nmod_mpolyun_content_last(b, Bn, ctx);
    nmod_mpolyun_divexact_last(An, a, ctx);
    nmod_mpolyun_divexact_last(Bn, b, ctx);
    nmod_polydr_gcd(c, a, b, ctx->ffinfo);
    nmod_polydr_gcd(g, nmod_mpolyun_leadcoeff_ref(An, ctx),
                       nmod_mpolyun_leadcoeff_ref(Bn, ctx), ctx->ffinfo);

    Alastdeg = nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = nmod_mpolyun_lastdeg(Bn, ctx);

    /* bound on the number of images */
    bound = 1 + FLINT_MIN(Alastdeg, Blastdeg)
              + nmod_polydr_degree(g, ctx->ffinfo);

    nmod_polydr_init(modulus, ctx->ffinfo);
    nmod_polydr_one(modulus, ctx->ffinfo);
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
        nmod_polydr_rem(geval, g, ffctx->fqctx->modulus, ctx->ffinfo);
        if (fq_nmod_is_zero(geval, ffctx->fqctx))
            goto outer_continue;

        /* make sure reduction does not kill either A or B */
        nmod_mpolyun_redto_fq_nmod_mpolyu(Aeval, An, ffctx, ctx);
        nmod_mpolyun_redto_fq_nmod_mpolyu(Beval, Bn, ffctx, ctx);
        if (Aeval->length == 0 || Beval->length == 0)
            goto outer_continue;

        FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Aeval, ffctx));
        FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Beval, ffctx));

        fq_nmod_mpolyu_gcdp_zippel_univar(Geval, Aeval, Beval, ffctx);

        if (fq_nmod_mpolyu_is_one(Geval, ffctx))
        {
            nmod_mpolyu_cvtfrom_poly_notmain(G, c, var, ctx);
            success = 1;
            goto finished;
        }

        FLINT_ASSERT(Geval->length > 0);

        if (nmod_polydr_degree(modulus, ctx->ffinfo) > 0)
        {
            if (Geval->exps[0] > H->exps[0])
            {
                goto outer_continue;
            }
            else if (Geval->exps[0] < H->exps[0])
            {
                nmod_polydr_one(modulus, ctx->ffinfo);                
            }
        }

        fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff_ref(Geval, ffctx), ffctx->fqctx);
        fq_nmod_mul(t, t, geval, ffctx->fqctx);
        fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, t, ffctx);

        if (nmod_polydr_degree(modulus, ctx->ffinfo) > 0)
        {
            changed = nmod_mpolyun_addinterp_fq_nmod_mpolyu(&lastdeg, H, Ht,
                                                   modulus, ctx, Geval, ffctx);
            nmod_polydr_mul(modulus, modulus, ffctx->fqctx->modulus, ctx->ffinfo);

            have_enough = nmod_polydr_degree(modulus, ctx->ffinfo) >= bound;

            if (changed && !have_enough)
            {
                goto outer_continue;
            }

            if (!changed || have_enough)
            {
                nmod_mpolyun_content_last(a, H, ctx);
                nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                nmod_mpolyun_divexact_last(Ht, a, ctx);
                nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);
                if (   nmod_mpolyu_divides(A, G, ctx)
                    && nmod_mpolyu_divides(B, G, ctx))
                {
                    success = 1;
                    goto finished;
                }
            }

            if (have_enough)
            {
                nmod_polydr_one(modulus, ctx->ffinfo);
                goto outer_continue;
            }
        }
        else
        {
            nmod_mpolyun_set_fq_nmod_mpolyu(H, ctx, Geval, ffctx);
            nmod_polydr_set(modulus, ffctx->fqctx->modulus, ctx->ffinfo);
        }

outer_continue:
        NULL;
    }

finished:

    nmod_polydr_clear(a, ctx->ffinfo);
    nmod_polydr_clear(b, ctx->ffinfo);
    nmod_polydr_clear(c, ctx->ffinfo);
    nmod_polydr_clear(g, ctx->ffinfo);
    nmod_polydr_clear(modulus, ctx->ffinfo);
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
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    slong degbound;
    slong bound;
    slong lastdeg;
    int success, changed, have_enough;
    nmod_mpolyun_t An, Bn, Hn, Ht;
    slong deg;
    fq_nmod_mpoly_ctx_t ffctx;
    fq_nmod_mpolyu_t Aff, Bff, Gff, Gform;
    nmod_polydr_t modulus, gamma, hc;
    fq_nmod_t t, gammaff;

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(B->bits == A->bits);

    success = nmod_mpolyu_gcdp_zippel(G, A, B, ctx->minfo->nvars - 1, ctx, zinfo, randstate);
    if (success)
    {
        return 1;
    }

    /* bivariate more comfortable separated */
    if (ctx->minfo->nvars == 1)
    {
        return nmod_mpolyu_gcdm_zippel_bivar(G, A, B, ctx, zinfo, randstate);
    }

    FLINT_ASSERT(ctx->minfo->nvars > 1);

    nmod_polydr_init(hc, ctx->ffinfo);
    nmod_polydr_init(modulus, ctx->ffinfo);

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

    nmod_polydr_init(gamma, ctx->ffinfo);
    nmod_polydr_gcd(gamma, nmod_mpolyun_leadcoeff_ref(An, ctx),
                           nmod_mpolyun_leadcoeff_ref(Bn, ctx), ctx->ffinfo);

    /* bound on the number of images */
    bound = 1 + nmod_polydr_degree(gamma, ctx->ffinfo)
              + FLINT_MIN(nmod_mpolyun_lastdeg(An, ctx),
                          nmod_mpolyun_lastdeg(Bn, ctx));

    /* degree bound on the gcd */
    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    nmod_polydr_one(modulus, ctx->ffinfo);

    nmod_mpolyun_init(Hn, A->bits, ctx);
    nmod_mpolyun_init(Ht, A->bits, ctx);

    deg = WORD(20)/(FLINT_BIT_COUNT(ctx->ffinfo->mod.n));
    deg = FLINT_MAX(WORD(2), deg);

    fq_nmod_mpoly_ctx_init_deg(ffctx, ctx->minfo->nvars, ORD_LEX,
                                                      ctx->ffinfo->mod.n, deg);

    fq_nmod_mpolyu_init(Aff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gff, A->bits, ffctx);
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
    fq_nmod_mpolyu_clear(Gform, ffctx);
    fq_nmod_clear(gammaff, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);

    fq_nmod_mpoly_ctx_change_modulus(ffctx, deg);

    fq_nmod_mpolyu_init(Aff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gform, A->bits, ffctx);
    fq_nmod_init(gammaff, ffctx->fqctx);
    fq_nmod_init(t, ffctx->fqctx);

    /* make sure reduction does not kill both lc(A) and lc(B) */
    nmod_polydr_rem(gammaff, gamma, ffctx->fqctx->modulus, ctx->ffinfo);
    if (fq_nmod_is_zero(gammaff, ffctx->fqctx))
        goto choose_prime_outer;

    /* make sure reduction does not kill either A or B */
    nmod_mpolyun_redto_fq_nmod_mpolyu(Aff, An, ffctx, ctx);
    nmod_mpolyun_redto_fq_nmod_mpolyu(Bff, Bn, ffctx, ctx);
    if (Aff->length == 0 || Bff->length == 0)
        goto choose_prime_outer;

    success = fq_nmod_mpolyu_gcdp_zippel(Gff, Aff, Bff, ctx->minfo->nvars - 2,
                                                      ffctx, zinfo, randstate);
    if (!success || Gff->exps[0] > degbound)
        goto choose_prime_outer;
    degbound = Gff->exps[0];

    if (Gff->length == 1 && Gff->exps[0] == 0)
    {
        FLINT_ASSERT(fq_nmod_mpoly_is_one(Gff->coeffs + 0, ffctx));
        FLINT_ASSERT(!fq_nmod_mpoly_is_zero(Gff->coeffs + 0, ffctx));
        nmod_mpolyu_one(G, ctx);
        
        success = 1;
        goto finished;
    }

    fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff_ref(Gff, ffctx), ffctx->fqctx);
    fq_nmod_mul(t, t, gammaff, ffctx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Gff, t, ffctx);

    fq_nmod_mpolyu_setform(Gform, Gff, ffctx);
    nmod_mpolyun_set_fq_nmod_mpolyu(Hn, ctx, Gff, ffctx);

    nmod_polydr_set(modulus, ffctx->fqctx->modulus, ctx->ffinfo);

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
    fq_nmod_clear(gammaff, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);

    fq_nmod_mpoly_ctx_change_modulus(ffctx, deg);

    fq_nmod_mpolyu_init(Aff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Bff, A->bits, ffctx);
    fq_nmod_mpolyu_init(Gff, A->bits, ffctx);
    fq_nmod_init(gammaff, ffctx->fqctx);
    fq_nmod_init(t, ffctx->fqctx);

    /* make sure reduction does not kill both lc(A) and lc(B) */
    nmod_polydr_rem(gammaff, gamma, ffctx->fqctx->modulus, ctx->ffinfo);
    if (fq_nmod_is_zero(gammaff, ffctx->fqctx))
        goto choose_prime_inner;

    /* make sure reduction does not kill either A or B */
    nmod_mpolyun_redto_fq_nmod_mpolyu(Aff, An, ffctx, ctx);
    nmod_mpolyun_redto_fq_nmod_mpolyu(Bff, Bn, ffctx, ctx);
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
            NULL;
    }

    if (fq_nmod_is_zero(fq_nmod_mpolyu_leadcoeff_ref(Gff, ffctx), ffctx->fqctx))
        goto choose_prime_inner;

    fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff_ref(Gff, ffctx), ffctx->fqctx);
    fq_nmod_mul(t, t, gammaff, ffctx->fqctx);
    fq_nmod_mpolyu_scalar_mul_fq_nmod(Gff, t, ffctx);

    changed = nmod_mpolyun_CRT_fq_nmod_mpolyu(&lastdeg, Hn, ctx, modulus, Gff, ffctx);
    nmod_polydr_mul(modulus, modulus, ffctx->fqctx->modulus, ctx->ffinfo);

    have_enough = nmod_polydr_degree(modulus, ctx->ffinfo) >= bound;

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
        if (   nmod_mpolyu_divides(A, G, ctx)
            && nmod_mpolyu_divides(B, G, ctx))
        {
            success = 1;
            goto finished;
        }
    }

    if (have_enough)
    {
        nmod_polydr_one(modulus, ctx->ffinfo);
        goto choose_prime_outer;
    }

    goto choose_prime_inner;

finished:

    nmod_polydr_clear(gamma, ctx->ffinfo);

    nmod_polydr_clear(hc, ctx->ffinfo);
    nmod_polydr_clear(modulus, ctx->ffinfo);

    nmod_mpolyun_clear(An, ctx);
    nmod_mpolyun_clear(Bn, ctx);
    nmod_mpolyun_clear(Hn, ctx);
    nmod_mpolyun_clear(Ht, ctx);

    fq_nmod_mpolyu_clear(Aff, ffctx);
    fq_nmod_mpolyu_clear(Bff, ffctx);
    fq_nmod_mpolyu_clear(Gff, ffctx);
    fq_nmod_mpolyu_clear(Gform, ffctx);
    fq_nmod_clear(gammaff, ffctx->fqctx);
    fq_nmod_clear(t, ffctx->fqctx);

    fq_nmod_mpoly_ctx_clear(ffctx);

    return success;
}


int nmod_mpolyu_gcd_zippel(
    nmod_mpolyu_t G,
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    int success = 0;
    slong i;
    slong ABminshift;
    nmod_mpoly_t content;
    nmod_mpolyu_t Abar, Bbar, Gbar;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    nmod_mpoly_init(content, ctx);
    nmod_mpolyu_init(Abar, A->bits, ctx);
    nmod_mpolyu_init(Bbar, A->bits, ctx);
    nmod_mpolyu_init(Gbar, A->bits, ctx);

    /* compute the content of GCD wrt non main variables */
    nmod_mpoly_set(content, A->coeffs + 0, ctx);
    for (i = 1; i < A->length; i++)
    {
        if (nmod_mpoly_is_one(content, ctx))
            break;
        success = _nmod_mpoly_gcd_zippel(content, content, A->coeffs + i, ctx, 1, randstate);
        if (!success)
            goto finished;
    }
    for (i = 0; i < B->length; i++)
    {
        if (nmod_mpoly_is_one(content, ctx))
            break;
        success = _nmod_mpoly_gcd_zippel(content, content, B->coeffs + i, ctx, 1, randstate);
        if (!success)
            goto finished;
    }

    nmod_mpolyu_divexact_mpoly(Abar, A, content, ctx);
    nmod_mpolyu_divexact_mpoly(Bbar, B, content, ctx);

    ABminshift = FLINT_MIN(Abar->exps[Abar->length - 1], Bbar->exps[Bbar->length - 1]);
    nmod_mpolyu_shift_right(Abar, Abar->exps[Abar->length - 1]);
    nmod_mpolyu_shift_right(Bbar, Bbar->exps[Bbar->length - 1]);

    success = nmod_mpolyu_gcdm_zippel(Gbar, Abar, Bbar, ctx, zinfo, randstate);
    if (!success)
        goto finished;

    nmod_mpolyu_shift_left(Gbar, ABminshift);
    nmod_mpolyu_mul_mpoly(G, Gbar, content, ctx);

    success = 1;

finished:

    nmod_mpolyu_clear(Abar, ctx);
    nmod_mpolyu_clear(Bbar, ctx);
    nmod_mpolyu_clear(Gbar, ctx);
    nmod_mpoly_clear(content, ctx);

    return success;
}



void nmod_mpoly_to_nmod_polydr_keepbits(nmod_polydr_t A, slong * Ashift,
               const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
{
    slong i, shift, off, N;
    slong _Ashift = 0, len = B->length;
    mp_limb_t * coeff = B->coeffs;
    ulong * exp = B->exps;
    mp_bitcnt_t bits = B->bits;

    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, var, bits, ctx->minfo);

    nmod_polydr_zero(A, ctx->ffinfo);
    if (len > 0)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        _Ashift = (exp[N*(len - 1)] >> shift) & mask;
        for (i = 0; i < len; i++)
        {
            ulong k = ((exp[N*i + off] >> shift) & mask) - _Ashift;
            FLINT_ASSERT(((slong)k) >= 0);
            nmod_polydr_set_coeff_ui(A, k, coeff[i], ctx->ffinfo);
        }
    }

    *Ashift = _Ashift;
}

void nmod_mpoly_from_nmod_polydr_keepbits(nmod_mpoly_t A, const nmod_polydr_t B,
                           slong Bshift, slong var, mp_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
{
    slong N;
    slong k;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * one;
    TMP_INIT;

    TMP_START;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(!nmod_polydr_is_zero(B, ctx->ffinfo));
    FLINT_ASSERT(Bshift >= 0);
    FLINT_ASSERT(Bshift + nmod_polydr_degree(B, ctx->ffinfo) >= 0);
    FLINT_ASSERT(1 + FLINT_BIT_COUNT(Bshift + nmod_polydr_degree(B, ctx->ffinfo)) <= bits);

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(one, var, bits, ctx->minfo);

    nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (k = nmod_polydr_degree(B, ctx->ffinfo); k >= 0; k--)
    {
        _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
        mpoly_monomial_mul_ui(Aexp + N*Alen, one, N, k + Bshift);
        Acoeff[Alen] = nmod_polydr_get_coeff_ui(B, k, ctx->ffinfo);
        Alen += Acoeff[Alen] != UWORD(0);
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _nmod_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}

int _nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx,
                                         int keepbits, flint_rand_t randstate)
{
    slong i;
    int ret, success = 0;
    mpoly_zipinfo_t zinfo;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu;
    mp_bitcnt_t new_bits;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT((!keepbits) || A->bits == B->bits);

    FLINT_ASSERT(!nmod_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!nmod_mpoly_is_zero(B, ctx));

    if (ctx->minfo->nvars == 1)
    {
        slong shiftA, shiftB;
        nmod_polydr_t a, b, g;
        nmod_polydr_init(a, ctx->ffinfo);
        nmod_polydr_init(b, ctx->ffinfo);
        nmod_polydr_init(g, ctx->ffinfo);
        nmod_mpoly_to_nmod_polydr_keepbits(a, &shiftA, A, 0, ctx);
        nmod_mpoly_to_nmod_polydr_keepbits(b, &shiftB, B, 0, ctx);
        nmod_polydr_gcd(g, a, b, ctx->ffinfo);
        nmod_mpoly_from_nmod_polydr_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        nmod_polydr_clear(a, ctx->ffinfo);
        nmod_polydr_clear(b, ctx->ffinfo);
        nmod_polydr_clear(g, ctx->ffinfo);
        return 1;
    }

    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    nmod_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    nmod_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        zinfo->perm[i] = i;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    nmod_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_mpolyu_init(Au, new_bits, uctx);
    nmod_mpolyu_init(Bu, new_bits, uctx);
    nmod_mpolyu_init(Gu, new_bits, uctx);

    nmod_mpoly_to_mpolyu_perm(Au, A, zinfo->perm, uctx, ctx);
    nmod_mpoly_to_mpolyu_perm(Bu, B, zinfo->perm, uctx, ctx);

    ret = nmod_mpolyu_gcd_zippel(Gu, Au, Bu, uctx, zinfo, randstate);
    if (ret) {
        nmod_mpoly_from_mpolyu_perm(G, Gu, keepbits, zinfo->perm, uctx, ctx);
        nmod_mpoly_make_monic(G, G, ctx);
        success = 1;
    }

    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Gu, uctx);
    nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    return success;
}


int nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    int success;
    flint_rand_t randstate;

    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
        {
            nmod_mpoly_zero(G, ctx);
        } else
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

    if (ctx->minfo->nvars == 1)
    {
        slong shiftA, shiftB;
        nmod_polydr_t a, b, g;
        nmod_polydr_init(a, ctx->ffinfo);
        nmod_polydr_init(b, ctx->ffinfo);
        nmod_polydr_init(g, ctx->ffinfo);
        nmod_mpoly_to_nmod_polydr_keepbits(a, &shiftA, A, 0, ctx);
        nmod_mpoly_to_nmod_polydr_keepbits(b, &shiftB, B, 0, ctx);
        nmod_polydr_gcd(g, a, b, ctx->ffinfo);
        nmod_mpoly_from_nmod_polydr_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        nmod_polydr_clear(a, ctx->ffinfo);
        nmod_polydr_clear(b, ctx->ffinfo);
        nmod_polydr_clear(g, ctx->ffinfo);
        return 1;
    }

    flint_randinit(randstate);
    success = _nmod_mpoly_gcd_zippel(G, A, B, ctx, 1, randstate);

    flint_randclear(randstate);

    return success;
}
