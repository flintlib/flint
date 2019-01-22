/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"




int fq_nmod_mpolyu_gcdm_zippel_bivar(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    fq_nmod_mpolyun_t An, Bn, H, Ht;
    fq_nmod_poly_t modulus, a, b, c, g;
    int success = 0, changed;
    slong Alastdeg;
    slong Blastdeg;
    slong lastdeg;
    slong var = 0;
    fq_nmod_embed_struct * embed;
    mp_limb_t p = ctx->fqctx->modulus->mod.n;
    slong m = nmod_poly_degree(ctx->fqctx->modulus);
    slong n;
    slong k;

    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(B->bits == A->bits);

    fq_nmod_mpolyun_init(An, A->bits, ctx);
    fq_nmod_mpolyun_init(Bn, A->bits, ctx);

    fq_nmod_mpolyu_cvtto_mpolyun(An, A, var, ctx);
    fq_nmod_mpolyu_cvtto_mpolyun(Bn, B, var, ctx);

    fq_nmod_poly_init(a, ctx->fqctx);
    fq_nmod_poly_init(b, ctx->fqctx);
    fq_nmod_poly_init(c, ctx->fqctx);
    fq_nmod_poly_init(g, ctx->fqctx);
    fq_nmod_mpolyun_content_last(a, An, ctx);
    fq_nmod_mpolyun_content_last(b, Bn, ctx);
    fq_nmod_mpolyun_divexact_last(An, a, ctx);
    fq_nmod_mpolyun_divexact_last(Bn, b, ctx);
    fq_nmod_poly_gcd(c, a, b, ctx->fqctx);
    fq_nmod_poly_gcd(g, fq_nmod_mpolyun_leadcoeff_ref(An, ctx),
                        fq_nmod_mpolyun_leadcoeff_ref(Bn, ctx), ctx->fqctx);

    Alastdeg = fq_nmod_mpolyun_lastdeg(An, ctx);
    Blastdeg = fq_nmod_mpolyun_lastdeg(Bn, ctx);

    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_mpolyun_init(H, A->bits, ctx);
    fq_nmod_mpolyun_init(Ht, A->bits, ctx);

    FLINT_ASSERT(An->bits == B->bits);
    FLINT_ASSERT(An->bits == G->bits);
    FLINT_ASSERT(An->length > 0);
    FLINT_ASSERT(Bn->length > 0);
    FLINT_ASSERT(An->exps[A->length - 1] == 0);
    FLINT_ASSERT(Bn->exps[B->length - 1] == 0);

    embed = (fq_nmod_embed_struct *) flint_malloc(m*sizeof(fq_nmod_embed_struct));

    /* n is the degree of the extension */
    n = (FLINT_BITS/2)/(m*FLINT_BIT_COUNT(p));
    for (n = FLINT_MAX(n, WORD(2)); n < 1000; n++)
    {
        nmod_poly_t ext_modulus;
        fq_nmod_ctx_t ext_fqctx;
        fq_nmod_mpoly_ctx_t ectx;
        fq_nmod_mpolyu_t Aeval, Beval, Geval;
        fq_nmod_t t, geval;

        /* init ectx with modulus of degree m*n */
        nmod_poly_init2(ext_modulus, p, m*n + 1);
        nmod_poly_randtest_sparse_irreducible(ext_modulus, randstate, m*n + 1);
        fq_nmod_ctx_init_modulus(ext_fqctx, ext_modulus, "$");
        fq_nmod_mpoly_ctx_init2(ectx, ctx->minfo->nvars, ORD_LEX, ext_fqctx);
        fq_nmod_ctx_clear(ext_fqctx);
        nmod_poly_clear(ext_modulus);

        _fq_nmod_embed_array_init(embed, ectx->fqctx, ctx->fqctx);

        fq_nmod_mpolyu_init(Aeval, A->bits, ectx);
        fq_nmod_mpolyu_init(Beval, A->bits, ectx);
        fq_nmod_mpolyu_init(Geval, A->bits, ectx);
        fq_nmod_init(geval, ectx->fqctx);
        fq_nmod_init(t, ectx->fqctx);

        for (k = 0; k < m; k++)
        {
            /* make sure reduction does not kill both lc(A) and lc(B) */
            _fq_nmod_embed_sm_to_lg(geval, g, embed + k);
            if (fq_nmod_is_zero(geval, ectx->fqctx))
                goto outer_continue;

            /* make sure reduction does not kill either A or B */
            fq_nmod_mpolyun_redto_fq_nmod_mpolyu(Aeval, An, ectx, ctx, embed + k);
            fq_nmod_mpolyun_redto_fq_nmod_mpolyu(Beval, Bn, ectx, ctx, embed + k);
            if (Aeval->length == 0 || Beval->length == 0)
                goto outer_continue;

            FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Aeval, ectx));
            FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(Beval, ectx));

            fq_nmod_mpolyu_gcdp_zippel_univar(Geval, Aeval, Beval, ectx);

            if (fq_nmod_mpolyu_is_one(Geval, ectx))
            {
                fq_nmod_mpolyu_cvtfrom_poly_notmain(G, c, var, ctx);

                fq_nmod_clear(geval, ectx->fqctx);
                fq_nmod_clear(t, ectx->fqctx);
                fq_nmod_mpolyu_clear(Aeval, ectx);
                fq_nmod_mpolyu_clear(Beval, ectx);
                fq_nmod_mpolyu_clear(Geval, ectx);
                _fq_nmod_embed_array_clear(embed, m);
                fq_nmod_mpoly_ctx_clear(ectx);

                success = 1;
                goto finished;
            }

            FLINT_ASSERT(Geval->length > 0);

            if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
            {
                if (Geval->exps[0] > H->exps[0])
                {
                    goto outer_continue;
                }
                else if (Geval->exps[0] < H->exps[0])
                {
                    fq_nmod_poly_one(modulus, ctx->fqctx);
                }
            }

            fq_nmod_inv(t, fq_nmod_mpolyu_leadcoeff_ref(Geval, ectx), ectx->fqctx);
            fq_nmod_mul(t, t, geval, ectx->fqctx);
            fq_nmod_mpolyu_scalar_mul_fq_nmod(Geval, t, ectx);

            if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
            {
                changed = fq_nmod_mpolyun_addinterp_lgprime(&lastdeg,
                                  H, Ht, modulus, ctx, Geval, ectx, embed + k);
                if (!changed)
                {
                    fq_nmod_mpolyun_content_last(a, H, ctx);
                    fq_nmod_mpolyun_mul_poly(Ht, H, c, ctx);
                    fq_nmod_mpolyun_divexact_last(Ht, a, ctx);

                    fq_nmod_mpolyu_cvtfrom_mpolyun(G, Ht, var, ctx);

                    if (!fq_nmod_mpolyu_divides(A, G, ctx))
                        goto outer_continue;
                    if (!fq_nmod_mpolyu_divides(B, G, ctx))
                        goto outer_continue;

                    fq_nmod_clear(geval, ectx->fqctx);
                    fq_nmod_clear(t, ectx->fqctx);
                    fq_nmod_mpolyu_clear(Aeval, ectx);
                    fq_nmod_mpolyu_clear(Beval, ectx);
                    fq_nmod_mpolyu_clear(Geval, ectx);
                    _fq_nmod_embed_array_clear(embed, m);
                    fq_nmod_mpoly_ctx_clear(ectx);

                    success = 1;
                    goto finished;
                }
            }
            else
            {
                fq_nmod_mpolyun_startinterp_lgprime(H, ctx, Geval, ectx, embed + k);
                lastdeg = fq_nmod_mpolyun_lastdeg(H, ctx);
            }
            fq_nmod_poly_mul(modulus, modulus, (embed + k)->h, ctx->fqctx);

            if (lastdeg > Alastdeg || lastdeg > Blastdeg)
            {
                fq_nmod_poly_one(modulus, ctx->fqctx);
                goto outer_continue;
            }

outer_continue:
            (void)(NULL);
        }

        fq_nmod_clear(geval, ectx->fqctx);
        fq_nmod_clear(t, ectx->fqctx);
        fq_nmod_mpolyu_clear(Aeval, ectx);
        fq_nmod_mpolyu_clear(Beval, ectx);
        fq_nmod_mpolyu_clear(Geval, ectx);
        _fq_nmod_embed_array_clear(embed, m);
        fq_nmod_mpoly_ctx_clear(ectx);
    }

    success = 0;

finished:

    fq_nmod_poly_clear(a, ctx->fqctx);
    fq_nmod_poly_clear(b, ctx->fqctx);
    fq_nmod_poly_clear(c, ctx->fqctx);
    fq_nmod_poly_clear(g, ctx->fqctx);
    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_mpolyun_clear(An, ctx);
    fq_nmod_mpolyun_clear(Bn, ctx);
    fq_nmod_mpolyun_clear(H, ctx);
    fq_nmod_mpolyun_clear(Ht, ctx);

    flint_free(embed);

    return success;
}



int fq_nmod_mpolyu_gcdm_zippel(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
/*
    slong degbound;
    slong coeff_deg_bound;
    slong lastdeg;
*/
    int success/*, changed*/;
/*
    slong deg = 4;
    nmod_mpolyun_t An, Bn, Hn, Ht;
    fq_nmod_mpoly_ctx_t ffctx;
    fq_nmod_mpolyu_t Aff, Bff, Gff, Gform;
    nmod_poly_t modulus, gamma, hc;
    fq_nmod_t t, gammaff;
*/
    FLINT_ASSERT(G->bits == A->bits);
    FLINT_ASSERT(B->bits == A->bits);


    success = fq_nmod_mpolyu_gcdp_zippel(G, A, B, ctx->minfo->nvars - 1,
                                                        ctx, zinfo, randstate);
    if (success)
        return 1;

    /* bivariate more comfortable separated */
    if (ctx->minfo->nvars == 1)
        return fq_nmod_mpolyu_gcdm_zippel_bivar(G, A, B, ctx, zinfo, randstate);



    return 0;
}


int fq_nmod_mpolyu_gcd_zippel(
    fq_nmod_mpolyu_t G,
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    int success = 0;
    slong i;
    slong ABminshift;
    fq_nmod_mpoly_t content;
    fq_nmod_mpolyu_t Abar, Bbar, Gbar;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    fq_nmod_mpoly_init(content, ctx);
    fq_nmod_mpolyu_init(Abar, A->bits, ctx);
    fq_nmod_mpolyu_init(Bbar, A->bits, ctx);
    fq_nmod_mpolyu_init(Gbar, A->bits, ctx);

    /* compute the content of GCD wrt non main variables */
    fq_nmod_mpoly_set(content, A->coeffs + 0, ctx);
    for (i = 1; i < A->length; i++)
    {
        if (fq_nmod_mpoly_is_one(content, ctx))
            break;
        success = _fq_nmod_mpoly_gcd_zippel(content, content, A->coeffs + i, ctx, 1, randstate);
        if (!success)
            goto finished;
    }
    for (i = 0; i < B->length; i++)
    {
        if (fq_nmod_mpoly_is_one(content, ctx))
            break;
        success = _fq_nmod_mpoly_gcd_zippel(content, content, B->coeffs + i, ctx, 1, randstate);
        if (!success)
            goto finished;
    }

    fq_nmod_mpolyu_divexact_mpoly(Abar, A, content, ctx);
    fq_nmod_mpolyu_divexact_mpoly(Bbar, B, content, ctx);

    ABminshift = FLINT_MIN(Abar->exps[Abar->length - 1], Bbar->exps[Bbar->length - 1]);
    fq_nmod_mpolyu_shift_right(Abar, Abar->exps[Abar->length - 1]);
    fq_nmod_mpolyu_shift_right(Bbar, Bbar->exps[Bbar->length - 1]);

    success = fq_nmod_mpolyu_gcdm_zippel(Gbar, Abar, Bbar, ctx, zinfo, randstate);
    if (!success)
        goto finished;

    fq_nmod_mpolyu_shift_left(Gbar, ABminshift);
    fq_nmod_mpolyu_mul_mpoly(G, Gbar, content, ctx);

    success = 1;

finished:

    fq_nmod_mpolyu_clear(Abar, ctx);
    fq_nmod_mpolyu_clear(Bbar, ctx);
    fq_nmod_mpolyu_clear(Gbar, ctx);
    fq_nmod_mpoly_clear(content, ctx);

    return success;
}


void fq_nmod_mpoly_to_fq_nmod_poly_keepbits(fq_nmod_poly_t A, slong * Ashift,
             const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, shift, off, N;
    slong _Ashift = 0, len = B->length;
    fq_nmod_struct * Bcoeff = B->coeffs;
    ulong * exp = B->exps;
    mp_bitcnt_t bits = B->bits;

    FLINT_ASSERT(bits <= FLINT_BITS);

    N = mpoly_words_per_exp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, var, bits, ctx->minfo);

    fq_nmod_poly_zero(A, ctx->fqctx);
    if (len > 0)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        _Ashift = (exp[N*(len - 1)] >> shift) & mask;
        for (i = 0; i < len; i++)
        {
            ulong k = ((exp[N*i + off] >> shift) & mask) - _Ashift;
            FLINT_ASSERT(((slong)k) >= 0);
            fq_nmod_poly_set_coeff(A, k, Bcoeff + i, ctx->fqctx);
        }
    }

    *Ashift = _Ashift;
}

void fq_nmod_mpoly_from_fq_nmod_poly_keepbits(fq_nmod_mpoly_t A,
         const fq_nmod_poly_t B, slong Bshift, slong var, mp_bitcnt_t bits,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;
    slong k;
    slong Alen;
    fq_nmod_struct * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * one;
    TMP_INIT;

    TMP_START;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(!fq_nmod_poly_is_zero(B, ctx->fqctx));
    FLINT_ASSERT(Bshift >= 0);
    FLINT_ASSERT(Bshift + fq_nmod_poly_degree(B, ctx->fqctx) >= 0);
    FLINT_ASSERT(1 + FLINT_BIT_COUNT(Bshift + fq_nmod_poly_degree(B, ctx->fqctx)) <= bits);

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(one, var, bits, ctx->minfo);

    fq_nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (k = fq_nmod_poly_degree(B, ctx->fqctx); k >= 0; k--)
    {
        _fq_nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N, ctx->fqctx);
        mpoly_monomial_mul_ui(Aexp + N*Alen, one, N, k + Bshift);
        fq_nmod_set(Acoeff + Alen, B->coeffs + k, ctx->fqctx);
        Alen += !fq_nmod_is_zero(Acoeff + Alen, ctx->fqctx);
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fq_nmod_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}




int _fq_nmod_mpoly_gcd_zippel(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                      const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx,
                                         int keepbits, flint_rand_t randstate)
{
    slong i;
    int ret, success = 0;
    mpoly_zipinfo_t zinfo;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyu_t Au, Bu, Gu;
    mp_bitcnt_t new_bits;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT((!keepbits) || A->bits == B->bits);

    FLINT_ASSERT(!fq_nmod_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fq_nmod_mpoly_is_zero(B, ctx));

    if (ctx->minfo->nvars == 1)
    {
        slong shiftA, shiftB;
        fq_nmod_poly_t a, b, g;
        fq_nmod_poly_init(a, ctx->fqctx);
        fq_nmod_poly_init(b, ctx->fqctx);
        fq_nmod_poly_init(g, ctx->fqctx);
        fq_nmod_mpoly_to_fq_nmod_poly_keepbits(a, &shiftA, A, 0, ctx);
        fq_nmod_mpoly_to_fq_nmod_poly_keepbits(b, &shiftB, B, 0, ctx);
        fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
        fq_nmod_mpoly_from_fq_nmod_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        fq_nmod_poly_clear(a, ctx->fqctx);
        fq_nmod_poly_clear(b, ctx->fqctx);
        fq_nmod_poly_clear(g, ctx->fqctx);
        return 1;
    }

    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    fq_nmod_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    fq_nmod_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        zinfo->perm[i] = i;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    fq_nmod_mpoly_ctx_init2(uctx, ctx->minfo->nvars - 1, ORD_LEX, ctx->fqctx);
    fq_nmod_mpolyu_init(Au, new_bits, uctx);
    fq_nmod_mpolyu_init(Bu, new_bits, uctx);
    fq_nmod_mpolyu_init(Gu, new_bits, uctx);

    fq_nmod_mpoly_to_mpolyu_perm(Au, A, zinfo->perm, uctx, ctx);
    fq_nmod_mpoly_to_mpolyu_perm(Bu, B, zinfo->perm, uctx, ctx);

    ret = fq_nmod_mpolyu_gcd_zippel(Gu, Au, Bu, uctx, zinfo, randstate);
    if (ret)
    {
        fq_nmod_mpoly_from_mpolyu_perm(G, Gu, keepbits, zinfo->perm, uctx, ctx);
        fq_nmod_mpoly_make_monic(G, G, ctx);
        success = 1;
    }

    fq_nmod_mpolyu_clear(Au, uctx);
    fq_nmod_mpolyu_clear(Bu, uctx);
    fq_nmod_mpolyu_clear(Gu, uctx);
    fq_nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    return success;
}

int fq_nmod_mpoly_gcd_zippel(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                        const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    flint_rand_t randstate;

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

    if (ctx->minfo->nvars == 1)
    {
        slong shiftA, shiftB;
        fq_nmod_poly_t a, b, g;
        fq_nmod_poly_init(a, ctx->fqctx);
        fq_nmod_poly_init(b, ctx->fqctx);
        fq_nmod_poly_init(g, ctx->fqctx);
        fq_nmod_mpoly_to_fq_nmod_poly_keepbits(a, &shiftA, A, 0, ctx);
        fq_nmod_mpoly_to_fq_nmod_poly_keepbits(b, &shiftB, B, 0, ctx);
        fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
        fq_nmod_mpoly_from_fq_nmod_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        fq_nmod_poly_clear(a, ctx->fqctx);
        fq_nmod_poly_clear(b, ctx->fqctx);
        fq_nmod_poly_clear(g, ctx->fqctx);
        return 1;
    }

    flint_randinit(randstate);
    success = _fq_nmod_mpoly_gcd_zippel(G, A, B, ctx, 1, randstate);
    flint_randclear(randstate);
    return success;
}
