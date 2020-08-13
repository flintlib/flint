/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/* from gcd.c */
void fq_nmod_mpoly_evals(
    fq_nmod_poly_struct * out,
    const int * ignore,
    const fq_nmod_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx);

void mpoly_gcd_info_set_estimates_fq_nmod_mpoly(
    mpoly_gcd_info_t I,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx);


/*********************** Easy when B is a monomial ***************************/
static int _try_monomial_gcd(
    fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
    fq_nmod_mpoly_t Abar, flint_bitcnt_t Abarbits,
    fq_nmod_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * minAfields, * minAdegs, * minBdegs;
    fq_nmod_mpoly_t _G, _Abar, _Bbar;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length == 1);

    fq_nmod_mpoly_init(_G, ctx);
    fq_nmod_mpoly_init(_Abar, ctx);
    fq_nmod_mpoly_init(_Bbar, ctx);

    TMP_START;

    /* get the field-wise minimum of A */
    minAfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(minAfields + i);
    mpoly_min_fields_fmpz(minAfields, A->exps, A->length, A->bits, ctx->minfo);

    /* unpack to get the min degrees of each variable in A */
    minAdegs = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(minAdegs + i);
    mpoly_get_monomial_ffmpz_unpacked_ffmpz(minAdegs, minAfields, ctx->minfo);

    /* get the degree of each variable in B */
    minBdegs = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(minBdegs + i);
    mpoly_get_monomial_ffmpz(minBdegs, B->exps, B->bits, ctx->minfo);

    /* compute the degree of each variable in G */
    _fmpz_vec_min_inplace(minBdegs, minAdegs, ctx->minfo->nvars);

    fq_nmod_mpoly_fit_length(_G, 1, ctx);
    fq_nmod_mpoly_fit_bits(_G, Gbits, ctx);
    _G->bits = Gbits;
    mpoly_set_monomial_ffmpz(_G->exps, minBdegs, Gbits, ctx->minfo);
    fq_nmod_one(_G->coeffs + 0, ctx->fqctx);
    _G->length = 1;

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(minAfields + i);
    }
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_clear(minAdegs + i);
        fmpz_clear(minBdegs + i);
    }

    TMP_END;

    fq_nmod_mpoly_divides(_Abar, A, _G, ctx);
    fq_nmod_mpoly_divides(_Bbar, B, _G, ctx);

    fq_nmod_mpoly_swap(G, _G, ctx);
    fq_nmod_mpoly_swap(Abar, _Abar, ctx);
    fq_nmod_mpoly_swap(Bbar, _Bbar, ctx);

    fq_nmod_mpoly_clear(_G, ctx);
    fq_nmod_mpoly_clear(_Abar, ctx);
    fq_nmod_mpoly_clear(_Bbar, ctx);

    return 1;
}


/********************** See if cofactors are monomials ***********************/
static int _try_monomial_cofactors(
    fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
    fq_nmod_mpoly_t Abar, flint_bitcnt_t Abarbits,
    fq_nmod_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong NA, NG;
    slong nvars = ctx->minfo->nvars;
    fmpz * Abarexps, * Bbarexps, * Texps;
    fq_nmod_t t1, t2, a0, b0;
    fq_nmod_mpoly_t T;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (A->length != B->length)
        return 0;

    fq_nmod_init(t1, ctx->fqctx);
    fq_nmod_init(t2, ctx->fqctx);

    for (i = A->length - 1; i > 0; i--)
    {
        fq_nmod_mul(t1, A->coeffs + 0, B->coeffs + i, ctx->fqctx);
        fq_nmod_mul(t2, B->coeffs + 0, A->coeffs + i, ctx->fqctx);
        success = fq_nmod_equal(t1, t2, ctx->fqctx);
        if (!success)
            goto cleanup;
    }

    TMP_START;

    fq_nmod_init(a0, ctx->fqctx);
    fq_nmod_init(b0, ctx->fqctx);
    fq_nmod_set(a0, A->coeffs + 0, ctx->fqctx);
    fq_nmod_set(b0, B->coeffs + 0, ctx->fqctx);

    Abarexps = (fmpz *) TMP_ALLOC(3*nvars*sizeof(fmpz));
    Bbarexps = Abarexps + 1*nvars;
    Texps    = Abarexps + 2*nvars;
    for (j = 0; j < nvars; j++)
    {
        fmpz_init(Abarexps + j);
        fmpz_init(Bbarexps + j);
        fmpz_init(Texps + j);
    }

    success = mpoly_monomial_cofactors(Abarexps, Bbarexps, A->exps, A->bits,
                                      B->exps, B->bits, A->length, ctx->minfo);
    if (!success)
        goto cleanup_tmp;

    fq_nmod_mpoly_init3(T, A->length, Gbits, ctx);
    NG = mpoly_words_per_exp(Gbits, ctx->minfo);
    NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    fq_nmod_inv(t1, A->coeffs + 0, ctx->fqctx);
    T->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ffmpz(Texps, A->exps + NA*i, A->bits, ctx->minfo);
        _fmpz_vec_sub(Texps, Texps, Abarexps, nvars);
        mpoly_set_monomial_ffmpz(T->exps + NG*i, Texps, Gbits, ctx->minfo);
        fq_nmod_mul(T->coeffs + i, A->coeffs + i, t1, ctx->fqctx);
    }
    fq_nmod_mpoly_swap(G, T, ctx);
    fq_nmod_mpoly_clear(T, ctx);

    fq_nmod_mpoly_fit_length(Abar, 1, ctx);
    fq_nmod_mpoly_fit_bits(Abar, Abarbits, ctx);
    Abar->bits = Abarbits;
    mpoly_set_monomial_ffmpz(Abar->exps, Abarexps, Abarbits, ctx->minfo);
    fq_nmod_swap(Abar->coeffs + 0, a0, ctx->fqctx);
    _fq_nmod_mpoly_set_length(Abar, 1, ctx);

    fq_nmod_mpoly_fit_length(Bbar, 1, ctx);
    fq_nmod_mpoly_fit_bits(Bbar, Bbarbits, ctx);
    Bbar->bits = Bbarbits;
    mpoly_set_monomial_ffmpz(Bbar->exps, Bbarexps, Bbarbits, ctx->minfo);
    fq_nmod_swap(Bbar->coeffs + 0, b0, ctx->fqctx);
    _fq_nmod_mpoly_set_length(Bbar, 1, ctx);

    success = 1;

cleanup_tmp:

    fq_nmod_clear(a0, ctx->fqctx);
    fq_nmod_clear(b0, ctx->fqctx);

    for (j = 0; j < nvars; j++)
    {
        fmpz_clear(Abarexps + j);
        fmpz_clear(Bbarexps + j);
        fmpz_clear(Texps + j);
    }

    TMP_END;

cleanup:

    fq_nmod_clear(t1, ctx->fqctx);
    fq_nmod_clear(t2, ctx->fqctx);

    return success;
}


/********* Assume B has length one when converted to univar format ***********/
static int _try_missing_var(
    fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
    fq_nmod_mpoly_t Abar, flint_bitcnt_t Abarbits,
    fq_nmod_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    slong var,
    const fq_nmod_mpoly_t A, ulong Ashift,
    const fq_nmod_mpoly_t B, ulong Bshift,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fq_nmod_mpoly_t tG, tAbar, tBbar;
    fq_nmod_mpoly_univar_t Ax;

    fq_nmod_mpoly_init(tG, ctx);
    fq_nmod_mpoly_init(tAbar, ctx);
    fq_nmod_mpoly_init(tBbar, ctx);
    fq_nmod_mpoly_univar_init(Ax, ctx);

    fq_nmod_mpoly_to_univar(Ax, A, var, ctx);

    FLINT_ASSERT(Ax->length > 0);
    success = _fq_nmod_mpoly_gcd(tG, Gbits, B, Ax->coeffs + 0, ctx);
    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _fq_nmod_mpoly_gcd(tG, Gbits, tG, Ax->coeffs + i, ctx);
        if (!success)
            goto cleanup;
    }

    _mpoly_gen_shift_left(tG->exps, tG->bits, tG->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

    success = fq_nmod_mpoly_divides(tAbar, A, tG, ctx);
    FLINT_ASSERT(success);
    success = fq_nmod_mpoly_divides(tBbar, B, tG, ctx);
    FLINT_ASSERT(success);

    fq_nmod_mpoly_swap(G, tG, ctx);
    fq_nmod_mpoly_swap(Abar, tAbar, ctx);
    fq_nmod_mpoly_swap(Bbar, tBbar, ctx);

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(tG, ctx);
    fq_nmod_mpoly_clear(tAbar, ctx);
    fq_nmod_mpoly_clear(tBbar, ctx);
    fq_nmod_mpoly_univar_clear(Ax, ctx);

    return success;
}


/******************* Test if B divides A or A divides B **********************/
/*
    Test if B divides A or A divides B
        TODO: incorporate deflation
*/
static int _try_divides(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A, int try_a,
    const fq_nmod_mpoly_t B, int try_b,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    fq_nmod_mpoly_t Q;

    fq_nmod_mpoly_init(Q, ctx);

    if (try_b && fq_nmod_mpoly_divides(Q, A, B, ctx))
    {
        fq_nmod_mpoly_set(G, B, ctx);
        fq_nmod_mpoly_swap(Abar, Q, ctx);
        fq_nmod_mpoly_one(Bbar, ctx);
        success = 1;
        goto cleanup;
    }

    if (try_a && fq_nmod_mpoly_divides(Q, B, A, ctx))
    {
        fq_nmod_mpoly_set(G, A, ctx);
        fq_nmod_mpoly_swap(Bbar, Q, ctx);
        fq_nmod_mpoly_one(Abar, ctx);
        success = 1;
        goto cleanup;
    }

    success = 0;

cleanup:

    fq_nmod_mpoly_clear(Q, ctx);

    return success;
}


/********************** Hit A and B with zippel ******************************/
static int _try_zippel(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    mpoly_zipinfo_t zinfo;
    flint_bitcnt_t wbits;
    flint_rand_t randstate;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    fq_nmod_mpoly_t Ac, Bc, Gc, Abarc, Bbarc;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    if (!I->can_use_zippel)
        return 0;

    FLINT_ASSERT(m >= WORD(2));
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    flint_randinit(randstate);

    /* interpolation will continue in m variables */
    mpoly_zipinfo_init(zinfo, m);

    /* uctx is context for Fq[y_1,...,y_{m-1}]*/
    fq_nmod_mpoly_ctx_init(uctx, m - 1, ORD_LEX, ctx->fqctx);

    /* fill in a valid zinfo->perm and degrees */
    for (i = 0; i < m; i++)
    {
        k = I->zippel_perm[i];
        zinfo->perm[i] = k;
        zinfo->Adegs[i] = I->Adeflate_deg[k];
        zinfo->Bdegs[i] = I->Bdeflate_deg[k];
        FLINT_ASSERT(I->Adeflate_deg[k] != 0);
        FLINT_ASSERT(I->Bdeflate_deg[k] != 0);
    }

    wbits = FLINT_MAX(A->bits, B->bits);

    fq_nmod_mpolyu_init(Au, wbits, uctx);
    fq_nmod_mpolyu_init(Bu, wbits, uctx);
    fq_nmod_mpolyu_init(Gu, wbits, uctx);
    fq_nmod_mpolyu_init(Abaru, wbits, uctx);
    fq_nmod_mpolyu_init(Bbaru, wbits, uctx);
    fq_nmod_mpoly_init3(Ac, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Bc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Gc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Abarc, 0, wbits, uctx);
    fq_nmod_mpoly_init3(Bbarc, 0, wbits, uctx);

    fq_nmod_mpoly_to_mpolyu_perm_deflate(Au, uctx, A, ctx, zinfo->perm,
                                                      I->Amin_exp, I->Gstride);
    fq_nmod_mpoly_to_mpolyu_perm_deflate(Bu, uctx, B, ctx, zinfo->perm,
                                                      I->Bmin_exp, I->Gstride);

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

    success = _fq_nmod_mpoly_gcd_cofactors(Gc, wbits, Abarc, wbits,
                                                   Bbarc, wbits, Ac, Bc, uctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);
    fq_nmod_mpolyu_mul_mpoly_inplace(Abaru, Abarc, uctx);
    fq_nmod_mpolyu_mul_mpoly_inplace(Bbaru, Bbarc, uctx);

    fq_nmod_mpoly_from_mpolyu_perm_inflate(G, I->Gbits, ctx, Gu, uctx,
                                         zinfo->perm, I->Gmin_exp, I->Gstride);
    fq_nmod_mpoly_from_mpolyu_perm_inflate(Abar, I->Abarbits, ctx, Abaru, uctx,
                                      zinfo->perm, I->Abarmin_exp, I->Gstride);
    fq_nmod_mpoly_from_mpolyu_perm_inflate(Bbar, I->Bbarbits, ctx, Bbaru, uctx,
                                      zinfo->perm, I->Bbarmin_exp, I->Gstride);
    success = 1;

cleanup:

    fq_nmod_mpolyu_clear(Au, uctx);
    fq_nmod_mpolyu_clear(Bu, uctx);
    fq_nmod_mpolyu_clear(Gu, uctx);
    fq_nmod_mpolyu_clear(Abaru, uctx);
    fq_nmod_mpolyu_clear(Bbaru, uctx);
    fq_nmod_mpoly_clear(Ac, uctx);
    fq_nmod_mpoly_clear(Bc, uctx);
    fq_nmod_mpoly_clear(Gc, uctx);
    fq_nmod_mpoly_clear(Abarc, uctx);
    fq_nmod_mpoly_clear(Bbarc, uctx);

    fq_nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    return success;
}


/*********************** Hit A and B with brown ******************************/
static int _try_brown(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    mpoly_gcd_info_t I,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong m = I->mvars;
    flint_bitcnt_t wbits;
    fq_nmod_mpoly_ctx_t nctx;
    fq_nmod_mpolyn_t An, Bn, Gn, Abarn, Bbarn;

    if (!I->can_use_brown)
        return 0;

    FLINT_ASSERT(m >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    wbits = FLINT_MAX(A->bits, B->bits);

    fq_nmod_mpoly_ctx_init(nctx, m, ORD_LEX, ctx->fqctx);
    fq_nmod_mpolyn_init(An, wbits, nctx);
    fq_nmod_mpolyn_init(Bn, wbits, nctx);
    fq_nmod_mpolyn_init(Gn, wbits, nctx);
    fq_nmod_mpolyn_init(Abarn, wbits, nctx);
    fq_nmod_mpolyn_init(Bbarn, wbits, nctx);

    fq_nmod_mpoly_to_mpolyn_perm_deflate(An, nctx, A, ctx,
                                       I->brown_perm, I->Amin_exp, I->Gstride);
    fq_nmod_mpoly_to_mpolyn_perm_deflate(Bn, nctx, B, ctx,
                                       I->brown_perm, I->Bmin_exp, I->Gstride);

    FLINT_ASSERT(An->bits == wbits);
    FLINT_ASSERT(Bn->bits == wbits);
    FLINT_ASSERT(An->length > 1);
    FLINT_ASSERT(Bn->length > 1);

    success = fq_nmod_mpolyn_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn,
                                                                  m - 1, nctx);
    if (!success)
    {
        fq_nmod_mpoly_to_mpolyn_perm_deflate(An, nctx, A, ctx,
                                       I->brown_perm, I->Amin_exp, I->Gstride);
        fq_nmod_mpoly_to_mpolyn_perm_deflate(Bn, nctx, B, ctx,
                                       I->brown_perm, I->Bmin_exp, I->Gstride);
        success = fq_nmod_mpolyn_gcd_brown_lgprime(Gn, Abarn, Bbarn, An, Bn,
                                                                  m - 1, nctx);
    }

    if (!success)
        goto cleanup;

    fq_nmod_mpoly_from_mpolyn_perm_inflate(G, I->Gbits, ctx, Gn, nctx,
                                       I->brown_perm, I->Gmin_exp, I->Gstride);
    fq_nmod_mpoly_from_mpolyn_perm_inflate(Abar, I->Abarbits, ctx, Abarn, nctx,
                                    I->brown_perm, I->Abarmin_exp, I->Gstride);
    fq_nmod_mpoly_from_mpolyn_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarn, nctx,
                                    I->brown_perm, I->Bbarmin_exp, I->Gstride);
    success = 1;

cleanup:

    fq_nmod_mpolyn_clear(An, nctx);
    fq_nmod_mpolyn_clear(Bn, nctx);
    fq_nmod_mpolyn_clear(Gn, nctx);
    fq_nmod_mpolyn_clear(Abarn, nctx);
    fq_nmod_mpolyn_clear(Bbarn, nctx);
    fq_nmod_mpoly_ctx_clear(nctx);

    return success;
}

static void _fq_nmod_vec_scalar_div_fq_nmod(
    fq_nmod_struct * A,
    const fq_nmod_struct * B,
    slong len,
    const fq_nmod_t c,
    const fq_nmod_ctx_t fqctx)
{
    slong i;
    fq_nmod_t cinv;

    fq_nmod_init(cinv, fqctx);
    fq_nmod_inv(cinv, c, fqctx);

    for (i = 0; i < len; i++)
        fq_nmod_mul(A + i, B + i, cinv, fqctx);

    fq_nmod_clear(cinv, fqctx);
}



/*
    The function must pack its answer into bits = Gbits <= FLINT_BITS
    Both A and B have to be packed into bits <= FLINT_BITS

    return is 1 for success, 0 for failure.
*/
int _fq_nmod_mpoly_gcd_cofactors(
    fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
    fq_nmod_mpoly_t Abar, flint_bitcnt_t Abarbits,
    fq_nmod_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong v_in_both;
    slong v_in_either;
    slong v_in_A_only;
    slong v_in_B_only;
    slong j;
    slong nvars = ctx->minfo->nvars;
    mpoly_gcd_info_t I;
#if WANT_ASSERT
    fq_nmod_mpoly_t T, Asave, Bsave;

    fq_nmod_mpoly_init(T, ctx);
    fq_nmod_mpoly_init(Asave, ctx);
    fq_nmod_mpoly_init(Bsave, ctx);
    fq_nmod_mpoly_set(Asave, A, ctx);
    fq_nmod_mpoly_set(Bsave, B, ctx);
#endif

    mpoly_gcd_info_init(I, nvars);
    I->Gbits = Gbits;
    I->Abarbits = Abarbits;
    I->Bbarbits = Bbarbits;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Gbits <= FLINT_BITS);
    FLINT_ASSERT(Abarbits <= FLINT_BITS);
    FLINT_ASSERT(Bbarbits <= FLINT_BITS);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    if (A->length == 1)
    {
        _try_monomial_gcd(G, I->Gbits, Bbar, I->Bbarbits, Abar, I->Abarbits,
                                                                    B, A, ctx);
        goto successful;
    }
    else if (B->length == 1)
    {
        _try_monomial_gcd(G, I->Gbits, Abar, I->Abarbits, Bbar, I->Bbarbits,
                                                                    A, B, ctx);
        goto successful;
    }

    /* entries of I are all now invalid */

    mpoly_gcd_info_limits(I->Amax_exp, I->Amin_exp, I->Alead_count,
                      I->Atail_count, A->exps, A->bits, A->length, ctx->minfo);
    mpoly_gcd_info_limits(I->Bmax_exp, I->Bmin_exp, I->Blead_count,
                      I->Btail_count, B->exps, B->bits, B->length, ctx->minfo);

    /* set ess(p) := p/term_content(p) */

    /* check if the cofactors could be monomials, i.e. ess(A) == ess(B) */
    for (j = 0; j < nvars; j++)
    {
        if (I->Amax_exp[j] - I->Amin_exp[j] != I->Bmax_exp[j] - I->Bmin_exp[j])
            goto skip_monomial_cofactors;
    }
    if (_try_monomial_cofactors(G, I->Gbits, Abar, I->Abarbits,
                                             Bbar, I->Bbarbits, A, B, ctx))
    {
        goto successful;
    }

skip_monomial_cofactors:

    mpoly_gcd_info_stride(I->Gstride,
            A->exps, A->bits, A->length, I->Amax_exp, I->Amin_exp,
            B->exps, B->bits, B->length, I->Bmax_exp, I->Bmin_exp, ctx->minfo);

    for (j = 0; j < nvars; j++)
    {
        ulong t = I->Gstride[j];

        if (t == 0)
        {
            FLINT_ASSERT(  I->Amax_exp[j] == I->Amin_exp[j]
                        || I->Bmax_exp[j] == I->Bmin_exp[j]);
        }
        else
        {
            FLINT_ASSERT((I->Amax_exp[j] - I->Amin_exp[j]) % t == 0);
            FLINT_ASSERT((I->Bmax_exp[j] - I->Bmin_exp[j]) % t == 0);
        }

        I->Adeflate_deg[j] = t == 0 ? 0 : (I->Amax_exp[j] - I->Amin_exp[j])/t;
        I->Bdeflate_deg[j] = t == 0 ? 0 : (I->Bmax_exp[j] - I->Bmin_exp[j])/t;

        t = FLINT_MIN(I->Amin_exp[j], I->Bmin_exp[j]);
        I->Gmin_exp[j] = t;
        I->Abarmin_exp[j] = I->Amin_exp[j] - t;
        I->Bbarmin_exp[j] = I->Bmin_exp[j] - t;
    }

    /*
        The following are now valid:
            I->Amax_exp, I->Amin_exp, I->Alead_count, I->Atail_count,
            I->Bmax_exp, I->Bmin_exp, I->Blead_count, I->Btail_count,
            I->Gstride
            I->Adeflate_deg
            I->Bdeflate_deg
            I->Gmin_exp
    */

    /* check if ess(A) and ess(B) have a variable v_in_both in common */
    v_in_both = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (I->Amax_exp[j] > I->Amin_exp[j] && I->Bmax_exp[j] > I->Bmin_exp[j])
        {
            v_in_both = j;
            break;
        }
    }
    if (v_in_both == -WORD(1))
    {
        /*
            The variables in ess(A) and ess(B) are disjoint.
            gcd is trivial to compute.
        */

calculate_trivial_gcd:

        if (Abar == B && Bbar == A)
        {
            fq_nmod_mpoly_set(Abar, B, ctx);
            fq_nmod_mpoly_set(Bbar, A, ctx);
            fq_nmod_mpoly_swap(Abar, Bbar, ctx);
        }
        else if (Abar == B && Bbar != A)
        {
            fq_nmod_mpoly_set(Bbar, B, ctx);
            fq_nmod_mpoly_set(Abar, A, ctx);
        }
        else
        {
            fq_nmod_mpoly_set(Abar, A, ctx);
            fq_nmod_mpoly_set(Bbar, B, ctx);
        }

        fq_nmod_mpoly_fit_length(G, 1, ctx);
        fq_nmod_mpoly_fit_bits(G, Gbits, ctx);
        G->bits = I->Gbits;
        mpoly_set_monomial_ui(G->exps, I->Gmin_exp, I->Gbits, ctx->minfo);
        fq_nmod_one(G->coeffs + 0, ctx->fqctx);
        _fq_nmod_mpoly_set_length(G, 1, ctx);

        mpoly_monomials_shift_right_ui(Abar->exps, Abar->bits, Abar->length,
                                                      I->Gmin_exp, ctx->minfo);
        mpoly_monomials_shift_right_ui(Bbar->exps, Bbar->bits, Bbar->length,
                                                      I->Gmin_exp, ctx->minfo);
        goto successful;
    }

    /* check if ess(A) and ess(B) depend on another variable v_in_either */
    FLINT_ASSERT(0 <= v_in_both);
    FLINT_ASSERT(v_in_both < nvars);

    v_in_either = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (j == v_in_both)
            continue;

        if (I->Amax_exp[j] > I->Amin_exp[j] || I->Bmax_exp[j] > I->Bmin_exp[j])
        {
            v_in_either = j;
            break;
        }
    }

    if (v_in_either == -WORD(1))
    {
        /*
            The ess(A) and ess(B) depend on only one variable v_in_both
            Calculate gcd using univariates
        */
        fq_nmod_poly_t a, b, g, t, r;

        fq_nmod_poly_init(a, ctx->fqctx);
        fq_nmod_poly_init(b, ctx->fqctx);
        fq_nmod_poly_init(g, ctx->fqctx);
        fq_nmod_poly_init(t, ctx->fqctx);
        fq_nmod_poly_init(r, ctx->fqctx);

        _fq_nmod_mpoly_to_fq_nmod_poly_deflate(a, A, v_in_both,
                                                 I->Amin_exp, I->Gstride, ctx);
        _fq_nmod_mpoly_to_fq_nmod_poly_deflate(b, B, v_in_both,
                                                 I->Bmin_exp, I->Gstride, ctx);
        fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
        _fq_nmod_mpoly_from_fq_nmod_poly_inflate(G, I->Gbits, g, v_in_both,
                                                 I->Gmin_exp, I->Gstride, ctx);
        fq_nmod_poly_divrem(t, r, a, g, ctx->fqctx);
        _fq_nmod_mpoly_from_fq_nmod_poly_inflate(Abar, I->Abarbits, t, v_in_both,
                                              I->Abarmin_exp, I->Gstride, ctx);
        fq_nmod_poly_divrem(t, r, b, g, ctx->fqctx);
        _fq_nmod_mpoly_from_fq_nmod_poly_inflate(Bbar, I->Bbarbits, t, v_in_both,
                                              I->Bbarmin_exp, I->Gstride, ctx);
        fq_nmod_poly_clear(a, ctx->fqctx);
        fq_nmod_poly_clear(b, ctx->fqctx);
        fq_nmod_poly_clear(g, ctx->fqctx);
        fq_nmod_poly_clear(t, ctx->fqctx);
        fq_nmod_poly_clear(r, ctx->fqctx);

        goto successful;
    }

    /* check if there is a variable in ess(A) that is not in ess(B) */
    v_in_A_only = -WORD(1);
    v_in_B_only = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (I->Amax_exp[j] > I->Amin_exp[j] && I->Bmax_exp[j] == I->Bmin_exp[j])
        {
            v_in_A_only = j;
            break;
        }
        if (I->Bmax_exp[j] > I->Bmin_exp[j] && I->Amax_exp[j] == I->Amin_exp[j])
        {
            v_in_B_only = j;
            break;
        }
    }
    if (v_in_A_only != -WORD(1))
    {
        success = _try_missing_var(G, I->Gbits,
                                   Abar, I->Abarbits,
                                   Bbar, I->Bbarbits,
                                   v_in_A_only,
                                   A, I->Amin_exp[v_in_A_only],
                                   B, I->Bmin_exp[v_in_A_only],
                                   ctx);
        goto cleanup;
    }
    if (v_in_B_only != -WORD(1))
    {
        success = _try_missing_var(G, I->Gbits,
                                   Bbar, I->Bbarbits,
                                   Abar, I->Abarbits,
                                   v_in_B_only,
                                   B, I->Bmin_exp[v_in_B_only],
                                   A, I->Amin_exp[v_in_B_only],
                                   ctx);
        goto cleanup;
    }

    /* all variable are now either
            missing from both ess(A) and ess(B), or
            present in both ess(A) and ess(B)
        and there are at least two in the latter case
    */

    mpoly_gcd_info_set_estimates_fq_nmod_mpoly(I, A, B, ctx);
    mpoly_gcd_info_set_perm(I, A->length, B->length, ctx->minfo);

    /* everything in I is valid now */

    {
        int gcd_is_trivial = 1;
        int try_a = I->Gdeflate_deg_bounds_are_nice;
        int try_b = I->Gdeflate_deg_bounds_are_nice;
        for (j = 0; j < nvars; j++)
        {
            if (I->Gdeflate_deg_bound[j] != 0)
            {
                gcd_is_trivial = 0;
            }

            if (I->Adeflate_deg[j] != I->Gdeflate_deg_bound[j]
                || I->Amin_exp[j] > I->Bmin_exp[j])
            {
                try_a = 0;
            }

            if (I->Bdeflate_deg[j] != I->Gdeflate_deg_bound[j]
                || I->Bmin_exp[j] > I->Amin_exp[j])
            {
                try_b = 0;
            }
        }

        if (gcd_is_trivial)
            goto calculate_trivial_gcd;

        if ((try_a || try_b) && _try_divides(G, Abar, Bbar,
                                                      A, try_a, B, try_b, ctx))
            goto successful;
    }

    mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
    mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);

    if (I->zippel_time_est < I->brown_time_est)
    {
        if (_try_zippel(G, Abar, Bbar, A, B, I, ctx))
            goto successful;

        if (_try_brown(G, Abar, Bbar, A, B, I, ctx))
            goto successful;
    }
    else
    {
        if (_try_brown(G, Abar, Bbar, A, B, I, ctx))
            goto successful;

        if (_try_zippel(G, Abar, Bbar, A, B, I, ctx))
            goto successful;
    }

    success = 0;
    goto cleanup;

successful:

    success = 1;

cleanup:

    if (success)
    {
        FLINT_ASSERT(G->length > 0);

        fq_nmod_mpoly_repack_bits_inplace(G, I->Gbits, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Abar, I->Abarbits, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Bbar, I->Bbarbits, ctx);

        if (!fq_nmod_is_one(G->coeffs + 0, ctx->fqctx))
        {
            _fq_nmod_vec_scalar_mul_fq_nmod(Abar->coeffs, Abar->coeffs,
                                      Abar->length, G->coeffs + 0, ctx->fqctx);
            _fq_nmod_vec_scalar_mul_fq_nmod(Bbar->coeffs, Bbar->coeffs,
                                      Bbar->length, G->coeffs + 0, ctx->fqctx);
            _fq_nmod_vec_scalar_div_fq_nmod(G->coeffs, G->coeffs,
                                         G->length, G->coeffs + 0, ctx->fqctx);
        }

#if WANT_ASSERT
        fq_nmod_mpoly_mul(T, G, Abar, ctx);
        FLINT_ASSERT(fq_nmod_mpoly_equal(T, Asave, ctx));
        fq_nmod_mpoly_mul(T, G, Bbar, ctx);
        FLINT_ASSERT(fq_nmod_mpoly_equal(T, Bsave, ctx));
#endif
    }

#if WANT_ASSERT
    fq_nmod_mpoly_clear(T, ctx);
    fq_nmod_mpoly_clear(Asave, ctx);
    fq_nmod_mpoly_clear(Bsave, ctx);
#endif

    mpoly_gcd_info_clear(I);

    return success;
}

int fq_nmod_mpoly_gcd_cofactors(
    fq_nmod_mpoly_t G,
    fq_nmod_mpoly_t Abar,
    fq_nmod_mpoly_t Bbar,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t Gbits;
    int success;
    fq_nmod_mpoly_t Anew, Bnew;

    if (A->length == 0)
    {
        if (B->length == 0)
        {
            fq_nmod_mpoly_zero(G, ctx);
            fq_nmod_mpoly_zero(Abar, ctx);
            fq_nmod_mpoly_zero(Bbar, ctx);
            return 1;
        }
        fq_nmod_mpoly_set(G, B, ctx);
        fq_nmod_mpoly_zero(Abar, ctx);
        fq_nmod_mpoly_one(Bbar, ctx);
        if (!fq_nmod_is_one(G->coeffs + 0, ctx->fqctx))
        {
            _fq_nmod_vec_scalar_mul_fq_nmod(Bbar->coeffs, Bbar->coeffs,
                                      Bbar->length, G->coeffs + 0, ctx->fqctx);
            _fq_nmod_vec_scalar_div_fq_nmod(G->coeffs, G->coeffs,
                                         G->length, G->coeffs + 0, ctx->fqctx);
        }
        return 1;
    }

    if (B->length == 0)
    {
        fq_nmod_mpoly_set(G, A, ctx);
        fq_nmod_mpoly_zero(Bbar, ctx);
        fq_nmod_mpoly_one(Abar, ctx);
        if (!fq_nmod_is_one(G->coeffs + 0, ctx->fqctx))
        {
            _fq_nmod_vec_scalar_mul_fq_nmod(Abar->coeffs, Abar->coeffs,
                                      Abar->length, G->coeffs + 0, ctx->fqctx);
            _fq_nmod_vec_scalar_div_fq_nmod(G->coeffs, G->coeffs,
                                         G->length, G->coeffs + 0, ctx->fqctx);
        }
        return 1;
    }

    Gbits = FLINT_MIN(A->bits, B->bits);

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
    {
        /* usual gcd's go right down here */
        return _fq_nmod_mpoly_gcd_cofactors(G, Gbits,
                                      Abar, A->bits, Bbar, B->bits, A, B, ctx);
    }

    fq_nmod_mpoly_init(Anew, ctx);
    fq_nmod_mpoly_init(Bnew, ctx);

    if (A->length == 1)
    {
        _try_monomial_gcd(G, Gbits, Bbar, B->bits, Abar, A->bits, B, A, ctx);
        success = 1;
        goto cleanup;
    }
    else if (B->length == 1)
    {
        _try_monomial_gcd(G, Gbits, Abar, A->bits, Bbar, B->bits, A, B, ctx);
        success = 1;
        goto cleanup;
    }
    else if (_try_monomial_cofactors(G, Gbits, Abar, A->bits, Bbar, B->bits,
                                                                    A, B, ctx))
    {
        success = 1;
        goto cleanup;
    }
    else
    {
        /*
            The gcd calculation is unusual.
            First see if both inputs fit into FLINT_BITS.
            Then, try deflation as a last resort.
        */

        slong k;
        fmpz * Ashift, * Astride;
        fmpz * Bshift, * Bstride;
        fmpz * Gshift, * Gstride;
        const fq_nmod_mpoly_struct * Ause, * Buse;

        Ause = A;
        if (A->bits > FLINT_BITS)
        {
            if (!fq_nmod_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx))
                goto could_not_repack;
            Ause = Anew;
        }

        Buse = B;
        if (B->bits > FLINT_BITS)
        {
            if (!fq_nmod_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx))
                goto could_not_repack;
            Buse = Bnew;
        }

        Gbits = FLINT_MIN(Ause->bits, Buse->bits);
        success = _fq_nmod_mpoly_gcd_cofactors(G, Gbits, Abar, Ause->bits,
                                        Bbar, Buse->bits, Ause, Buse, ctx);
        goto cleanup;


could_not_repack:

        /*
            One of A or B could not be repacked into FLINT_BITS. See if
            they both fit into FLINT_BITS after deflation.
        */

        Ashift  = _fmpz_vec_init(ctx->minfo->nvars);
        Astride = _fmpz_vec_init(ctx->minfo->nvars);
        Bshift  = _fmpz_vec_init(ctx->minfo->nvars);
        Bstride = _fmpz_vec_init(ctx->minfo->nvars);
        Gshift  = _fmpz_vec_init(ctx->minfo->nvars);
        Gstride = _fmpz_vec_init(ctx->minfo->nvars);

        fq_nmod_mpoly_deflation(Ashift, Astride, A, ctx);
        fq_nmod_mpoly_deflation(Bshift, Bstride, B, ctx);
        _fmpz_vec_min(Gshift, Ashift, Bshift, ctx->minfo->nvars);
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_gcd(Gstride + k, Astride + k, Bstride + k);
        }

        success = 0;

        fq_nmod_mpoly_deflate(Anew, A, Ashift, Gstride, ctx);
        if (Anew->bits > FLINT_BITS)
        {
            if (!fq_nmod_mpoly_repack_bits(Anew, Anew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        fq_nmod_mpoly_deflate(Bnew, B, Bshift, Gstride, ctx);
        if (Bnew->bits > FLINT_BITS)
        {
            if (!fq_nmod_mpoly_repack_bits(Bnew, Bnew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        Gbits = FLINT_MIN(Anew->bits, Bnew->bits);
        success = _fq_nmod_mpoly_gcd_cofactors(G, Gbits, Abar, Anew->bits,
                                            Bbar, Bnew->bits, Anew, Bnew, ctx);
        if (!success)
            goto deflate_cleanup;

        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_sub(Ashift + k, Ashift + k, Gshift + k);
            fmpz_sub(Bshift + k, Bshift + k, Gshift + k);
            FLINT_ASSERT(fmpz_sgn(Ashift + k) >= 0);
            FLINT_ASSERT(fmpz_sgn(Bshift + k) >= 0);
        }

        fq_nmod_mpoly_inflate(G, G, Gshift, Gstride, ctx);
        fq_nmod_mpoly_inflate(Abar, Abar, Ashift, Gstride, ctx);
        fq_nmod_mpoly_inflate(Bbar, Bbar, Bshift, Gstride, ctx);

        FLINT_ASSERT(G->length > 0);
        if (!fq_nmod_is_one(G->coeffs + 0, ctx->fqctx))
        {
            _fq_nmod_vec_scalar_mul_fq_nmod(Abar->coeffs, Abar->coeffs,
                                      Abar->length, G->coeffs + 0, ctx->fqctx);
            _fq_nmod_vec_scalar_mul_fq_nmod(Bbar->coeffs, Bbar->coeffs,
                                       Bbar->length, G->coeffs + 0, ctx->fqctx);
            _fq_nmod_vec_scalar_div_fq_nmod(G->coeffs, G->coeffs,
                                         G->length, G->coeffs + 0, ctx->fqctx);
        }

        success = 1;

deflate_cleanup:

        _fmpz_vec_clear(Ashift, ctx->minfo->nvars);
        _fmpz_vec_clear(Astride, ctx->minfo->nvars);
        _fmpz_vec_clear(Bshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Bstride, ctx->minfo->nvars);
        _fmpz_vec_clear(Gshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Gstride, ctx->minfo->nvars);

cleanup:

        fq_nmod_mpoly_clear(Anew, ctx);
        fq_nmod_mpoly_clear(Bnew, ctx);

        return success;
    }
}
