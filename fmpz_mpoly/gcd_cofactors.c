/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


/* from gcd.c */
void fmpz_mpoly_evals(
    nmod_poly_struct * out,
    const int * ignore,
    const fmpz_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles);

void mpoly_gcd_info_set_estimates_fmpz_mpoly(
    mpoly_gcd_info_t I,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles);


/*********************** Easy when B is a monomial ***************************/
static void _try_monomial_gcd(
    fmpz_mpoly_t G, flint_bitcnt_t Gbits,
    fmpz_mpoly_t Abar, flint_bitcnt_t Abarbits,
    fmpz_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t g;
    fmpz * minAfields, * minAdegs, * minBdegs;
    fmpz_mpoly_t _G, _Abar, _Bbar;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length == 1);

    fmpz_mpoly_init(_G, ctx);
    fmpz_mpoly_init(_Abar, ctx);
    fmpz_mpoly_init(_Bbar, ctx);

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

    fmpz_mpoly_fit_length(_G, 1, ctx);
    fmpz_mpoly_fit_bits(_G, Gbits, ctx);
    _G->bits = Gbits;
    mpoly_set_monomial_ffmpz(_G->exps, minBdegs, Gbits, ctx->minfo);

    fmpz_init_set(g, B->coeffs + 0);
    _fmpz_vec_content_chained(g, A->coeffs, A->length);
    fmpz_swap(_G->coeffs + 0, g);
    fmpz_clear(g);

    _fmpz_mpoly_set_length(_G, 1, ctx);

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

    fmpz_mpoly_divides_threaded(_Abar, A, _G, ctx, 0);
    fmpz_mpoly_divides_threaded(_Bbar, B, _G, ctx, 0);

    fmpz_mpoly_swap(G, _G, ctx);
    fmpz_mpoly_swap(Abar, _Abar, ctx);
    fmpz_mpoly_swap(Bbar, _Bbar, ctx);

    fmpz_mpoly_clear(_G, ctx);
    fmpz_mpoly_clear(_Abar, ctx);
    fmpz_mpoly_clear(_Bbar, ctx);
}


/********************** See if cofactors are monomials ***********************/
static int _try_monomial_cofactors(
    fmpz_mpoly_t G, flint_bitcnt_t Gbits,
    fmpz_mpoly_t Abar, flint_bitcnt_t Abarbits,
    fmpz_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong NA, NG;
    slong nvars = ctx->minfo->nvars;
    fmpz * Abarexps, * Bbarexps, * Texps;
    fmpz_t t1, t2;
    fmpz_t gA, gB;
    fmpz_mpoly_t T;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (A->length != B->length)
        return 0;

    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init_set(gA, A->coeffs + 0);
    fmpz_init_set(gB, B->coeffs + 0);

    for (i = A->length - 1; i > 0; i--)
    {
        fmpz_mul(t1, A->coeffs + 0, B->coeffs + i);
        fmpz_mul(t2, B->coeffs + 0, A->coeffs + i);
        success = fmpz_equal(t1, t2);
        if (!success)
            goto cleanup;

        fmpz_gcd(gA, gA, A->coeffs + i);
        fmpz_gcd(gB, gB, B->coeffs + i);
    }

    TMP_START;

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

    /* put A's cofactor coefficient in t1 */
    fmpz_gcd(t2, gA, gB);
    fmpz_divexact(t1, gA, t2);
    if (fmpz_sgn(A->coeffs + 0) < 0)
        fmpz_neg(t1, t1);

    /* put B's cofactor coefficient in t2 */
    fmpz_divexact(gA, A->coeffs + 0, t1);
    fmpz_divexact(t2, B->coeffs + 0, gA);

    fmpz_mpoly_init3(T, A->length, Gbits, ctx);
    NG = mpoly_words_per_exp(Gbits, ctx->minfo);
    NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    T->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ffmpz(Texps, A->exps + NA*i, A->bits, ctx->minfo);
        _fmpz_vec_sub(Texps, Texps, Abarexps, nvars);
        mpoly_set_monomial_ffmpz(T->exps + NG*i, Texps, Gbits, ctx->minfo);
        fmpz_divexact(T->coeffs + i, A->coeffs + i, t1);
    }
    fmpz_mpoly_swap(G, T, ctx);
    fmpz_mpoly_clear(T, ctx);

    fmpz_mpoly_fit_length(Abar, 1, ctx);
    fmpz_mpoly_fit_bits(Abar, Abarbits, ctx);
    Abar->bits = Abarbits;
    mpoly_set_monomial_ffmpz(Abar->exps, Abarexps, Abarbits, ctx->minfo);
    fmpz_swap(Abar->coeffs + 0, t1);
    _fmpz_mpoly_set_length(Abar, 1, ctx);

    fmpz_mpoly_fit_length(Bbar, 1, ctx);
    fmpz_mpoly_fit_bits(Bbar, Bbarbits, ctx);
    Bbar->bits = Bbarbits;
    mpoly_set_monomial_ffmpz(Bbar->exps, Bbarexps, Bbarbits, ctx->minfo);
    fmpz_swap(Bbar->coeffs + 0, t2);
    _fmpz_mpoly_set_length(Bbar, 1, ctx);

    success = 1;

cleanup_tmp:

    for (j = 0; j < nvars; j++)
    {
        fmpz_clear(Abarexps + j);
        fmpz_clear(Bbarexps + j);
        fmpz_clear(Texps + j);
    }

    TMP_END;

cleanup:

    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(gA);
    fmpz_clear(gB);

    return success;
}


/********* Assume B has length one when converted to univar format ***********/
static int _try_missing_var(
    fmpz_mpoly_t G, flint_bitcnt_t Gbits,
    fmpz_mpoly_t Abar, flint_bitcnt_t Abarbits,
    fmpz_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    slong var,
    const fmpz_mpoly_t A, ulong Ashift,
    const fmpz_mpoly_t B, ulong Bshift,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fmpz_mpoly_t tG, tAbar, tBbar;
    fmpz_mpoly_univar_t Ax;

    fmpz_mpoly_init(tG, ctx);
    fmpz_mpoly_init(tAbar, ctx);
    fmpz_mpoly_init(tBbar, ctx);
    fmpz_mpoly_univar_init(Ax, ctx);

#if WANT_ASSERT
    fmpz_mpoly_to_univar(Ax, B, var, ctx);
    FLINT_ASSERT(Ax->length == 1);
#endif

    fmpz_mpoly_to_univar(Ax, A, var, ctx);

    FLINT_ASSERT(Ax->length > 0);
    success = _fmpz_mpoly_gcd(tG, Gbits, B, Ax->coeffs + 0, ctx, NULL, 0);
    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _fmpz_mpoly_gcd(tG, Gbits, tG, Ax->coeffs + i, ctx, NULL, 0);
        if (!success)
            goto cleanup;
    }

    _mpoly_gen_shift_left(tG->exps, tG->bits, tG->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

    success = fmpz_mpoly_divides_threaded(tAbar, A, tG, ctx, 0);
    FLINT_ASSERT(success);
    success = fmpz_mpoly_divides_threaded(tBbar, B, tG, ctx, 0);
    FLINT_ASSERT(success);

    fmpz_mpoly_swap(G, tG, ctx);
    fmpz_mpoly_swap(Abar, tAbar, ctx);
    fmpz_mpoly_swap(Bbar, tBbar, ctx);

    success = 1;

cleanup:

    fmpz_mpoly_clear(tG, ctx);
    fmpz_mpoly_clear(tAbar, ctx);
    fmpz_mpoly_clear(tBbar, ctx);
    fmpz_mpoly_univar_clear(Ax, ctx);

    return success;
}


/******************* Test if B divides A or A divides B **********************/
static int _try_divides(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A, int try_a,
    const fmpz_mpoly_t B, int try_b,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_t cA, cB, cG;
    fmpz_mpoly_t Q;
    fmpz_mpoly_t AA, BB;
    slong AA_alloc, BB_alloc;

    *AA = *A;
    *BB = *B;
    fmpz_init(cA);
    fmpz_init(cB);
    fmpz_init(cG);
    fmpz_mpoly_init(Q, ctx);

    _fmpz_vec_content(cA, A->coeffs, A->length);
    _fmpz_vec_content(cB, B->coeffs, B->length);
    fmpz_gcd(cG, cA, cB);

    AA_alloc = 0;
    if (!fmpz_is_one(cA))
    {
        AA_alloc = A->length;
        AA->coeffs = _fmpz_vec_init(A->length);
        _fmpz_vec_scalar_divexact_fmpz(AA->coeffs, A->coeffs, A->length, cA);
        FLINT_ASSERT(AA_alloc > 0);
    }

    BB_alloc = 0;
    if (!fmpz_is_one(cB))
    {
        BB_alloc = B->length;
        BB->coeffs = _fmpz_vec_init(B->length);
        _fmpz_vec_scalar_divexact_fmpz(BB->coeffs, B->coeffs, B->length, cB);
        FLINT_ASSERT(BB_alloc > 0);
    }

    fmpz_divexact(cA, cA, cG);
    fmpz_divexact(cB, cB, cG);

    if (try_b && fmpz_mpoly_divides_threaded(Q, AA, BB, ctx, 1))
    {
        fmpz_mpoly_scalar_divexact_fmpz(G, B, cB, ctx);
        fmpz_mpoly_swap(Abar, Q, ctx);
        _fmpz_vec_scalar_mul_fmpz(Abar->coeffs, Abar->coeffs, Abar->length, cA);
        fmpz_mpoly_set_fmpz(Bbar, cB, ctx);
        success = 1;
        goto cleanup;
    }

    if (try_a && fmpz_mpoly_divides_threaded(Q, BB, AA, ctx, 1))
    {
        fmpz_mpoly_scalar_divexact_fmpz(G, A, cA, ctx);
        fmpz_mpoly_swap(Bbar, Q, ctx);
        _fmpz_vec_scalar_mul_fmpz(Bbar->coeffs, Bbar->coeffs, Bbar->length, cB);
        fmpz_mpoly_set_fmpz(Abar, cA, ctx);
        success = 1;
        goto cleanup;
    }

    success = 0;

cleanup:

    if (AA_alloc > 0)
        _fmpz_vec_clear(AA->coeffs, AA_alloc);

    if (BB_alloc > 0)
        _fmpz_vec_clear(BB->coeffs, BB_alloc);

    fmpz_mpoly_clear(Q, ctx);
    fmpz_clear(cA);
    fmpz_clear(cB);
    fmpz_clear(cG);

    return success;
}


/********************** Hit A and B with zippel ******************************/
static int _try_zippel(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    mpoly_zipinfo_t zinfo;
    flint_bitcnt_t wbits;
    flint_rand_t randstate;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    fmpz_mpoly_t Ac, Bc, Gc, Abarc, Bbarc;

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

    /* uctx is context for Z[y_1,...,y_{m-1}]*/
    fmpz_mpoly_ctx_init(uctx, m - 1, ORD_LEX);

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

    fmpz_mpolyu_init(Au, wbits, uctx);
    fmpz_mpolyu_init(Bu, wbits, uctx);
    fmpz_mpolyu_init(Gu, wbits, uctx);
    fmpz_mpolyu_init(Abaru, wbits, uctx);
    fmpz_mpolyu_init(Bbaru, wbits, uctx);
    fmpz_mpoly_init3(Ac, 0, wbits, uctx);
    fmpz_mpoly_init3(Bc, 0, wbits, uctx);
    fmpz_mpoly_init3(Gc, 0, wbits, uctx);
    fmpz_mpoly_init3(Abarc, 0, wbits, uctx);
    fmpz_mpoly_init3(Bbarc, 0, wbits, uctx);

    fmpz_mpoly_to_mpolyu_perm_deflate(Au, uctx, A, ctx, zinfo->perm,
                                I->Amin_exp, I->Gstride, I->Amax_exp, NULL, 0);
    fmpz_mpoly_to_mpolyu_perm_deflate(Bu, uctx, B, ctx, zinfo->perm,
                                I->Bmin_exp, I->Gstride, I->Bmax_exp, NULL, 0);

    FLINT_ASSERT(Au->bits == wbits);
    FLINT_ASSERT(Bu->bits == wbits);
    FLINT_ASSERT(Au->length > 1);
    FLINT_ASSERT(Bu->length > 1);

    success = fmpz_mpolyu_content_mpoly(Ac, Au, uctx, NULL, 0);
    success = success && fmpz_mpolyu_content_mpoly(Bc, Bu, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_divexact_mpoly_inplace(Au, Ac, uctx);
    fmpz_mpolyu_divexact_mpoly_inplace(Bu, Bc, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */
    success = fmpz_mpolyu_gcdm_zippel(Gu, Abaru, Bbaru, Au, Bu,
                                                       uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    success = _fmpz_mpoly_gcd_cofactors(Gc, wbits, Abarc, wbits, Bbarc, wbits,
                                                        Ac, Bc, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);
    fmpz_mpolyu_mul_mpoly_inplace(Abaru, Abarc, uctx);
    fmpz_mpolyu_mul_mpoly_inplace(Bbaru, Bbarc, uctx);

    fmpz_mpoly_from_mpolyu_perm_inflate(G, I->Gbits, ctx, Gu, uctx,
                                         zinfo->perm, I->Gmin_exp, I->Gstride);
    fmpz_mpoly_from_mpolyu_perm_inflate(Abar, I->Abarbits, ctx, Abaru, uctx,
                                      zinfo->perm, I->Abarmin_exp, I->Gstride);
    fmpz_mpoly_from_mpolyu_perm_inflate(Bbar, I->Bbarbits, ctx, Bbaru, uctx,
                                      zinfo->perm, I->Bbarmin_exp, I->Gstride);
    success = 1;

cleanup:

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Gu, uctx);
    fmpz_mpolyu_clear(Abaru, uctx);
    fmpz_mpolyu_clear(Bbaru, uctx);
    fmpz_mpoly_clear(Ac, uctx);
    fmpz_mpoly_clear(Bc, uctx);
    fmpz_mpoly_clear(Gc, uctx);
    fmpz_mpoly_clear(Abarc, uctx);
    fmpz_mpoly_clear(Bbarc, uctx);

    fmpz_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    return success;
}


/************************ Hit A and B with bma *******************************/
typedef struct
{
    const fmpz_mpoly_struct * P;
    fmpz_mpoly_struct * Pcontent;
    fmpz_mpolyu_struct * Puu;
    const slong * perm;
    const ulong * shift, * stride, * maxexps;
    const fmpz_mpoly_ctx_struct * ctx;
    const fmpz_mpoly_ctx_struct * uctx;
    const thread_pool_handle * handles;
    slong num_handles;
    int success;
}
_convertuu_arg_struct;

typedef _convertuu_arg_struct _convertuu_arg_t[1];

static void _worker_convertuu(void * varg)
{
    _convertuu_arg_struct * arg = (_convertuu_arg_struct *) varg;

    fmpz_mpoly_to_mpolyuu_perm_deflate(arg->Puu, arg->uctx, arg->P, arg->ctx,
                             arg->perm, arg->shift, arg->stride, arg->maxexps,
                                               arg->handles, arg->num_handles);

    arg->success = fmpz_mpolyu_content_mpoly(arg->Pcontent, arg->Puu,
                                                           arg->uctx, NULL, 0);
    if (arg->success)
    {
        fmpz_mpolyu_divexact_mpoly_inplace(arg->Puu, arg->Pcontent, arg->uctx);
    }
}

static int _try_bma(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Auu, Buu, Guu, Abaruu, Bbaruu;
    fmpz_mpoly_t Ac, Bc, Gc, Abarc, Bbarc, Gamma;
    slong max_minor_degree;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!I->can_use_bma)
        return 0;

    FLINT_ASSERT(m >= WORD(3));

    /* uctx is context for Z[y_2,...,y_{m - 1}]*/
    fmpz_mpoly_ctx_init(uctx, m - 2, ORD_LEX);

    max_minor_degree = 0;
    for (i = 2; i < m; i++)
    {
        k = I->bma_perm[i];
        max_minor_degree = FLINT_MAX(max_minor_degree, I->Adeflate_deg[k]);
        max_minor_degree = FLINT_MAX(max_minor_degree, I->Bdeflate_deg[k]);
    }

    wbits = 1 + FLINT_BIT_COUNT(max_minor_degree);
    wbits = FLINT_MAX(MPOLY_MIN_BITS, wbits);
    wbits = mpoly_fix_bits(wbits, uctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    fmpz_mpolyu_init(Auu, wbits, uctx);
    fmpz_mpolyu_init(Buu, wbits, uctx);
    fmpz_mpolyu_init(Guu, wbits, uctx);
    fmpz_mpolyu_init(Abaruu, wbits, uctx);
    fmpz_mpolyu_init(Bbaruu, wbits, uctx);
    fmpz_mpoly_init3(Ac, 0, wbits, uctx);
    fmpz_mpoly_init3(Bc, 0, wbits, uctx);
    fmpz_mpoly_init3(Gc, 0, wbits, uctx);
    fmpz_mpoly_init3(Abarc, 0, wbits, uctx);
    fmpz_mpoly_init3(Bbarc, 0, wbits, uctx);
    fmpz_mpoly_init3(Gamma, 0, wbits, uctx);

    /* convert to bivariate format and remove content from A and B */
    if (num_handles > 0)
    {
        slong s = mpoly_divide_threads(num_handles, A->length, B->length);
        _convertuu_arg_t arg;

        FLINT_ASSERT(s >= 0);
        FLINT_ASSERT(s < num_handles);

        arg->ctx = ctx;
        arg->uctx = uctx;
        arg->P = B;
        arg->Puu = Buu;
        arg->Pcontent = Bc;
        arg->perm = I->bma_perm;
        arg->shift = I->Bmin_exp;
        arg->stride = I->Gstride;
        arg->maxexps = I->Bmax_exp;
        arg->handles = handles + (s + 1);
        arg->num_handles = num_handles - (s + 1);

        thread_pool_wake(global_thread_pool, handles[s], _worker_convertuu, arg);

        fmpz_mpoly_to_mpolyuu_perm_deflate(Auu, uctx, A, ctx,
                          I->bma_perm, I->Amin_exp, I->Gstride, I->Amax_exp,
                                                               handles + 0, s);
        success = fmpz_mpolyu_content_mpoly(Ac, Auu, uctx, handles + 0, s);
        if (success)
        {
            fmpz_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);
        }

        thread_pool_wait(global_thread_pool, handles[s]);

        success = success && arg->success;
        if (!success)
            goto cleanup;
    }
    else
    {
        fmpz_mpoly_to_mpolyuu_perm_deflate(Auu, uctx, A, ctx,
                   I->bma_perm, I->Amin_exp, I->Gstride, I->Amax_exp, NULL, 0);
        fmpz_mpoly_to_mpolyuu_perm_deflate(Buu, uctx, B, ctx,
                   I->bma_perm, I->Bmin_exp, I->Gstride, I->Bmax_exp, NULL, 0);

        success = fmpz_mpolyu_content_mpoly(Ac, Auu, uctx, NULL, 0);
        success = success && fmpz_mpolyu_content_mpoly(Bc, Buu, uctx, NULL, 0);
        if (!success)
            goto cleanup;

        fmpz_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);
        fmpz_mpolyu_divexact_mpoly_inplace(Buu, Bc, uctx);
    }

    FLINT_ASSERT(Auu->length > 1);
    FLINT_ASSERT(Buu->length > 1);

    _fmpz_mpoly_gcd(Gamma, wbits, Auu->coeffs + 0,
                                  Buu->coeffs + 0, uctx, handles, num_handles);
    if (!success)
        goto cleanup;

    success = (num_handles > 0)
           ? fmpz_mpolyuu_gcd_berlekamp_massey_threaded(Guu, Abaruu, Bbaruu,
                                  Auu, Buu, Gamma, uctx, handles, num_handles)
           : fmpz_mpolyuu_gcd_berlekamp_massey(Guu, Abaruu, Bbaruu,
                                                        Auu, Buu, Gamma, uctx);
    if (!success)
        goto cleanup;

    success = _fmpz_mpoly_gcd_cofactors(Gc, wbits, Abarc, wbits, Bbarc, wbits,
                                           Ac, Bc, uctx, handles, num_handles);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_mul_mpoly_inplace(Guu, Gc, uctx);
    fmpz_mpolyu_mul_mpoly_inplace(Abaruu, Abarc, uctx);
    fmpz_mpolyu_mul_mpoly_inplace(Bbaruu, Bbarc, uctx);

    fmpz_mpoly_from_mpolyuu_perm_inflate(G, I->Gbits, ctx, Guu, uctx,
                                         I->bma_perm, I->Gmin_exp, I->Gstride);
    fmpz_mpoly_from_mpolyuu_perm_inflate(Abar, I->Abarbits, ctx, Abaruu, uctx,
                                      I->bma_perm, I->Abarmin_exp, I->Gstride);
    fmpz_mpoly_from_mpolyuu_perm_inflate(Bbar, I->Bbarbits, ctx, Bbaruu, uctx,
                                      I->bma_perm, I->Bbarmin_exp, I->Gstride);
    success = 1;

cleanup:

    fmpz_mpolyu_clear(Auu, uctx);
    fmpz_mpolyu_clear(Buu, uctx);
    fmpz_mpolyu_clear(Guu, uctx);
    fmpz_mpolyu_clear(Abaruu, uctx);
    fmpz_mpolyu_clear(Bbaruu, uctx);
    fmpz_mpoly_clear(Ac, uctx);
    fmpz_mpoly_clear(Bc, uctx);
    fmpz_mpoly_clear(Gc, uctx);
    fmpz_mpoly_clear(Abarc, uctx);
    fmpz_mpoly_clear(Bbarc, uctx);
    fmpz_mpoly_clear(Gamma, uctx);

    fmpz_mpoly_ctx_clear(uctx);

    return success;
}


/*********************** Hit A and B with brown ******************************/
typedef struct
{
    fmpz_mpoly_struct * Pl;
    const fmpz_mpoly_ctx_struct * lctx;
    const fmpz_mpoly_struct * P;
    const fmpz_mpoly_ctx_struct * ctx;
    const slong * perm;
    const ulong * shift, * stride, * maxexps;
    const thread_pool_handle * handles;
    slong num_handles;
}
_convertl_arg_struct;

typedef _convertl_arg_struct _convertl_arg_t[1];

static void _worker_convertu(void * varg)
{
    _convertl_arg_struct * arg = (_convertl_arg_struct *) varg;

    fmpz_mpoly_to_mpoly_perm_deflate(arg->Pl, arg->lctx, arg->P, arg->ctx,
                                         arg->perm, arg->shift, arg->stride,
                                               arg->handles, arg->num_handles);
}

static int _try_brown(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles, slong num_handles)
{
    int success;
    slong m = I->mvars;
    flint_bitcnt_t wbits;
    fmpz_mpoly_ctx_t lctx;
    fmpz_mpoly_t Al, Bl, Gl, Abarl, Bbarl;

    if (!I->can_use_brown)
        return 0;

    FLINT_ASSERT(m >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    wbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_ctx_init(lctx, m, ORD_LEX);
    fmpz_mpoly_init3(Al, 0, wbits, lctx);
    fmpz_mpoly_init3(Bl, 0, wbits, lctx);
    fmpz_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Bbarl, 0, wbits, lctx);

    if (num_handles > 0)
    {
        slong s = mpoly_divide_threads(num_handles, A->length, B->length);
        _convertl_arg_t arg;

        FLINT_ASSERT(s >= 0);
        FLINT_ASSERT(s < num_handles);

        arg->Pl = Bl;
        arg->lctx = lctx;
        arg->P = B;
        arg->ctx = ctx;
        arg->perm = I->brown_perm;
        arg->shift = I->Bmin_exp;
        arg->stride = I->Gstride;
        arg->maxexps = I->Bmax_exp;
        arg->handles = handles + (s + 1);
        arg->num_handles = num_handles - (s + 1);

        thread_pool_wake(global_thread_pool, handles[s], _worker_convertu, arg);

        fmpz_mpoly_to_mpoly_perm_deflate(Al, lctx, A, ctx,
                                    I->brown_perm, I->Amin_exp, I->Gstride,
                                                               handles + 0, s);

        thread_pool_wait(global_thread_pool, handles[s]);
    }
    else
    {
        fmpz_mpoly_to_mpoly_perm_deflate(Al, lctx, A, ctx,
                              I->brown_perm, I->Amin_exp, I->Gstride, NULL, 0);
        fmpz_mpoly_to_mpoly_perm_deflate(Bl, lctx, B, ctx,
                              I->brown_perm, I->Bmin_exp, I->Gstride, NULL, 0);
    }

    FLINT_ASSERT(Al->bits == wbits);
    FLINT_ASSERT(Bl->bits == wbits);
    FLINT_ASSERT(Al->length > 1);
    FLINT_ASSERT(Bl->length > 1);

    success = (num_handles > 0)
           ? fmpz_mpolyl_gcd_brown_threaded(Gl, Abarl, Bbarl, Al, Bl, lctx, I,
                                                         handles, num_handles)
           : fmpz_mpolyl_gcd_brown(Gl, Abarl, Bbarl, Al, Bl, lctx, I);

    if (!success)
        goto cleanup;

    fmpz_mpoly_from_mpoly_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                       I->brown_perm, I->Gmin_exp, I->Gstride);
    fmpz_mpoly_from_mpoly_perm_inflate(Abar, I->Abarbits, ctx, Abarl, lctx,
                                    I->brown_perm, I->Abarmin_exp, I->Gstride);
    fmpz_mpoly_from_mpoly_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarl, lctx,
                                    I->brown_perm, I->Bbarmin_exp, I->Gstride);

    success = 1;

cleanup:

    fmpz_mpoly_clear(Al, lctx);
    fmpz_mpoly_clear(Bl, lctx);
    fmpz_mpoly_clear(Gl, lctx);
    fmpz_mpoly_clear(Abarl, lctx);
    fmpz_mpoly_clear(Bbarl, lctx);
    fmpz_mpoly_ctx_clear(lctx);

    return success;
}


/*
    The function must pack successful answers into the corresponding bits.
    Both A and B have to be packed into bits <= FLINT_BITS.
    return is 1 for success, 0 for failure.
*/
int _fmpz_mpoly_gcd_cofactors(
    fmpz_mpoly_t G, flint_bitcnt_t Gbits,
    fmpz_mpoly_t Abar, flint_bitcnt_t Abarbits,
    fmpz_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int success;
    slong v_in_both;
    slong v_in_either;
    slong v_in_A_only;
    slong v_in_B_only;
    slong j;
    slong nvars = ctx->minfo->nvars;
    mpoly_gcd_info_t I;

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
            I->Abarmin_exp
            I->Bbarmin_exp
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
        fmpz_t cG;

calculate_trivial_gcd:

        fmpz_init(cG);
        _fmpz_vec_content_chained(cG, A->coeffs, A->length);
        _fmpz_vec_content_chained(cG, B->coeffs, B->length);

        if (Abar == B && Bbar == A)
        {
            fmpz_mpoly_scalar_divexact_fmpz(Abar, B, cG, ctx);
            fmpz_mpoly_scalar_divexact_fmpz(Bbar, A, cG, ctx);
            fmpz_mpoly_swap(Abar, Bbar, ctx);
        }
        else if (Abar == B && Bbar != A)
        {
            fmpz_mpoly_scalar_divexact_fmpz(Bbar, B, cG, ctx);
            fmpz_mpoly_scalar_divexact_fmpz(Abar, A, cG, ctx);
        }
        else
        {
            fmpz_mpoly_scalar_divexact_fmpz(Abar, A, cG, ctx);
            fmpz_mpoly_scalar_divexact_fmpz(Bbar, B, cG, ctx);
        }

        fmpz_mpoly_fit_length(G, 1, ctx);
        fmpz_mpoly_fit_bits(G, I->Gbits, ctx);
        G->bits = I->Gbits;
        mpoly_set_monomial_ui(G->exps, I->Gmin_exp, I->Gbits, ctx->minfo);
        fmpz_swap(G->coeffs + 0, cG);
        _fmpz_mpoly_set_length(G, 1, ctx);

        mpoly_monomials_shift_right_ui(Abar->exps, Abar->bits, Abar->length,
                                                      I->Gmin_exp, ctx->minfo);
        mpoly_monomials_shift_right_ui(Bbar->exps, Bbar->bits, Bbar->length,
                                                      I->Gmin_exp, ctx->minfo);

        fmpz_clear(cG);

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
        fmpz_poly_t a, b, g, t;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        fmpz_poly_init(t);

        _fmpz_mpoly_to_fmpz_poly_deflate(a, A, v_in_both,
                                                 I->Amin_exp, I->Gstride, ctx);
        _fmpz_mpoly_to_fmpz_poly_deflate(b, B, v_in_both,
                                                 I->Bmin_exp, I->Gstride, ctx);
        fmpz_poly_gcd(g, a, b);
        _fmpz_mpoly_from_fmpz_poly_inflate(G, I->Gbits, g, v_in_both,
                                                 I->Gmin_exp, I->Gstride, ctx);
        fmpz_poly_div(t, a, g);
        _fmpz_mpoly_from_fmpz_poly_inflate(Abar, I->Abarbits, t, v_in_both,
                                              I->Abarmin_exp, I->Gstride, ctx);
        fmpz_poly_div(t, b, g);
        _fmpz_mpoly_from_fmpz_poly_inflate(Bbar, I->Bbarbits, t, v_in_both,
                                              I->Bbarmin_exp, I->Gstride, ctx);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);
        fmpz_poly_clear(t);

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

    /*
        all variable are now either
            missing from both ess(A) and ess(B), or
            present in both ess(A) and ess(B)
        and there are at least two in the latter case
    */

    mpoly_gcd_info_set_estimates_fmpz_mpoly(I, A, B, ctx, handles, num_handles);
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
        {
            goto successful;
        }
    }

    mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
    mpoly_gcd_info_measure_bma(I, A->length, B->length, ctx->minfo);

    if (I->mvars == 2)
    {
        /* TODO: bivariate heuristic here */

        if (_try_brown(G, Abar, Bbar, A, B, I, ctx, handles, num_handles))
            goto successful;
    }
    else if (I->can_use_brown && I->can_use_bma
            && I->bma_time_est < I->brown_time_est
            && (I->mvars*(I->Adensity + I->Bdensity) < 1
                || I->bma_time_est < 0.01*I->brown_time_est))
    {
        if (_try_bma(G, Abar, Bbar, A, B, I, ctx, handles, num_handles))
            goto successful;

        if (_try_brown(G, Abar, Bbar, A, B, I, ctx, handles, num_handles))
            goto successful;
    }
    else
    {
        if (_try_brown(G, Abar, Bbar, A, B, I, ctx, handles, num_handles))
            goto successful;

        if (_try_bma(G, Abar, Bbar, A, B, I, ctx, handles, num_handles))
            goto successful;
    }

    mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);

    if (_try_zippel(G, Abar, Bbar, A, B, I, ctx))
        goto successful;

    success = 0;
    goto cleanup;

successful:

    success = 1;

cleanup:

    if (success)
    {
        FLINT_ASSERT(G->length > 0);

        fmpz_mpoly_repack_bits_inplace(G, I->Gbits, ctx);
        fmpz_mpoly_repack_bits_inplace(Abar, I->Abarbits, ctx);
        fmpz_mpoly_repack_bits_inplace(Bbar, I->Bbarbits, ctx);

        if (fmpz_sgn(G->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, G, ctx);
            fmpz_mpoly_neg(Abar, Abar, ctx);
            fmpz_mpoly_neg(Bbar, Bbar, ctx);
        }
    }

    mpoly_gcd_info_clear(I);

    return success;
}


int fmpz_mpoly_gcd_cofactors_threaded(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    slong thread_limit)
{
    slong i;
    flint_bitcnt_t Gbits;
    int success;
    thread_pool_handle * handles;
    slong num_handles;
    fmpz_mpoly_t Anew, Bnew;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (B->length == 0)
        {
            fmpz_mpoly_zero(G, ctx);
            fmpz_mpoly_zero(Abar, ctx);
            fmpz_mpoly_zero(Bbar, ctx);
            return 1;
        }
        fmpz_mpoly_set(G, B, ctx);
        fmpz_mpoly_zero(Abar, ctx);
        fmpz_mpoly_one(Bbar, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, G, ctx);
            fmpz_mpoly_neg(Bbar, Bbar, ctx);
        }
        return 1;
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        fmpz_mpoly_set(G, A, ctx);
        fmpz_mpoly_zero(Bbar, ctx);
        fmpz_mpoly_one(Abar, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, G, ctx);
            fmpz_mpoly_neg(Abar, Abar, ctx);
        }
        return 1;
    }

    Gbits = FLINT_MIN(A->bits, B->bits);

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
    {
        /* usual gcd's go right down here */

        /* get workers */
        handles = NULL;
        num_handles = 0;
        if (global_thread_pool_initialized)
        {
            slong max_num_handles = thread_pool_get_size(global_thread_pool);
            max_num_handles = FLINT_MIN(thread_limit - 1, max_num_handles);
            if (max_num_handles > 0)
            {
                handles = (thread_pool_handle *) flint_malloc(
                                   max_num_handles*sizeof(thread_pool_handle));
                num_handles = thread_pool_request(global_thread_pool,
                                                     handles, max_num_handles);
            }
        }

        success = _fmpz_mpoly_gcd_cofactors(G, Gbits,
                Abar, A->bits, Bbar, B->bits, A, B, ctx, handles, num_handles);

        for (i = 0; i < num_handles; i++)
        {
            thread_pool_give_back(global_thread_pool, handles[i]);
        }
        if (handles)
        {
            flint_free(handles);
        }

        return success;
    }

    fmpz_mpoly_init(Anew, ctx);
    fmpz_mpoly_init(Bnew, ctx);

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
        const fmpz_mpoly_struct * Ause, * Buse;

        Ause = A;
        if (A->bits > FLINT_BITS)
        {
            if (!fmpz_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx))
                goto could_not_repack;
            Ause = Anew;
        }

        Buse = B;
        if (B->bits > FLINT_BITS)
        {
            if (!fmpz_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx))
                goto could_not_repack;
            Buse = Bnew;
        }

        Gbits = FLINT_MIN(Ause->bits, Buse->bits);
        success = _fmpz_mpoly_gcd_cofactors(G, Gbits, Abar, Ause->bits,
                                   Bbar, Buse->bits, Ause, Buse, ctx, NULL, 0);
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

        fmpz_mpoly_deflation(Ashift, Astride, A, ctx);
        fmpz_mpoly_deflation(Bshift, Bstride, B, ctx);
        _fmpz_vec_min(Gshift, Ashift, Bshift, ctx->minfo->nvars);
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_gcd(Gstride + k, Astride + k, Bstride + k);
        }

        success = 0;

        fmpz_mpoly_deflate(Anew, A, Ashift, Gstride, ctx);
        if (Anew->bits > FLINT_BITS)
        {
            if (!fmpz_mpoly_repack_bits(Anew, Anew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        fmpz_mpoly_deflate(Bnew, B, Bshift, Gstride, ctx);
        if (Bnew->bits > FLINT_BITS)
        {
            if (!fmpz_mpoly_repack_bits(Bnew, Bnew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        Gbits = FLINT_MIN(Anew->bits, Bnew->bits);
        success = _fmpz_mpoly_gcd_cofactors(G, Gbits, Abar, Anew->bits,
                                   Bbar, Bnew->bits, Anew, Bnew, ctx, NULL, 0);

        if (!success)
            goto deflate_cleanup;

        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_sub(Ashift + k, Ashift + k, Gshift + k);
            fmpz_sub(Bshift + k, Bshift + k, Gshift + k);
            FLINT_ASSERT(fmpz_sgn(Ashift + k) >= 0);
            FLINT_ASSERT(fmpz_sgn(Bshift + k) >= 0);
        }

        fmpz_mpoly_inflate(G, G, Gshift, Gstride, ctx);
        fmpz_mpoly_inflate(Abar, Abar, Ashift, Gstride, ctx);
        fmpz_mpoly_inflate(Bbar, Bbar, Bshift, Gstride, ctx);

        /* inflation may have changed the lc */
        FLINT_ASSERT(G->length > 0);
        if (fmpz_sgn(G->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, G, ctx);
            fmpz_mpoly_neg(Abar, Abar, ctx);
            fmpz_mpoly_neg(Bbar, Bbar, ctx);
        }

deflate_cleanup:

        _fmpz_vec_clear(Ashift, ctx->minfo->nvars);
        _fmpz_vec_clear(Astride, ctx->minfo->nvars);
        _fmpz_vec_clear(Bshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Bstride, ctx->minfo->nvars);
        _fmpz_vec_clear(Gshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Gstride, ctx->minfo->nvars);
    }

cleanup:

    fmpz_mpoly_clear(Anew, ctx);
    fmpz_mpoly_clear(Bnew, ctx);

    return success;

}

int fmpz_mpoly_gcd_cofactors(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_gcd_cofactors_threaded(G, Abar, Bbar, A, B,
                                              ctx, MPOLY_DEFAULT_THREAD_LIMIT);
}

