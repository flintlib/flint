/*
    Copyright (C) 2018, 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


/* from gcd.c */
void nmod_mpoly_evals(
    nmod_poly_struct * out,
    const int * ignore,
    const nmod_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles);

void mpoly_gcd_info_set_estimates_nmod_mpoly(
    mpoly_gcd_info_t I,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles);


/*********************** Easy when B is a monomial ***************************/
static void _try_monomial_gcd(
    nmod_mpoly_t G, flint_bitcnt_t Gbits,
    nmod_mpoly_t Abar, flint_bitcnt_t Abarbits,
    nmod_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * minAfields, * minAdegs, * minBdegs;
    nmod_mpoly_t _G, _Abar, _Bbar;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length == 1);

    nmod_mpoly_init(_G, ctx);
    nmod_mpoly_init(_Abar, ctx);
    nmod_mpoly_init(_Bbar, ctx);

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

    nmod_mpoly_fit_length(_G, 1, ctx);
    nmod_mpoly_fit_bits(_G, Gbits, ctx);
    _G->bits = Gbits;
    mpoly_set_monomial_ffmpz(_G->exps, minBdegs, Gbits, ctx->minfo);
    _G->coeffs[0] = UWORD(1);
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

    _nmod_mpoly_divides_threaded_pool(_Abar, A, _G, ctx, NULL, 0);
    _nmod_mpoly_divides_threaded_pool(_Bbar, B, _G, ctx, NULL, 0);

    nmod_mpoly_swap(G, _G, ctx);
    nmod_mpoly_swap(Abar, _Abar, ctx);
    nmod_mpoly_swap(Bbar, _Bbar, ctx);

    nmod_mpoly_clear(_G, ctx);
    nmod_mpoly_clear(_Abar, ctx);
    nmod_mpoly_clear(_Bbar, ctx);
}


/********************** See if cofactors are monomials ***********************/
static int _try_monomial_cofactors(
    nmod_mpoly_t G, flint_bitcnt_t Gbits,
    nmod_mpoly_t Abar, flint_bitcnt_t Abarbits,
    nmod_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong NA, NG;
    slong nvars = ctx->minfo->nvars;
    fmpz * Abarexps, * Bbarexps, * Texps;
    mp_limb_t a0, b0, a0inv;
    nmod_mpoly_t T;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (A->length != B->length)
        return 0;

    a0 = A->coeffs[0];
    b0 = B->coeffs[0];

    for (i = A->length - 1; i > 0; i--)
    {
        success = (nmod_mul(a0, B->coeffs[i], ctx->ffinfo->mod)
                == nmod_mul(b0, A->coeffs[i], ctx->ffinfo->mod));
        if (!success)
            goto cleanup;
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

    nmod_mpoly_init3(T, A->length, Gbits, ctx);
    NG = mpoly_words_per_exp(Gbits, ctx->minfo);
    NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    a0inv = nmod_inv(A->coeffs[0], ctx->ffinfo->mod);
    T->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ffmpz(Texps, A->exps + NA*i, A->bits, ctx->minfo);
        _fmpz_vec_sub(Texps, Texps, Abarexps, nvars);
        mpoly_set_monomial_ffmpz(T->exps + NG*i, Texps, Gbits, ctx->minfo);
        T->coeffs[i] = nmod_mul(A->coeffs[i], a0inv, ctx->ffinfo->mod);
    }
    nmod_mpoly_swap(G, T, ctx);
    nmod_mpoly_clear(T, ctx);

    nmod_mpoly_fit_length(Abar, 1, ctx);
    nmod_mpoly_fit_bits(Abar, Abarbits, ctx);
    Abar->bits = Abarbits;
    mpoly_set_monomial_ffmpz(Abar->exps, Abarexps, Abarbits, ctx->minfo);
    Abar->coeffs[0] = a0;
    _nmod_mpoly_set_length(Abar, 1, ctx);

    nmod_mpoly_fit_length(Bbar, 1, ctx);
    nmod_mpoly_fit_bits(Bbar, Bbarbits, ctx);
    Bbar->bits = Bbarbits;
    mpoly_set_monomial_ffmpz(Bbar->exps, Bbarexps, Bbarbits, ctx->minfo);
    Bbar->coeffs[0] = b0;
    _nmod_mpoly_set_length(Bbar, 1, ctx);

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

    return success;
}


/********* Assume B has length one when converted to univar format ***********/
static int _try_missing_var(
    nmod_mpoly_t G, flint_bitcnt_t Gbits,
    nmod_mpoly_t Abar, flint_bitcnt_t Abarbits,
    nmod_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    slong var,
    const nmod_mpoly_t A, ulong Ashift,
    const nmod_mpoly_t B, ulong Bshift,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    nmod_mpoly_t tG, tAbar, tBbar;
    nmod_mpoly_univar_t Ax;

    nmod_mpoly_init(tG, ctx);
    nmod_mpoly_init(tAbar, ctx);
    nmod_mpoly_init(tBbar, ctx);
    nmod_mpoly_univar_init(Ax, ctx);

#if WANT_ASSERT
    nmod_mpoly_to_univar(Ax, B, var, ctx);
    FLINT_ASSERT(Ax->length == 1);
#endif

    nmod_mpoly_to_univar(Ax, A, var, ctx);

    FLINT_ASSERT(Ax->length > 0);
    success = _nmod_mpoly_gcd_threaded_pool(tG, Gbits, B, Ax->coeffs + 0,
                                                                 ctx, NULL, 0);
    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _nmod_mpoly_gcd_threaded_pool(tG, Gbits, tG, Ax->coeffs + i,
                                                                 ctx, NULL, 0);
        if (!success)
            goto cleanup;
    }

    _mpoly_gen_shift_left(tG->exps, tG->bits, tG->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

    success = _nmod_mpoly_divides_threaded_pool(tAbar, A, tG, ctx, NULL, 0);
    FLINT_ASSERT(success);
    success = _nmod_mpoly_divides_threaded_pool(tBbar, B, tG, ctx, NULL, 0);
    FLINT_ASSERT(success);

    nmod_mpoly_swap(G, tG, ctx);
    nmod_mpoly_swap(Abar, tAbar, ctx);
    nmod_mpoly_swap(Bbar, tBbar, ctx);

    success = 1;

cleanup:

    nmod_mpoly_clear(tG, ctx);
    nmod_mpoly_clear(tAbar, ctx);
    nmod_mpoly_clear(tBbar, ctx);
    nmod_mpoly_univar_clear(Ax, ctx);

    return success;
}


/******************* Test if B divides A or A divides B **********************/
/*
    Test if B divides A or A divides B
        TODO: incorporate deflation
*/
static int _try_divides(
    nmod_mpoly_t G,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t A, int try_a,
    const nmod_mpoly_t B, int try_b,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int success;
    nmod_mpoly_t Q;

    nmod_mpoly_init(Q, ctx);

    if (try_b && _nmod_mpoly_divides_threaded_pool(Q, A, B,
                                                    ctx, handles, num_handles))
    {
        nmod_mpoly_set(G, B, ctx);
        nmod_mpoly_swap(Abar, Q, ctx);
        nmod_mpoly_one(Bbar, ctx);
        success = 1;
        goto cleanup;
    }

    if (try_a && _nmod_mpoly_divides_threaded_pool(Q, B, A,
                                                    ctx, handles, num_handles))
    {
        nmod_mpoly_set(G, A, ctx);
        nmod_mpoly_swap(Bbar, Q, ctx);
        nmod_mpoly_one(Abar, ctx);
        success = 1;
        goto cleanup;
    }

    success = 0;

cleanup:

    nmod_mpoly_clear(Q, ctx);

    return success;
}


/********************** Hit A and B with zippel ******************************/
static int _try_zippel(
    nmod_mpoly_t G,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    mpoly_zipinfo_t zinfo;
    flint_bitcnt_t wbits;
    flint_rand_t randstate;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    nmod_mpoly_t Ac, Bc, Gc, Abarc, Bbarc;

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
    nmod_mpoly_ctx_init(uctx, m - 1, ORD_LEX, ctx->ffinfo->mod.n);

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

    nmod_mpolyu_init(Au, wbits, uctx);
    nmod_mpolyu_init(Bu, wbits, uctx);
    nmod_mpolyu_init(Gu, wbits, uctx);
    nmod_mpolyu_init(Abaru, wbits, uctx);
    nmod_mpolyu_init(Bbaru, wbits, uctx);
    nmod_mpoly_init3(Ac, 0, wbits, uctx);
    nmod_mpoly_init3(Bc, 0, wbits, uctx);
    nmod_mpoly_init3(Gc, 0, wbits, uctx);
    nmod_mpoly_init3(Abarc, 0, wbits, uctx);
    nmod_mpoly_init3(Bbarc, 0, wbits, uctx);

    nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(Au, uctx, A, ctx,
                        zinfo->perm, I->Amin_exp, I->Gstride, NULL, 0);
    nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(Bu, uctx, B, ctx,
                        zinfo->perm, I->Bmin_exp, I->Gstride, NULL, 0);

    FLINT_ASSERT(Au->bits == wbits);
    FLINT_ASSERT(Bu->bits == wbits);
    FLINT_ASSERT(Au->length > 1);
    FLINT_ASSERT(Bu->length > 1);

    success = nmod_mpolyu_content_mpoly_threaded_pool(Ac, Au, uctx, NULL, 0);
    success = success && nmod_mpolyu_content_mpoly_threaded_pool(Bc, Bu, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    nmod_mpolyu_divexact_mpoly_inplace(Au, Ac, uctx);
    nmod_mpolyu_divexact_mpoly_inplace(Bu, Bc, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */
    success = nmod_mpolyu_gcdm_zippel(Gu, Abaru, Bbaru, Au, Bu,
                                                       uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    success = _nmod_mpoly_gcd_cofactors_threaded_pool(Gc, wbits, Abarc, wbits, Bbarc, wbits,
                                                        Ac, Bc, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    nmod_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);
    nmod_mpolyu_mul_mpoly_inplace(Abaru, Abarc, uctx);
    nmod_mpolyu_mul_mpoly_inplace(Bbaru, Bbarc, uctx);

    nmod_mpoly_from_mpolyu_perm_inflate(G, I->Gbits, ctx, Gu, uctx,
                                         zinfo->perm, I->Gmin_exp, I->Gstride);
    nmod_mpoly_from_mpolyu_perm_inflate(Abar, I->Abarbits, ctx, Abaru, uctx,
                                      zinfo->perm, I->Abarmin_exp, I->Gstride);
    nmod_mpoly_from_mpolyu_perm_inflate(Bbar, I->Bbarbits, ctx, Bbaru, uctx,
                                      zinfo->perm, I->Bbarmin_exp, I->Gstride);
    success = 1;

cleanup:

    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Gu, uctx);
    nmod_mpolyu_clear(Abaru, uctx);
    nmod_mpolyu_clear(Bbaru, uctx);
    nmod_mpoly_clear(Ac, uctx);
    nmod_mpoly_clear(Bc, uctx);
    nmod_mpoly_clear(Gc, uctx);
    nmod_mpoly_clear(Abarc, uctx);
    nmod_mpoly_clear(Bbarc, uctx);

    nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    return success;
}


/*********************** Hit A and B with brown ******************************/
typedef struct
{
    nmod_mpolyn_struct * Pn;
    const nmod_mpoly_ctx_struct * nctx;
    const nmod_mpoly_struct * P;
    const nmod_mpoly_ctx_struct * ctx;
    const slong * perm;
    const ulong * shift, * stride;
    const thread_pool_handle * handles;
    slong num_handles;
}
_convertn_arg_struct;

typedef _convertn_arg_struct _convertn_arg_t[1];

static void _worker_convertn(void * varg)
{
    _convertn_arg_struct * arg = (_convertn_arg_struct *) varg;

    nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(arg->Pn, arg->nctx, arg->P, arg->ctx,
           arg->perm, arg->shift, arg->stride, arg->handles, arg->num_handles);
}

static int _try_brown(
    nmod_mpoly_t G,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    mpoly_gcd_info_t I,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int success;
    slong m = I->mvars;
    flint_bitcnt_t wbits;
    nmod_mpoly_ctx_t nctx;
    nmod_mpolyn_t An, Bn, Gn, Abarn, Bbarn;
    nmod_poly_stack_t Sp;

    if (!I->can_use_brown)
        return 0;

    FLINT_ASSERT(m >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    wbits = FLINT_MAX(A->bits, B->bits);

    nmod_mpoly_ctx_init(nctx, m, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_poly_stack_init(Sp, wbits, nctx);
    nmod_mpolyn_init(An, wbits, nctx);
    nmod_mpolyn_init(Bn, wbits, nctx);
    nmod_mpolyn_init(Gn, wbits, nctx);
    nmod_mpolyn_init(Abarn, wbits, nctx);
    nmod_mpolyn_init(Bbarn, wbits, nctx);

    if (num_handles > 0)
    {
        slong s = mpoly_divide_threads(num_handles, A->length, B->length);
        _convertn_arg_t arg;

        FLINT_ASSERT(s >= 0);
        FLINT_ASSERT(s < num_handles);

        arg->Pn = Bn;
        arg->nctx = nctx;
        arg->P = B;
        arg->ctx = ctx;
        arg->perm = I->brown_perm;
        arg->shift = I->Bmin_exp;
        arg->stride = I->Gstride;
        arg->handles = handles + (s + 1);
        arg->num_handles = num_handles - (s + 1);

        thread_pool_wake(global_thread_pool, handles[s], 0, _worker_convertn, arg);

        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(An, nctx, A, ctx,
                                   I->brown_perm, I->Amin_exp, I->Gstride,
                                                               handles + 0, s);

        thread_pool_wait(global_thread_pool, handles[s]);
    }
    else
    {
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(An, nctx, A, ctx,
                              I->brown_perm, I->Amin_exp, I->Gstride, NULL, 0);
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(Bn, nctx, B, ctx,
                              I->brown_perm, I->Bmin_exp, I->Gstride, NULL, 0);
    }

    FLINT_ASSERT(An->bits == wbits);
    FLINT_ASSERT(Bn->bits == wbits);
    FLINT_ASSERT(An->length > 1);
    FLINT_ASSERT(Bn->length > 1);

    success = (num_handles > 0)
        ? nmod_mpolyn_gcd_brown_smprime_threaded_pool(Gn, Abarn, Bbarn, An, Bn,
                                          m - 1, nctx, I, handles, num_handles)
        : nmod_mpolyn_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn,
                                                           m - 1, nctx, I, Sp);

    if (!success)
    {
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(An, nctx, A, ctx,
                              I->brown_perm, I->Amin_exp, I->Gstride, NULL, 0);
        nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(Bn, nctx, B, ctx,
                              I->brown_perm, I->Bmin_exp, I->Gstride, NULL, 0);
        success = nmod_mpolyn_gcd_brown_lgprime(Gn, Abarn, Bbarn, An, Bn,
                                                                  m - 1, nctx);
    }

    if (!success)
        goto cleanup;

    nmod_mpoly_from_mpolyn_perm_inflate(G, I->Gbits, ctx, Gn, nctx,
                                       I->brown_perm, I->Gmin_exp, I->Gstride);
    nmod_mpoly_from_mpolyn_perm_inflate(Abar, I->Abarbits, ctx, Abarn, nctx,
                                    I->brown_perm, I->Abarmin_exp, I->Gstride);
    nmod_mpoly_from_mpolyn_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarn, nctx,
                                    I->brown_perm, I->Bbarmin_exp, I->Gstride);
    success = 1;

cleanup:

    nmod_mpolyn_clear(An, nctx);
    nmod_mpolyn_clear(Bn, nctx);
    nmod_mpolyn_clear(Gn, nctx);
    nmod_mpolyn_clear(Abarn, nctx);
    nmod_mpolyn_clear(Bbarn, nctx);
    nmod_poly_stack_clear(Sp);
    nmod_mpoly_ctx_clear(nctx);

    return success;
}



/*
    The function must pack its answer into bits = Gbits <= FLINT_BITS
    Both A and B have to be packed into bits <= FLINT_BITS

    return is 1 for success, 0 for failure.
*/
int _nmod_mpoly_gcd_cofactors_threaded_pool(
    nmod_mpoly_t G, flint_bitcnt_t Gbits,
    nmod_mpoly_t Abar, flint_bitcnt_t Abarbits,
    nmod_mpoly_t Bbar, flint_bitcnt_t Bbarbits,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
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
#if WANT_ASSERT
    nmod_mpoly_t T, Asave, Bsave;

    nmod_mpoly_init(T, ctx);
    nmod_mpoly_init(Asave, ctx);
    nmod_mpoly_init(Bsave, ctx);
    nmod_mpoly_set(Asave, A, ctx);
    nmod_mpoly_set(Bsave, B, ctx);
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
            nmod_mpoly_set(Abar, B, ctx);
            nmod_mpoly_set(Bbar, A, ctx);
            nmod_mpoly_swap(Abar, Bbar, ctx);
        }
        else if (Abar == B && Bbar != A)
        {
            nmod_mpoly_set(Bbar, B, ctx);
            nmod_mpoly_set(Abar, A, ctx);
        }
        else
        {
            nmod_mpoly_set(Abar, A, ctx);
            nmod_mpoly_set(Bbar, B, ctx);
        }

        nmod_mpoly_fit_length(G, 1, ctx);
        nmod_mpoly_fit_bits(G, I->Gbits, ctx);
        G->bits = I->Gbits;
        mpoly_set_monomial_ui(G->exps, I->Gmin_exp, I->Gbits, ctx->minfo);
        G->coeffs[0] = UWORD(1);
        G->length = 1;

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
        nmod_poly_t a, b, g, t;

        nmod_poly_init_mod(a, ctx->ffinfo->mod);
        nmod_poly_init_mod(b, ctx->ffinfo->mod);
        nmod_poly_init_mod(g, ctx->ffinfo->mod);
        nmod_poly_init_mod(t, ctx->ffinfo->mod);

        _nmod_mpoly_to_nmod_poly_deflate(a, A, v_in_both,
                                                 I->Amin_exp, I->Gstride, ctx);
        _nmod_mpoly_to_nmod_poly_deflate(b, B, v_in_both,
                                                 I->Bmin_exp, I->Gstride, ctx);
        nmod_poly_gcd(g, a, b);
        _nmod_mpoly_from_nmod_poly_inflate(G, I->Gbits, g, v_in_both,
                                                 I->Gmin_exp, I->Gstride, ctx);
        nmod_poly_div(t, a, g);
        _nmod_mpoly_from_nmod_poly_inflate(Abar, I->Abarbits, t, v_in_both,
                                              I->Abarmin_exp, I->Gstride, ctx);
        nmod_poly_div(t, b, g);
        _nmod_mpoly_from_nmod_poly_inflate(Bbar, I->Bbarbits, t, v_in_both,
                                              I->Bbarmin_exp, I->Gstride, ctx);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);
        nmod_poly_clear(t);

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

    mpoly_gcd_info_set_estimates_nmod_mpoly(I, A, B, ctx, handles, num_handles);
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
                                A, try_a, B, try_b, ctx, handles, num_handles))
        {
            goto successful;
        }
    }

    mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
    mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);

    if (I->zippel_time_est < I->brown_time_est)
    {
        if (_try_zippel(G, Abar, Bbar, A, B, I, ctx))
            goto successful;

        if (_try_brown(G, Abar, Bbar, A, B, I, ctx, handles, num_handles))
            goto successful;
    }
    else
    {
        if (_try_brown(G, Abar, Bbar, A, B, I, ctx, handles, num_handles))
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

        nmod_mpoly_repack_bits_inplace(G, I->Gbits, ctx);
        nmod_mpoly_repack_bits_inplace(Abar, I->Abarbits, ctx);
        nmod_mpoly_repack_bits_inplace(Bbar, I->Bbarbits, ctx);

        if (G->coeffs[0] != 1)
        {
            _nmod_vec_scalar_mul_nmod(Abar->coeffs, Abar->coeffs,
                                 Abar->length, G->coeffs[0], ctx->ffinfo->mod);
            _nmod_vec_scalar_mul_nmod(Bbar->coeffs, Bbar->coeffs,
                                 Bbar->length, G->coeffs[0], ctx->ffinfo->mod);
            _nmod_vec_scalar_mul_nmod(G->coeffs, G->coeffs, G->length,
                   nmod_inv(G->coeffs[0], ctx->ffinfo->mod), ctx->ffinfo->mod);
        }

#if WANT_ASSERT
        nmod_mpoly_mul(T, G, Abar, ctx);
        FLINT_ASSERT(nmod_mpoly_equal(T, Asave, ctx));
        nmod_mpoly_mul(T, G, Bbar, ctx);
        FLINT_ASSERT(nmod_mpoly_equal(T, Bsave, ctx));
#endif
    }

#if WANT_ASSERT
    nmod_mpoly_clear(T, ctx);
    nmod_mpoly_clear(Asave, ctx);
    nmod_mpoly_clear(Bsave, ctx);
#endif

    mpoly_gcd_info_clear(I);

    return success;
}


int nmod_mpoly_gcd_cofactors(
    nmod_mpoly_t G,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t Gbits;
    int success;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit;
    nmod_mpoly_t Anew, Bnew;

    thread_limit = FLINT_MIN(A->length, B->length)/256;

    if (A->length == 0)
    {
        if (B->length == 0)
        {
            nmod_mpoly_zero(G, ctx);
            nmod_mpoly_zero(Abar, ctx);
            nmod_mpoly_zero(Bbar, ctx);
            return 1;
        }
        nmod_mpoly_set(G, B, ctx);
        nmod_mpoly_zero(Abar, ctx);
        nmod_mpoly_one(Bbar, ctx);
        if (G->coeffs[0] != 1)
        {
            _nmod_vec_scalar_mul_nmod(Bbar->coeffs, Bbar->coeffs,
                                 Bbar->length, G->coeffs[0], ctx->ffinfo->mod);
            _nmod_vec_scalar_mul_nmod(G->coeffs, G->coeffs, G->length,
                   nmod_inv(G->coeffs[0], ctx->ffinfo->mod), ctx->ffinfo->mod);
        }
        return 1;
    }

    if (B->length == 0)
    {
        nmod_mpoly_set(G, A, ctx);
        nmod_mpoly_zero(Bbar, ctx);
        nmod_mpoly_one(Abar, ctx);
        if (G->coeffs[0] != 1)
        {
            _nmod_vec_scalar_mul_nmod(Abar->coeffs, Abar->coeffs,
                                 Abar->length, G->coeffs[0], ctx->ffinfo->mod);
            _nmod_vec_scalar_mul_nmod(G->coeffs, G->coeffs, G->length,
                   nmod_inv(G->coeffs[0], ctx->ffinfo->mod), ctx->ffinfo->mod);
        }
        return 1;
    }

    Gbits = FLINT_MIN(A->bits, B->bits);

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
    {
        num_handles = flint_request_threads(&handles, thread_limit);
        success = _nmod_mpoly_gcd_cofactors_threaded_pool(G, Gbits,
                Abar, A->bits, Bbar, B->bits, A, B, ctx, handles, num_handles);
        flint_give_back_threads(handles, num_handles);
        return success;
    }

    nmod_mpoly_init(Anew, ctx);
    nmod_mpoly_init(Bnew, ctx);

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
        const nmod_mpoly_struct * Ause, * Buse;

        Ause = A;
        if (A->bits > FLINT_BITS)
        {
            if (!nmod_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx))
                goto could_not_repack;
            Ause = Anew;
        }

        Buse = B;
        if (B->bits > FLINT_BITS)
        {
            if (!nmod_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx))
                goto could_not_repack;
            Buse = Bnew;
        }

        num_handles = flint_request_threads(&handles, thread_limit);
        Gbits = FLINT_MIN(Ause->bits, Buse->bits);
        success = _nmod_mpoly_gcd_cofactors_threaded_pool(G, Gbits,
                              Abar, Ause->bits, Bbar, Buse->bits, Ause, Buse,
                                                    ctx, handles, num_handles);
        flint_give_back_threads(handles, num_handles);

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

        nmod_mpoly_deflation(Ashift, Astride, A, ctx);
        nmod_mpoly_deflation(Bshift, Bstride, B, ctx);
        _fmpz_vec_min(Gshift, Ashift, Bshift, ctx->minfo->nvars);
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_gcd(Gstride + k, Astride + k, Bstride + k);
        }

        success = 0;

        nmod_mpoly_deflate(Anew, A, Ashift, Gstride, ctx);
        if (Anew->bits > FLINT_BITS)
        {
            if (!nmod_mpoly_repack_bits(Anew, Anew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        nmod_mpoly_deflate(Bnew, B, Bshift, Gstride, ctx);
        if (Bnew->bits > FLINT_BITS)
        {
            if (!nmod_mpoly_repack_bits(Bnew, Bnew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        num_handles = flint_request_threads(&handles, thread_limit);
        Gbits = FLINT_MIN(Anew->bits, Bnew->bits);
        success = _nmod_mpoly_gcd_cofactors_threaded_pool(G, Gbits, Abar, Anew->bits,
                      Bbar, Bnew->bits, Anew, Bnew, ctx, handles, num_handles);
        flint_give_back_threads(handles, num_handles);

        if (!success)
            goto deflate_cleanup;

        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_sub(Ashift + k, Ashift + k, Gshift + k);
            fmpz_sub(Bshift + k, Bshift + k, Gshift + k);
            FLINT_ASSERT(fmpz_sgn(Ashift + k) >= 0);
            FLINT_ASSERT(fmpz_sgn(Bshift + k) >= 0);
        }

        nmod_mpoly_inflate(G, G, Gshift, Gstride, ctx);
        nmod_mpoly_inflate(Abar, Abar, Ashift, Gstride, ctx);
        nmod_mpoly_inflate(Bbar, Bbar, Bshift, Gstride, ctx);

        FLINT_ASSERT(G->length > 0);
        if (G->coeffs[0] != 1)
        {
            _nmod_vec_scalar_mul_nmod(Abar->coeffs, Abar->coeffs,
                                 Abar->length, G->coeffs[0], ctx->ffinfo->mod);
            _nmod_vec_scalar_mul_nmod(Bbar->coeffs, Bbar->coeffs,
                                 Bbar->length, G->coeffs[0], ctx->ffinfo->mod);
            _nmod_vec_scalar_mul_nmod(G->coeffs, G->coeffs, G->length,
                   nmod_inv(G->coeffs[0], ctx->ffinfo->mod), ctx->ffinfo->mod);
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

        nmod_mpoly_clear(Anew, ctx);
        nmod_mpoly_clear(Bnew, ctx);

        return success;
    }
}

