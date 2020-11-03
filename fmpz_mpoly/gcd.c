/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    For each j, set out[j] to the evaluation of A at x_i = alpha[i] (i != j)
    i.e. if nvars = 3
        out[0] = A(x, alpha[1], alpha[2])
        out[1] = A(alpha[0], x, alpha[2])
        out[2] = A(alpha[0], alpha[1], x)

    If ignore[j] is nonzero, then out[j] need not be calculated, probably
    because we shouldn't calculate it in dense form.
*/
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
    slong num_handles)
{
    slong i, j;
    slong nvars = ctx->minfo->nvars;
    slong total_limit, total_length;
    int use_direct_LUT;
    ulong varexp;
    ulong mask;
    slong * offsets, * shifts;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * Aexp = A->exps;
    fmpz * Acoeff = A->coeffs;
    mp_limb_t meval;
    mp_limb_t t;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        nmod_poly_zero(out + j);
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    /*
        two cases:
        (1) the Amax_exp[j] are small enough to calculate a direct LUT
        (2) use a LUT for exponents that are powers of two
    */

    total_limit = A->length/256;
    total_limit = FLINT_MAX(WORD(9999), total_limit);
    total_length = 0;
    use_direct_LUT = 1;
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        total_length += Amax_exp[j] + 1;
        if ((ulong) total_length > (ulong) total_limit)
            use_direct_LUT = 0;
    }

    if (use_direct_LUT)
    {
        slong off;
        mp_limb_t * LUT, ** LUTvalue, ** LUTvalueinv;

        /* value of powers of alpha[j] */
        LUT = (mp_limb_t *) flint_malloc(2*total_length*sizeof(mp_limb_t));

        /* pointers into LUT */
        LUTvalue    = (mp_limb_t **) flint_malloc(nvars*sizeof(mp_limb_t *));
        LUTvalueinv = (mp_limb_t **) flint_malloc(nvars*sizeof(mp_limb_t *));

        off = 0;
        for (j = 0; j < nvars; j++)
        {
            ulong k;
            mp_limb_t alphainvj = nmod_inv(alpha[j], (out + 0)->mod);

            LUTvalue[j] = LUT + off;
            LUTvalueinv[j] = LUT + total_length + off;
            LUTvalue[j][0] = 1;
            LUTvalueinv[j][0] = 1;
            for (k = 0; k < Amax_exp[j]; k++)
            {
                LUTvalue[j][k + 1] = nmod_mul(LUTvalue[j][k], alpha[j],
                                                               (out + 0)->mod);
                LUTvalueinv[j][k + 1] = nmod_mul(LUTvalueinv[j][k], alphainvj,
                                                               (out + 0)->mod);
            }

            off += Amax_exp[j] + 1;
        }
        FLINT_ASSERT(off == total_length);

        for (i = 0; i < A->length; i++)
        {
            meval = fmpz_fdiv_ui(Acoeff + i, out->mod.n);

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;
                FLINT_ASSERT(varexp <= Amax_exp[j]);
                meval = nmod_mul(meval, LUTvalue[j][varexp], (out + 0)->mod);
            }

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;

                if (ignore[j])
                    continue;

                t = nmod_mul(meval, LUTvalueinv[j][varexp], (out + j)->mod);

                FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j])
                                  || (varexp - Amin_exp[j]) % Astride[j] == 0);

                varexp = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                           (varexp - Amin_exp[j])/Astride[j];

                t = nmod_add(t, nmod_poly_get_coeff_ui(out + j, varexp),
                                                               (out + j)->mod);
                nmod_poly_set_coeff_ui(out + j, varexp, t);
            }
        }

        flint_free(LUT);
        flint_free(LUTvalue);
        flint_free(LUTvalueinv);
    }
    else
    {
        slong LUTlen;
        ulong * LUTmask;
        slong * LUToffset, * LUTvar;
        mp_limb_t * LUTvalue, * LUTvalueinv;
        mp_limb_t * vieval;
        mp_limb_t t, xpoweval, xinvpoweval;

        LUToffset   = (slong *) flint_malloc(N*FLINT_BITS*sizeof(slong));
        LUTmask     = (ulong *) flint_malloc(N*FLINT_BITS*sizeof(ulong));
        LUTvalue    = (mp_limb_t *) flint_malloc(N*FLINT_BITS*sizeof(mp_limb_t));
        LUTvar      = (slong *) flint_malloc(N*FLINT_BITS*sizeof(slong));
        LUTvalueinv = (mp_limb_t *) flint_malloc(N*FLINT_BITS*sizeof(mp_limb_t));

        vieval = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));

        LUTlen = 0;
        for (j = nvars - 1; j >= 0; j--)
        {
            flint_bitcnt_t bits = FLINT_BIT_COUNT(Amax_exp[j]);
            xpoweval = alpha[j]; /* xpoweval = alpha[j]^(2^i) */
            xinvpoweval = nmod_inv(xpoweval, (out + 0)->mod); /* alpha[j]^(-2^i) */
            for (i = 0; i < bits; i++)
            {
                LUToffset[LUTlen] = offsets[j];
                LUTmask[LUTlen] = (UWORD(1) << (shifts[j] + i));
                LUTvalue[LUTlen] = xpoweval;
                LUTvalueinv[LUTlen] = xinvpoweval;
                LUTvar[LUTlen] = j;
                LUTlen++;
                xpoweval = nmod_mul(xpoweval, xpoweval, (out + 0)->mod);
                xinvpoweval = nmod_mul(xinvpoweval, xinvpoweval, (out + 0)->mod);
            }

            vieval[j] = 1;
        }
        FLINT_ASSERT(LUTlen < N*FLINT_BITS);

        for (i = 0; i < A->length; i++)
        {
            meval = fmpz_fdiv_ui(Acoeff + i, (out + 0)->mod.n);

            for (j = 0; j < LUTlen; j++)
            {
                if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
                {
                    meval = nmod_mul(meval, LUTvalue[j], (out + 0)->mod);
                    vieval[LUTvar[j]] = nmod_mul(vieval[LUTvar[j]],
                                               LUTvalueinv[j], (out + 0)->mod);
                }
            }

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;

                FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j])
                                  || (varexp - Amin_exp[j]) % Astride[j] == 0);

                varexp = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                           (varexp - Amin_exp[j])/Astride[j];

                t = nmod_mul(meval, vieval[j], (out + j)->mod);
                t = nmod_add(t, nmod_poly_get_coeff_ui(out + j, varexp),
                                                               (out + j)->mod);
                nmod_poly_set_coeff_ui(out + j, varexp, t);
                vieval[j] = 1;
            }
        }

        flint_free(LUToffset);
        flint_free(LUTmask);
        flint_free(LUTvalue);
        flint_free(LUTvar);
        flint_free(LUTvalueinv);

        flint_free(vieval);
    }

    flint_free(offsets);
    flint_free(shifts);
}


void mpoly_gcd_info_set_estimates_fmpz_mpoly(
    mpoly_gcd_info_t I,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int try_count = 0;
    slong i, j;
    nmod_poly_t Geval;
    nmod_poly_struct * Aevals, * Bevals;
    mp_limb_t p = UWORD(1) << (FLINT_BITS - 1);
    mp_limb_t * alpha;
    flint_rand_t randstate;
    slong ignore_limit;
    int * ignore;

    flint_randinit(randstate);

    ignore = (int *) flint_malloc(ctx->minfo->nvars*sizeof(int));
    alpha = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
    Aevals = (nmod_poly_struct *) flint_malloc(
                                   ctx->minfo->nvars*sizeof(nmod_poly_struct));
    Bevals = (nmod_poly_struct *) flint_malloc(
                                   ctx->minfo->nvars*sizeof(nmod_poly_struct));

    nmod_poly_init(Geval, p);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        nmod_poly_init(Aevals + j, p);
        nmod_poly_init(Bevals + j, p);
    }

    ignore_limit = A->length/4096 + B->length/4096;
    ignore_limit = FLINT_MAX(WORD(9999), ignore_limit);
    I->Gdeflate_deg_bounds_are_nice = 1;
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        if (   I->Adeflate_deg[j] > ignore_limit
            || I->Bdeflate_deg[j] > ignore_limit)
        {
            ignore[j] = 1;
            I->Gdeflate_deg_bounds_are_nice = 0;
        }
        else
        {
            ignore[j] = 0;
        }
    }

try_again:

    if (++try_count > 10)
    {
        I->Gdeflate_deg_bounds_are_nice = 0;
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = (I->Gdeflate_deg_bound[j] + 1)/2;
        }

        goto cleanup;
    }

    p = n_nextprime(p, 1);
    nmod_init(&Geval->mod, p);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        alpha[j] = n_urandint(randstate, p - 1) + 1;
        nmod_init(&(Aevals + j)->mod, p);
        nmod_init(&(Bevals + j)->mod, p);
    }

    fmpz_mpoly_evals(Aevals, ignore, A, I->Amin_exp, I->Amax_exp, I->Gstride,
                                             alpha, ctx, handles, num_handles);
    fmpz_mpoly_evals(Bevals, ignore, B, I->Bmin_exp, I->Bmax_exp, I->Gstride,
                                             alpha, ctx, handles, num_handles);

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        if (ignore[j])
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = (I->Gdeflate_deg_bound[j] + 1)/2;
        }
        else
        {
            if (   I->Adeflate_deg[j] != nmod_poly_degree(Aevals + j)
                || I->Bdeflate_deg[j] != nmod_poly_degree(Bevals + j))
            {
                goto try_again;
            }

            nmod_poly_gcd(Geval, Aevals + j, Bevals + j);

            I->Gterm_count_est[j] = 0;
            I->Gdeflate_deg_bound[j] = nmod_poly_degree(Geval);
            for (i = I->Gdeflate_deg_bound[j]; i >= 0; i--)
            {
                I->Gterm_count_est[j] += (Geval->coeffs[i] != 0);
            }
        }
    }

cleanup:

    nmod_poly_clear(Geval);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        nmod_poly_clear(Aevals + j);
        nmod_poly_clear(Bevals + j);
    }

    flint_free(ignore);
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(Bevals);

    flint_randclear(randstate);

    return;
}


/*********************** Easy when B is a monomial ***************************/
static int _try_monomial_gcd(
    fmpz_mpoly_t G, flint_bitcnt_t Gbits,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t g;
    fmpz * minAfields, * minAdegs, * minBdegs;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length == 1);

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

    /* compute the coefficient of G */
    fmpz_init_set(g, B->coeffs + 0);
    _fmpz_vec_content_chained(g, A->coeffs, A->length);

    /* write G */
    fmpz_mpoly_fit_length(G, 1, ctx);
    fmpz_mpoly_fit_bits(G, Gbits, ctx);
    G->bits = Gbits;
    mpoly_set_monomial_ffmpz(G->exps, minBdegs, Gbits, ctx->minfo);
    fmpz_swap(G->coeffs + 0, g);
    fmpz_clear(g);
    _fmpz_mpoly_set_length(G, 1, ctx);

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

    return 1;
}


/********************** See if cofactors are monomials ***********************/
static int _try_monomial_cofactors(
    fmpz_mpoly_t G, flint_bitcnt_t Gbits,
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
    slong var,
    const fmpz_mpoly_t A, ulong Ashift,
    const fmpz_mpoly_t B, ulong Bshift,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fmpz_mpoly_t tG;
    fmpz_mpoly_univar_t Ax;

    fmpz_mpoly_init(tG, ctx);
    fmpz_mpoly_univar_init(Ax, ctx);

    fmpz_mpoly_to_univar(Ax, A, var, ctx);

    FLINT_ASSERT(Ax->length > 0);
    success = _fmpz_mpoly_gcd_threaded_pool(tG, Gbits, B, Ax->coeffs + 0, ctx, NULL, 0);
    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _fmpz_mpoly_gcd_threaded_pool(tG, Gbits, tG, Ax->coeffs + i, ctx, NULL, 0);
        if (!success)
            goto cleanup;
    }

    fmpz_mpoly_swap(G, tG, ctx);
    _mpoly_gen_shift_left(G->exps, G->bits, G->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

cleanup:

    fmpz_mpoly_clear(tG, ctx);
    fmpz_mpoly_univar_clear(Ax, ctx);

    return success;
}


/************************ See if B divides A ********************************/
static int _try_divides(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t BB,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int success = 0;
    fmpz_mpoly_t Q, B, M;

    fmpz_mpoly_init(Q, ctx);
    fmpz_mpoly_init(B, ctx);
    fmpz_mpoly_init(M, ctx);

    /* BB = M*B */
    fmpz_mpoly_term_content(M, BB, ctx);
    fmpz_mpoly_divides(B, BB, M, ctx);

    if ((num_handles > 0) ? _fmpz_mpoly_divides_heap_threaded_pool(Q, A, B,
                                                     ctx, handles, num_handles)
                          : fmpz_mpoly_divides_monagan_pearce(Q, A, B, ctx))
    {
        /* gcd(Q*B, M*B) */
        _try_monomial_gcd(G, A->bits, Q, M, ctx);
        fmpz_mpoly_mul(G, G, B, ctx);
        success = 1;
    }

    fmpz_mpoly_clear(Q, ctx);
    fmpz_mpoly_clear(B, ctx);
    fmpz_mpoly_clear(M, ctx);

    return success;
}


/********************** Hit A and B with zippel ******************************/
static int _try_zippel(
    fmpz_mpoly_t G,
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
    fmpz_mpoly_t Ac, Bc, Gc;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    if (!(I->can_use & MPOLY_GCD_USE_ZIPPEL))
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

    fmpz_mpoly_to_mpolyu_perm_deflate_threaded_pool(Au, uctx, A, ctx, zinfo->perm,
                                I->Amin_exp, I->Gstride, I->Amax_exp, NULL, 0);
    fmpz_mpoly_to_mpolyu_perm_deflate_threaded_pool(Bu, uctx, B, ctx, zinfo->perm,
                                I->Bmin_exp, I->Gstride, I->Bmax_exp, NULL, 0);

    success = fmpz_mpolyu_content_mpoly_threaded_pool(Ac, Au, uctx, NULL, 0);
    success = success &&
              fmpz_mpolyu_content_mpoly_threaded_pool(Bc, Bu, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_divexact_mpoly_inplace(Au, Ac, uctx);
    fmpz_mpolyu_divexact_mpoly_inplace(Bu, Bc, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */
    success = fmpz_mpolyu_gcdm_zippel(Gu, Abaru, Bbaru, Au, Bu,
                                                       uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    success = _fmpz_mpoly_gcd_threaded_pool(Gc, wbits, Ac, Bc, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);

    fmpz_mpoly_from_mpolyu_perm_inflate(G, I->Gbits, ctx, Gu, uctx,
                                         zinfo->perm, I->Gmin_exp, I->Gstride);
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

    fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(arg->Puu, arg->uctx, arg->P, arg->ctx,
                             arg->perm, arg->shift, arg->stride, arg->maxexps,
                                               arg->handles, arg->num_handles);

    arg->success = fmpz_mpolyu_content_mpoly_threaded_pool(arg->Pcontent,
                          arg->Puu, arg->uctx, arg->handles, arg->num_handles);
    if (arg->success)
    {
        fmpz_mpolyu_divexact_mpoly_inplace(arg->Puu, arg->Pcontent, arg->uctx);
    }
}

static int _try_bma(
    fmpz_mpoly_t G,
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
    fmpz_mpoly_t Ac, Bc, Gc, Gamma;
    slong max_minor_degree;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!(I->can_use & MPOLY_GCD_USE_ZIPPEL2))
        return 0;

    FLINT_ASSERT(m >= WORD(3));

    /* uctx is context for Z[y_2,...,y_{m - 1}]*/
    fmpz_mpoly_ctx_init(uctx, m - 2, ORD_LEX);

    max_minor_degree = 0;
    for (i = 2; i < m; i++)
    {
        k = I->zippel2_perm[i];
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
        arg->perm = I->zippel2_perm;
        arg->shift = I->Bmin_exp;
        arg->stride = I->Gstride;
        arg->maxexps = I->Bmax_exp;
        arg->handles = handles + (s + 1);
        arg->num_handles = num_handles - (s + 1);

        thread_pool_wake(global_thread_pool, handles[s], 0, _worker_convertuu, arg);

        fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Auu, uctx, A, ctx,
                        I->zippel2_perm, I->Amin_exp, I->Gstride, I->Amax_exp,
                                                               handles + 0, s);
        success = fmpz_mpolyu_content_mpoly_threaded_pool(Ac, Auu,
                                                         uctx, handles + 0, s);
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
        fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Auu, uctx, A, ctx,
               I->zippel2_perm, I->Amin_exp, I->Gstride, I->Amax_exp, NULL, 0);
        fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Buu, uctx, B, ctx,
               I->zippel2_perm, I->Bmin_exp, I->Gstride, I->Bmax_exp, NULL, 0);

        success = fmpz_mpolyu_content_mpoly_threaded_pool(Ac, Auu, uctx, NULL, 0);
        success = success &&
                  fmpz_mpolyu_content_mpoly_threaded_pool(Bc, Buu, uctx, NULL, 0);
        if (!success)
            goto cleanup;

        fmpz_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);
        fmpz_mpolyu_divexact_mpoly_inplace(Buu, Bc, uctx);
    }

    FLINT_ASSERT(Auu->bits == wbits);
    FLINT_ASSERT(Buu->bits == wbits);
    FLINT_ASSERT(Auu->length > 1);
    FLINT_ASSERT(Buu->length > 1);
    FLINT_ASSERT(Ac->bits == wbits);
    FLINT_ASSERT(Bc->bits == wbits);

    if (fmpz_mpolyu_equal_upto_unit(Auu, Buu, uctx))
    {
        fmpz_mpolyu_swap(Guu, Auu, uctx);
        goto done;
    }

    success = _fmpz_mpoly_gcd_threaded_pool(Gamma, wbits, Auu->coeffs + 0,
                                 Buu->coeffs + 0, uctx, handles, num_handles);
    if (!success)
        goto cleanup;

    success = (num_handles > 0)
           ? fmpz_mpolyuu_gcd_berlekamp_massey_threaded_pool(Guu, Abaruu, Bbaruu,
                                  Auu, Buu, Gamma, uctx, handles, num_handles)
           : fmpz_mpolyuu_gcd_berlekamp_massey(Guu, Abaruu, Bbaruu,
                                                        Auu, Buu, Gamma, uctx);
    if (!success)
        goto cleanup;
done:

    success = _fmpz_mpoly_gcd_threaded_pool(Gc, wbits, Ac, Bc, uctx, handles, num_handles);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_mul_mpoly_inplace(Guu, Gc, uctx);

    fmpz_mpoly_from_mpolyuu_perm_inflate(G, I->Gbits, ctx, Guu, uctx,
                                     I->zippel2_perm, I->Gmin_exp, I->Gstride);
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

    fmpz_mpoly_to_mpoly_perm_deflate_threaded_pool(arg->Pl, arg->lctx, arg->P, arg->ctx,
                                         arg->perm, arg->shift, arg->stride,
                                               arg->handles, arg->num_handles);
}

static int _try_brown(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int success;
    slong m = I->mvars;
    flint_bitcnt_t wbits;
    fmpz_mpoly_ctx_t lctx;
    fmpz_mpoly_t Al, Bl, Gl, Abarl, Bbarl;

    if (!(I->can_use & MPOLY_GCD_USE_BROWN))
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

    /* convert to univariate format */
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

        thread_pool_wake(global_thread_pool, handles[s], 0, _worker_convertu, arg);

        fmpz_mpoly_to_mpoly_perm_deflate_threaded_pool(Al, lctx, A, ctx,
                                    I->brown_perm, I->Amin_exp, I->Gstride,
                                                               handles + 0, s);

        thread_pool_wait(global_thread_pool, handles[s]);
    }
    else
    {
        fmpz_mpoly_to_mpoly_perm_deflate_threaded_pool(Al, lctx, A, ctx,
                              I->brown_perm, I->Amin_exp, I->Gstride, NULL, 0);
        fmpz_mpoly_to_mpoly_perm_deflate_threaded_pool(Bl, lctx, B, ctx,
                              I->brown_perm, I->Bmin_exp, I->Gstride, NULL, 0);
    }

    FLINT_ASSERT(Al->bits == wbits);
    FLINT_ASSERT(Bl->bits == wbits);
    FLINT_ASSERT(Al->length > 1);
    FLINT_ASSERT(Bl->length > 1);

    success = (num_handles > 0) ?
                fmpz_mpolyl_gcd_brown_threaded_pool(Gl, Abarl, Bbarl, Al, Bl,
                                               lctx, I, handles, num_handles) :
                fmpz_mpolyl_gcd_brown(Gl, Abarl, Bbarl, Al, Bl, lctx, I);

    if (!success)
        goto cleanup;

    fmpz_mpoly_from_mpoly_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                       I->brown_perm, I->Gmin_exp, I->Gstride);
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
    The function must pack a successful answer into bits = Gbits <= FLINT_BITS.
    Both A and B have to be packed into bits <= FLINT_BITS.
    return is 1 for success, 0 for failure.
*/
int _fmpz_mpoly_gcd_threaded_pool(
    fmpz_mpoly_t G, flint_bitcnt_t Gbits,
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

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Gbits <= FLINT_BITS);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    if (A->length == 1)
    {
        return _try_monomial_gcd(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _try_monomial_gcd(G, Gbits, A, B, ctx);
    }

    mpoly_gcd_info_init(I, nvars);

    /* entries of I are all now invalid */

    I->Gbits = Gbits;

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
    if (_try_monomial_cofactors(G, I->Gbits, A, B, ctx))
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
        I->Gmin_exp[j] = FLINT_MIN(I->Amin_exp[j], I->Bmin_exp[j]);
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
        fmpz_t gA, gB;

calculate_trivial_gcd:

        fmpz_init(gA);
        fmpz_init(gB);
        _fmpz_vec_content(gA, A->coeffs, A->length);
        _fmpz_vec_content(gB, B->coeffs, B->length);

        fmpz_mpoly_fit_length(G, 1, ctx);
        fmpz_mpoly_fit_bits(G, Gbits, ctx);
        G->bits = Gbits;
        mpoly_set_monomial_ui(G->exps, I->Gmin_exp, Gbits, ctx->minfo);
        fmpz_gcd(G->coeffs + 0, gA, gB);
        _fmpz_mpoly_set_length(G, 1, ctx);

        fmpz_clear(gA);
        fmpz_clear(gB);

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
        fmpz_poly_t a, b, g;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        _fmpz_mpoly_to_fmpz_poly_deflate(a, A, v_in_both,
                                                 I->Amin_exp, I->Gstride, ctx);
        _fmpz_mpoly_to_fmpz_poly_deflate(b, B, v_in_both,
                                                 I->Bmin_exp, I->Gstride, ctx);
        fmpz_poly_gcd(g, a, b);
        _fmpz_mpoly_from_fmpz_poly_inflate(G, Gbits, g, v_in_both,
                                                 I->Gmin_exp, I->Gstride, ctx);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);

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
                                   v_in_A_only,
                                   A, I->Amin_exp[v_in_A_only],
                                   B, I->Bmin_exp[v_in_A_only],
                                   ctx);
        goto cleanup;
    }
    if (v_in_B_only != -WORD(1))
    {
        success = _try_missing_var(G, I->Gbits,
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

    /* check divisibility A/B and B/A */
    {
        int gcd_is_trivial = 1;
        int try_a = I->Gdeflate_deg_bounds_are_nice;
        int try_b = I->Gdeflate_deg_bounds_are_nice;

        for (j = 0; j < nvars; j++)
        {
            if (I->Gdeflate_deg_bound[j] != 0)
                gcd_is_trivial = 0;

            if (I->Adeflate_deg[j] != I->Gdeflate_deg_bound[j])
                try_a = 0;

            if (I->Bdeflate_deg[j] != I->Gdeflate_deg_bound[j])
                try_b = 0;
        }

        if (gcd_is_trivial)
            goto calculate_trivial_gcd;

        if (try_a && _try_divides(G, B, A, ctx, handles, num_handles))
            goto successful;

        if (try_b && _try_divides(G, A, B, ctx, handles, num_handles))
            goto successful;
    }

    mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
    mpoly_gcd_info_measure_bma(I, A->length, B->length, ctx->minfo);

    if (I->mvars < 3)
    {
        /* TODO: bivariate heuristic here */

        if (_try_brown(G, A, B, I, ctx, handles, num_handles))
            goto successful;
    }
    else if ((I->can_use & MPOLY_GCD_USE_BROWN) &&
             (I->can_use & MPOLY_GCD_USE_ZIPPEL2) &&
             I->zippel2_time < I->brown_time &&
             (I->mvars*(I->Adensity + I->Bdensity) < 1 ||
              I->zippel2_time < 0.01*I->brown_time))
    {
        if (_try_bma(G, A, B, I, ctx, handles, num_handles))
            goto successful;

        if (_try_brown(G, A, B, I, ctx, handles, num_handles))
            goto successful;
    }
    else
    {
        if (_try_brown(G, A, B, I, ctx, handles, num_handles))
            goto successful;

        if (_try_bma(G, A, B, I, ctx, handles, num_handles))
            goto successful;
    }

    mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);

    if (_try_zippel(G, A, B, I, ctx))
        goto successful;

    success = 0;
    goto cleanup;

successful:

    success = 1;

cleanup:

    mpoly_gcd_info_clear(I);

    if (success)
    {
        fmpz_mpoly_repack_bits_inplace(G, Gbits, ctx);
        FLINT_ASSERT(G->length > 0);
        if (fmpz_sgn(G->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, G, ctx);
    }

    return success;
}


int fmpz_mpoly_gcd(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t Gbits;
    int success;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit;

    thread_limit = FLINT_MIN(A->length, B->length)/256;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (B->length == 0)
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }
        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_fmpz(A, ctx))
    {
        fmpz_t t;
        fmpz_init_set(t, A->coeffs);
        _fmpz_vec_content_chained(t, B->coeffs, B->length);
        fmpz_mpoly_set_fmpz(G, t, ctx);
        fmpz_clear(t);
        return 1;
    }

    if (fmpz_mpoly_is_fmpz(B, ctx))
    {
        fmpz_t t;
        fmpz_init_set(t, B->coeffs);
        _fmpz_vec_content_chained(t, A->coeffs, A->length);
        fmpz_mpoly_set_fmpz(G, t, ctx);
        fmpz_clear(t);
        return 1;
    }

    Gbits = FLINT_MIN(A->bits, B->bits);

    if (A->length == 1)
    {
        return _try_monomial_gcd(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _try_monomial_gcd(G, Gbits, A, B, ctx);
    }

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
    {
        num_handles = flint_request_threads(&handles, thread_limit);
        success = _fmpz_mpoly_gcd_threaded_pool(G, Gbits, A, B, ctx, handles, num_handles);
        flint_give_back_threads(handles, num_handles);
        return success;
    }

    if (_try_monomial_cofactors(G, Gbits, A, B, ctx))
    {
        return 1;
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
        fmpz_mpoly_t Anew, Bnew;
        const fmpz_mpoly_struct * Ause, * Buse;

        fmpz_mpoly_init(Anew, ctx);
        fmpz_mpoly_init(Bnew, ctx);

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

        num_handles = flint_request_threads(&handles, thread_limit);
        Gbits = FLINT_MIN(Ause->bits, Buse->bits);
        success = _fmpz_mpoly_gcd_threaded_pool(G, Gbits, Ause, Buse, ctx,
                                                         handles, num_handles);
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

        num_handles = flint_request_threads(&handles, thread_limit);
        Gbits = FLINT_MIN(Anew->bits, Bnew->bits);
        success = _fmpz_mpoly_gcd_threaded_pool(G, Gbits, Anew, Bnew, ctx,
                                                         handles, num_handles);
        flint_give_back_threads(handles, num_handles);

        if (success)
        {
            fmpz_mpoly_inflate(G, G, Gshift, Gstride, ctx);

            /* inflation may have changed the lc */
            FLINT_ASSERT(G->length > 0);
            if (fmpz_sgn(G->coeffs + 0) < 0)
            {
                fmpz_mpoly_neg(G, G, ctx);
            }
        }

deflate_cleanup:

        _fmpz_vec_clear(Ashift, ctx->minfo->nvars);
        _fmpz_vec_clear(Astride, ctx->minfo->nvars);
        _fmpz_vec_clear(Bshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Bstride, ctx->minfo->nvars);
        _fmpz_vec_clear(Gshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Gstride, ctx->minfo->nvars);

cleanup:

        fmpz_mpoly_clear(Anew, ctx);
        fmpz_mpoly_clear(Bnew, ctx);

        return success;
    }
}

