/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/*
    For each j, set out[j] to the evaluation of A at x_i = alpha[i] (i != j)
    i.e. if nvars = 3
        out[0] = A(x, alpha[1], alpha[2])
        out[1] = A(alpha[0], x, alpha[2])
        out[2] = A(alpha[0], alpha[1], x)

    If ignore[j] is nonzero, then out[j] need not be calculated, probably
    because we shouldn't calculate it in dense form.
*/
void fq_nmod_mpoly_evals(
    fq_nmod_poly_struct * out,
    const int * ignore,
    const fq_nmod_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx)
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
    fq_nmod_struct * Acoeff = A->coeffs;
    fq_nmod_t meval, t, t2;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    fq_nmod_init(meval, ctx->fqctx);
    fq_nmod_init(t, ctx->fqctx);
    fq_nmod_init(t2, ctx->fqctx);

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        fq_nmod_poly_zero(out + j, ctx->fqctx);
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
        fq_nmod_struct * LUT, ** LUTvalue, ** LUTvalueinv;

        /* value of powers of alpha[j] */
        LUT = (fq_nmod_struct *) flint_malloc(2*total_length*sizeof(fq_nmod_struct));
        for (j = 0; j < 2*total_length; j++)
            fq_nmod_init(LUT + j, ctx->fqctx);

        /* pointers into LUT */
        LUTvalue    = (fq_nmod_struct **) flint_malloc(nvars*sizeof(fq_nmod_struct *));
        LUTvalueinv = (fq_nmod_struct **) flint_malloc(nvars*sizeof(fq_nmod_struct *));

        off = 0;
        for (j = 0; j < nvars; j++)
        {
            ulong k;

            fq_nmod_inv(t2, alpha + j, ctx->fqctx);

            LUTvalue[j] = LUT + off;
            LUTvalueinv[j] = LUT + total_length + off;
            fq_nmod_one(LUTvalue[j] + 0, ctx->fqctx);
            fq_nmod_one(LUTvalueinv[j] + 0, ctx->fqctx);
            for (k = 0; k < Amax_exp[j]; k++)
            {
                fq_nmod_mul(LUTvalue[j] + k + 1, LUTvalue[j] + k,
                                                        alpha + j, ctx->fqctx);
                fq_nmod_mul(LUTvalueinv[j] + k + 1, LUTvalueinv[j] + k,
                                                               t2, ctx->fqctx);
            }

            off += Amax_exp[j] + 1;
        }
        FLINT_ASSERT(off == total_length);

        for (i = 0; i < A->length; i++)
        {
            fq_nmod_set(meval, Acoeff + i, ctx->fqctx);

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;
                FLINT_ASSERT(varexp <= Amax_exp[j]);
                fq_nmod_mul(t, meval, LUTvalue[j] + varexp, ctx->fqctx);
                fq_nmod_swap(meval, t, ctx->fqctx);
            }

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;

                if (ignore[j])
                    continue;

                fq_nmod_mul(t, meval, LUTvalueinv[j] + varexp, ctx->fqctx);

                FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j])
                                  || (varexp - Amin_exp[j]) % Astride[j] == 0);

                varexp = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                           (varexp - Amin_exp[j])/Astride[j];

                fq_nmod_poly_get_coeff(t2, out + j, varexp, ctx->fqctx);
                fq_nmod_add(t, t, t2, ctx->fqctx);
                fq_nmod_poly_set_coeff(out + j, varexp, t, ctx->fqctx);
            }
        }

        for (j = 0; j < 2*total_length; j++)
            fq_nmod_clear(LUT + j, ctx->fqctx);

        flint_free(LUT);
        flint_free(LUTvalue);
        flint_free(LUTvalueinv);
    }
    else
    {
        slong LUTlen;
        ulong * LUTmask;
        slong * LUToffset, * LUTvar;
        fq_nmod_struct * LUTvalue, * LUTvalueinv;
        fq_nmod_struct * vieval;
        fq_nmod_t xpoweval, xinvpoweval;

        fq_nmod_init(xpoweval, ctx->fqctx);
        fq_nmod_init(xinvpoweval, ctx->fqctx);

        LUToffset   = (slong *) flint_malloc(N*FLINT_BITS*sizeof(slong));
        LUTmask     = (ulong *) flint_malloc(N*FLINT_BITS*sizeof(ulong));
        LUTvalue    = (fq_nmod_struct *) flint_malloc(N*FLINT_BITS*sizeof(fq_nmod_struct));
        LUTvar      = (slong *) flint_malloc(N*FLINT_BITS*sizeof(slong));
        LUTvalueinv = (fq_nmod_struct *) flint_malloc(N*FLINT_BITS*sizeof(fq_nmod_struct));
        for (j = 0; j < N*FLINT_BITS; j++)
        {
            fq_nmod_init(LUTvalue + j, ctx->fqctx);
            fq_nmod_init(LUTvalueinv + j, ctx->fqctx);
        }

        vieval = (fq_nmod_struct *) flint_malloc(nvars*sizeof(fq_nmod_struct));
        for (j = 0; j < nvars; j++)
        {
            fq_nmod_init(vieval + j, ctx->fqctx);
        }

        LUTlen = 0;
        for (j = nvars - 1; j >= 0; j--)
        {
            flint_bitcnt_t bits = FLINT_BIT_COUNT(Amax_exp[j]);
            fq_nmod_set(xpoweval, alpha + j, ctx->fqctx); /* xpoweval = alpha[j]^(2^i) */
            fq_nmod_inv(xinvpoweval, xpoweval, ctx->fqctx); /* alpha[j]^(-2^i) */
            for (i = 0; i < bits; i++)
            {
                LUToffset[LUTlen] = offsets[j];
                LUTmask[LUTlen] = (UWORD(1) << (shifts[j] + i));
                fq_nmod_set(LUTvalue + LUTlen, xpoweval, ctx->fqctx);
                fq_nmod_set(LUTvalueinv + LUTlen, xinvpoweval, ctx->fqctx);
                LUTvar[LUTlen] = j;
                LUTlen++;
                fq_nmod_mul(xpoweval, xpoweval, xpoweval, ctx->fqctx);
                fq_nmod_mul(xinvpoweval, xinvpoweval, xinvpoweval, ctx->fqctx);
            }

            fq_nmod_one(vieval + j, ctx->fqctx);
        }
        FLINT_ASSERT(LUTlen < N*FLINT_BITS);

        for (i = 0; i < A->length; i++)
        {
            fq_nmod_set(meval, Acoeff + i, ctx->fqctx);

            for (j = 0; j < LUTlen; j++)
            {
                if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
                {
                    fq_nmod_mul(meval, meval, LUTvalue + j, ctx->fqctx);
                    fq_nmod_mul(vieval + LUTvar[j], vieval + LUTvar[j],
                                                  LUTvalueinv + j, ctx->fqctx);
                }
            }

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;

                FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j])
                                  || (varexp - Amin_exp[j]) % Astride[j] == 0);

                varexp = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                           (varexp - Amin_exp[j])/Astride[j];

                fq_nmod_mul(t, meval, vieval + j, ctx->fqctx);
                fq_nmod_poly_get_coeff(t2, out + j, varexp, ctx->fqctx);
                fq_nmod_add(t, t, t2, ctx->fqctx);
                fq_nmod_poly_set_coeff(out + j, varexp, t, ctx->fqctx);
                fq_nmod_one(vieval + j, ctx->fqctx);
            }
        }

        for (j = 0; j < N*FLINT_BITS; j++)
        {
            fq_nmod_clear(LUTvalue + j, ctx->fqctx);
            fq_nmod_clear(LUTvalueinv + j, ctx->fqctx);
        }
        flint_free(LUToffset);
        flint_free(LUTmask);
        flint_free(LUTvalue);
        flint_free(LUTvar);
        flint_free(LUTvalueinv);

        for (j = 0; j < nvars; j++)
        {
            fq_nmod_clear(vieval + j, ctx->fqctx);
        }
        flint_free(vieval);

        fq_nmod_clear(xpoweval, ctx->fqctx);
        fq_nmod_clear(xinvpoweval, ctx->fqctx);
    }

    flint_free(offsets);
    flint_free(shifts);

    fq_nmod_clear(meval, ctx->fqctx);
    fq_nmod_clear(t, ctx->fqctx);
    fq_nmod_clear(t2, ctx->fqctx);
}


void mpoly_gcd_info_set_estimates_fq_nmod_mpoly(
    mpoly_gcd_info_t I,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int try_count = 0;
    slong i, j;
    fq_nmod_poly_t Geval;
    fq_nmod_poly_struct * Aevals, * Bevals;
    fq_nmod_struct * alpha;
    flint_rand_t randstate;
    slong ignore_limit;
    int * ignore;

    flint_randinit(randstate);

    ignore = (int *) flint_malloc(ctx->minfo->nvars*sizeof(int));
    alpha = (fq_nmod_struct *) flint_malloc(
                                     ctx->minfo->nvars*sizeof(fq_nmod_struct));
    Aevals = (fq_nmod_poly_struct *) flint_malloc(
                                ctx->minfo->nvars*sizeof(fq_nmod_poly_struct));
    Bevals = (fq_nmod_poly_struct *) flint_malloc(
                                ctx->minfo->nvars*sizeof(fq_nmod_poly_struct));

    fq_nmod_poly_init(Geval, ctx->fqctx);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        fq_nmod_init(alpha + j, ctx->fqctx);
        fq_nmod_poly_init(Aevals + j, ctx->fqctx);
        fq_nmod_poly_init(Bevals + j, ctx->fqctx);
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
            I->Gterm_count_est[j] = 1 + I->Gdeflate_deg_bound[j]/2;
        }

        goto cleanup;
    }

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        fq_nmod_randtest_not_zero(alpha + j, randstate, ctx->fqctx);
    }


    fq_nmod_mpoly_evals(Aevals, ignore, A, I->Amin_exp, I->Amax_exp,
                                                       I->Gstride, alpha, ctx);
    fq_nmod_mpoly_evals(Bevals, ignore, B, I->Bmin_exp, I->Bmax_exp,
                                                       I->Gstride, alpha, ctx);

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        if (ignore[j])
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = 1 + I->Gdeflate_deg_bound[j]/2;
        }
        else
        {
            if (   I->Adeflate_deg[j] != fq_nmod_poly_degree(Aevals + j, ctx->fqctx)
                || I->Bdeflate_deg[j] != fq_nmod_poly_degree(Bevals + j, ctx->fqctx))
            {
                goto try_again;
            }

            fq_nmod_poly_gcd(Geval, Aevals + j, Bevals + j, ctx->fqctx);

            I->Gterm_count_est[j] = 0;
            I->Gdeflate_deg_bound[j] = fq_nmod_poly_degree(Geval, ctx->fqctx);
            for (i = I->Gdeflate_deg_bound[j]; i >= 0; i--)
            {
                I->Gterm_count_est[j] += !fq_nmod_is_zero(Geval->coeffs + i, ctx->fqctx);
            }
        }
    }

cleanup:

    fq_nmod_poly_clear(Geval, ctx->fqctx);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        fq_nmod_clear(alpha + j, ctx->fqctx);
        fq_nmod_poly_clear(Aevals + j, ctx->fqctx);
        fq_nmod_poly_clear(Bevals + j, ctx->fqctx);
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
    fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
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

    fq_nmod_mpoly_fit_length(G, 1, ctx);
    fq_nmod_mpoly_fit_bits(G, Gbits, ctx);
    G->bits = Gbits;
    mpoly_set_monomial_ffmpz(G->exps, minBdegs, Gbits, ctx->minfo);
    fq_nmod_one(G->coeffs + 0, ctx->fqctx);
    G->length = 1;

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
    fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong NA, NG;
    slong nvars = ctx->minfo->nvars;
    fmpz * Abarexps, * Bbarexps, * Texps;
    fq_nmod_t t1, t2;
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

    fq_nmod_clear(t1, ctx->fqctx);
    fq_nmod_clear(t2, ctx->fqctx);

    return success;
}


/********* Assume B has length one when converted to univar format ***********/
static int _try_missing_var(
    fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
    slong var,
    const fq_nmod_mpoly_t A, ulong Ashift,
    const fq_nmod_mpoly_t B, ulong Bshift,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fq_nmod_mpoly_t tG;
    fq_nmod_mpoly_univar_t Ax;

    fq_nmod_mpoly_init(tG, ctx);
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

    fq_nmod_mpoly_swap(G, tG, ctx);
    _mpoly_gen_shift_left(G->exps, G->bits, G->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

cleanup:

    fq_nmod_mpoly_clear(tG, ctx);
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
        success = 1;
        goto cleanup;
    }

    if (try_a && fq_nmod_mpoly_divides(Q, B, A, ctx))
    {
        fq_nmod_mpoly_set(G, A, ctx);
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
    fq_nmod_mpoly_t Ac, Bc, Gc;

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

    success = _fq_nmod_mpoly_gcd(Gc, wbits, Ac, Bc, uctx);
    if (!success)
        goto cleanup;

    fq_nmod_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);

    fq_nmod_mpoly_from_mpolyu_perm_inflate(G, I->Gbits, ctx, Gu, uctx,
                                         zinfo->perm, I->Gmin_exp, I->Gstride);
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

    fq_nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    return success;
}


/*********************** Hit A and B with brown ******************************/
static int _try_brown(
    fq_nmod_mpoly_t G,
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


/*
    The function must pack its answer into bits = Gbits <= FLINT_BITS
    Both A and B have to be packed into bits <= FLINT_BITS

    return is 1 for success, 0 for failure.
*/
int _fq_nmod_mpoly_gcd(
    fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
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

calculate_trivial_gcd:

        fq_nmod_mpoly_fit_length(G, 1, ctx);
        fq_nmod_mpoly_fit_bits(G, Gbits, ctx);
        G->bits = Gbits;
        mpoly_set_monomial_ui(G->exps, I->Gmin_exp, Gbits, ctx->minfo);
        fq_nmod_one(G->coeffs + 0, ctx->fqctx);
        _fq_nmod_mpoly_set_length(G, 1, ctx);

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
        fq_nmod_poly_t a, b, g;

        fq_nmod_poly_init(a, ctx->fqctx);
        fq_nmod_poly_init(b, ctx->fqctx);
        fq_nmod_poly_init(g, ctx->fqctx);
        _fq_nmod_mpoly_to_fq_nmod_poly_deflate(a, A, v_in_both,
                                                 I->Amin_exp, I->Gstride, ctx);
        _fq_nmod_mpoly_to_fq_nmod_poly_deflate(b, B, v_in_both,
                                                 I->Bmin_exp, I->Gstride, ctx);
        fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
        _fq_nmod_mpoly_from_fq_nmod_poly_inflate(G, Gbits, g, v_in_both,
                                                 I->Gmin_exp, I->Gstride, ctx);
        fq_nmod_poly_clear(a, ctx->fqctx);
        fq_nmod_poly_clear(b, ctx->fqctx);
        fq_nmod_poly_clear(g, ctx->fqctx);

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

    /* all variable are now either
            missing from both ess(A) and ess(B), or
            present in both ess(A) and ess(B)
        and there are at least two in the latter case
    */

    mpoly_gcd_info_set_estimates_fq_nmod_mpoly(I, A, B, ctx);
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

        if ((try_a || try_b) && _try_divides(G, A, try_a, B, try_b, ctx))
            goto successful;
    }

    mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
    mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);

    if (I->zippel_time_est < I->brown_time_est)
    {
        if (_try_zippel(G, A, B, I, ctx))
            goto successful;

        if (_try_brown(G, A, B, I, ctx))
            goto successful;
    }
    else
    {
        if (_try_brown(G, A, B, I, ctx))
            goto successful;

        if (_try_zippel(G, A, B, I, ctx))
            goto successful;
    }

    success = 0;
    goto cleanup;

successful:

    success = 1;

cleanup:

    mpoly_gcd_info_clear(I);

    if (success)
    {
        fq_nmod_mpoly_repack_bits_inplace(G, Gbits, ctx);
        fq_nmod_mpoly_make_monic(G, G, ctx);
    }

    return success;
}

int fq_nmod_mpoly_gcd(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t Gbits;

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        if (fq_nmod_mpoly_is_zero(B, ctx))
            fq_nmod_mpoly_zero(G, ctx);
        else
            fq_nmod_mpoly_make_monic(G, B, ctx);
        return 1;
    }

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    Gbits = FLINT_MIN(A->bits, B->bits);

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
    {
        /* usual gcd's go right down here */
        return _fq_nmod_mpoly_gcd(G, Gbits, A, B, ctx);
    }

    if (A->length == 1)
    {
        return _try_monomial_gcd(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _try_monomial_gcd(G, Gbits, A, B, ctx);
    }
    else if (_try_monomial_cofactors(G, Gbits, A, B, ctx))
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

        int success;
        int useAnew = 0;
        int useBnew = 0;
        slong k;
        fmpz * Ashift, * Astride;
        fmpz * Bshift, * Bstride;
        fmpz * Gshift, * Gstride;
        fq_nmod_mpoly_t Anew;
        fq_nmod_mpoly_t Bnew;

        fq_nmod_mpoly_init(Anew, ctx);
        fq_nmod_mpoly_init(Bnew, ctx);

        if (A->bits > FLINT_BITS)
        {
            useAnew = fq_nmod_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx);
            if (!useAnew)
                goto could_not_repack;
        }

        if (B->bits > FLINT_BITS)
        {
            useBnew = fq_nmod_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx);
            if (!useBnew)
                goto could_not_repack;
        }

        success = _fq_nmod_mpoly_gcd(G, FLINT_BITS, useAnew ? Anew : A,
                                                    useBnew ? Bnew : B, ctx);
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

        success = _fq_nmod_mpoly_gcd(G, FLINT_BITS, Anew, Bnew, ctx);

        if (success)
        {
            fq_nmod_mpoly_inflate(G, G, Gshift, Gstride, ctx);
            fq_nmod_mpoly_make_monic(G, G, ctx);
        }

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
