/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "long_extras.h"

/*
    Convert D to A if the degrees of A are <= expected_deg
    If not, return 0 and set A to 0.
*/
static int _from_dense(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    slong * Adeg_bounds,
    slong * expect_deg,
    fmpz_mod_poly_t D,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int ret;
    slong off, j, k, N;
    flint_bitcnt_t bits;
    slong nvars = ctx->minfo->nvars;
    slong Alen;
    ulong topmask, outrange;
    ulong * exps, * pcurexp, * pexps, * rangemask;
    TMP_INIT;

    FLINT_ASSERT(nvars <= FLINT_BITS);

    TMP_START;
    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    /* find bits needed for the result */
    off = 1;
    for (j = 0; j < nvars; j++)
    {
        FLINT_ASSERT(expect_deg[j] >= 0);
        off *= Adeg_bounds[j];
        exps[j] = expect_deg[j];
    }

    bits = mpoly_exp_bits_required_ui(exps, ctx->minfo);
    bits = FLINT_MAX(Abits, bits);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    N = mpoly_words_per_exp(bits, ctx->minfo);

    /* we are going to push back terms manually */
    fmpz_mod_mpoly_fit_length_reset_bits(A, 0, bits, ctx);
    Alen = 0;

    /* find exponent vector for all variables */
    pexps = (ulong *) TMP_ALLOC(N*nvars*sizeof(ulong));
    for (k = 0; k < nvars; k++)
        mpoly_gen_monomial_sp(pexps + k*N, k, bits, ctx->minfo);

    /* get most significant exponent in exps and its vector in ptempexp */
    off--;
    pcurexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    rangemask = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    outrange = 0;
    mpoly_monomial_zero(pcurexp, N);
    k = off;
    for (j = nvars - 1; j >= 0; j--) 
    {
        exps[j] = k % Adeg_bounds[j];
        rangemask[j] = UWORD(1) << j;
        outrange ^= (outrange ^ FLINT_SIGN_EXT(expect_deg[j] - exps[j])) & rangemask[j];
        k = k / Adeg_bounds[j];
        mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + N*j, N);
    }

    /* scan down through the exponents */
    topmask = 0;
    for (; off >= 0; off--)
    {
        if (off < D->length && !fmpz_is_zero(D->coeffs + off))
        {
            if (outrange)
            {
                _fmpz_mod_mpoly_set_length(A, 0, ctx);
                ret = 0;
                goto cleanup;
            }

            _fmpz_mod_mpoly_fit_length(&A->coeffs, &A->coeffs_alloc,
                                       &A->exps, &A->exps_alloc, N, Alen + 1);
            fmpz_swap(A->coeffs + Alen, D->coeffs + off);
            mpoly_monomial_set(A->exps + N*Alen, pcurexp, N);
            topmask |= (A->exps + N*Alen)[N - 1];
            Alen++;
        }

        j = nvars - 1;
        do {
            --exps[j];
            outrange ^= (outrange ^ FLINT_SIGN_EXT(expect_deg[j] - exps[j])) & rangemask[j];
            if (FLINT_SIGN_EXT(exps[j]) != 0)
            {
                FLINT_ASSERT(off == 0 || j > 0);
                FLINT_ASSERT(exps[j] == -UWORD(1));
                exps[j] = Adeg_bounds[j] - 1;
                outrange ^= (outrange ^ FLINT_SIGN_EXT(expect_deg[j] - exps[j])) & rangemask[j];
                mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + N*j, N);
            }
            else
            {
                mpoly_monomial_sub_mp(pcurexp, pcurexp, pexps + N*j, N);
                break;
            }
        } while (--j >= 0);
    }
    _fmpz_mod_mpoly_set_length(A, Alen, ctx);

    /* sort the exponents if needed */
    if (ctx->minfo->ord != ORD_LEX)
    {
        flint_bitcnt_t pos;
        mpoly_get_cmpmask(pcurexp, N, bits, ctx->minfo);
        pos = FLINT_BIT_COUNT(topmask);
        if (N == 1)
            _fmpz_mod_mpoly_radix_sort1(A->coeffs, A->exps, 0, A->length,
                                                     pos, pcurexp[0], topmask);
        else
            _fmpz_mod_mpoly_radix_sort(A->coeffs, A->exps, 0, A->length,
                                         (N - 1)*FLINT_BITS + pos, N, pcurexp);
    }

    ret = 1;

cleanup:

    TMP_END;
    return ret;
}

int _fmpz_mod_mpoly_divides_dense_maxfields(
    fmpz_mod_mpoly_t Q,
    const fmpz_mod_mpoly_t A, fmpz * maxAfields,
    const fmpz_mod_mpoly_t B, fmpz * maxBfields,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    slong nvars = ctx->minfo->nvars;
    fmpz_mod_poly_t Ad, Bd, Qd, Rd;
    slong * Abounds, * Bbounds, * Qbounds, * Edegs;
    slong prod_deg;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(0 < ctx->minfo->nvars);
    FLINT_ASSERT(ctx->minfo->nvars <= FLINT_BITS);

    TMP_START;

    /*
        for each variable v we need to pack to degree deg_v(A)
        except for the outermost variable
    */
    Abounds = TMP_ARRAY_ALLOC(4*nvars, slong);
    Bbounds = Abounds + nvars;
    Qbounds = Bbounds + nvars;
    Edegs   = Qbounds + nvars;

    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Abounds, maxAfields, ctx->minfo);
    mpoly_get_monomial_ui_unpacked_ffmpz((ulong *)Bbounds, maxBfields, ctx->minfo);

    prod_deg = 1;
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        /* if divides, expected degrees */
        Edegs[i] = Abounds[i] - Bbounds[i];

        if (Abounds[i] < Bbounds[i])
        {
            success = 0;
            fmpz_mod_mpoly_zero(Q, ctx);
            goto cleanup;
        }

        if (i != 0)
        {
            /* variable of index i is not the outermost */
            Qbounds[i] = Abounds[i] + 1;
            Bbounds[i] = Abounds[i] + 1;
        }
        else
        {
            /* variable of index i is the outermost */
            Qbounds[i] = Abounds[i] - Bbounds[i] + 1;
            Bbounds[i] = Bbounds[i] + 1;
        }

        if (z_add_checked(&Abounds[i], Abounds[i], 1) ||
            z_mul_checked(&prod_deg, prod_deg, Abounds[i]))
        {
            success = -1;
            fmpz_mod_mpoly_zero(Q, ctx);
            goto cleanup;
        }
    }

    _fmpz_mod_mpoly_init_dense_mock(Ad, A, Abounds, ctx);
    _fmpz_mod_mpoly_init_dense_mock(Bd, B, Bbounds, ctx);
    fmpz_mod_poly_init(Qd, ctx->ffinfo);
    fmpz_mod_poly_init(Rd, ctx->ffinfo);
    fmpz_mod_poly_divrem(Qd, Rd, Ad, Bd, ctx->ffinfo);
    success = fmpz_mod_poly_is_zero(Rd, ctx->ffinfo) &&
              _from_dense(Q, A->bits, Qbounds, Edegs, Qd, ctx);
    fmpz_mod_poly_clear(Qd, ctx->ffinfo);
    fmpz_mod_poly_clear(Rd, ctx->ffinfo);
    _fmpz_mod_mpoly_clear_dense_mock(Ad);
    _fmpz_mod_mpoly_clear_dense_mock(Bd);

cleanup:

    TMP_END;
    return success;
}

/*
    return  -1 : function failed
             0 : B does not divide A
             1 : B divides A and quotient is in Q
*/
int fmpz_mod_mpoly_divides_dense(
    fmpz_mod_mpoly_t Q,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fmpz * maxAfields, * maxBfields;
    TMP_INIT;

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        if (!fmpz_mod_mpoly_is_zero(A, ctx) &&
            !fmpz_is_one(fmpz_mod_mpoly_ctx_modulus(ctx)))
        {
            flint_throw(FLINT_DIVZERO, "fmpz_mod_mpoly_divides_dense: divide by zero");
        }

        fmpz_mod_mpoly_zero(Q, ctx);
        return 1;
    }

    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        fmpz_mod_mpoly_zero(Q, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS ||
        ctx->minfo->nvars > FLINT_BITS ||
        ctx->minfo->nvars < 1)
    {
        return -1;
    }

    TMP_START;

    maxAfields = TMP_ARRAY_ALLOC(2*ctx->minfo->nfields, fmpz);
    maxBfields = maxAfields + ctx->minfo->nfields;
    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_init(maxAfields + i);

    mpoly_max_fields_fmpz(maxAfields, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);

    success = _fmpz_mod_mpoly_divides_dense_maxfields(Q,
                                            A, maxAfields, B, maxBfields, ctx);

    for (i = 0; i < 2*ctx->minfo->nfields; i++)
        fmpz_clear(maxAfields + i);

    TMP_END;
    return success;
}
