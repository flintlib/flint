/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"
#include "long_extras.h"

static void _clear_dense_mock(
    fmpz_mod_poly_t D)
{
    flint_free(D->coeffs);
}

static void _init_dense_mock(
    fmpz_mod_poly_t D,
    const fmpz_mod_mpolyn_t A,
    const slong * Adeg_bounds,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N, i, j, k, off, degb_prod;
    ulong * exps;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    degb_prod = 1;
    for (i = 0; i <= nvars; i++)
        degb_prod *= Adeg_bounds[i];

    D->alloc = degb_prod;
    D->coeffs = (fmpz *) flint_calloc(degb_prod, sizeof(fmpz));

    TMP_START;
    exps = TMP_ARRAY_ALLOC(nvars + 1, ulong);
    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);

    D->length = 0;
    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ui(exps, A->exps + N*i, A->bits, ctx->minfo);
        off = exps[0];
        for (j = 1; j < nvars; j++)
            off = exps[j] + Adeg_bounds[j]*off;

        off *= Adeg_bounds[nvars];

        D->length = FLINT_MAX(D->length, off + A->coeffs[i].length);
        FLINT_ASSERT(D->length <= D->alloc);
        for (k = 0; k < A->coeffs[i].length; k++)
            D->coeffs[off + k] = A->coeffs[i].coeffs[k]; /* shallow copy */
    }

    TMP_END;
}

/*
    Convert D to A if the degrees of A are <= expected_deg
    If not, return 0.
*/
static int _from_dense(
    fmpz_mod_mpolyn_t A,
    slong * Adeg_bounds,
    slong * expect_deg,
    fmpz_mod_poly_t D,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int ret;
    slong off, j, k;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong outrange;
    ulong * exps, * pcurexp, * pexps, * rangemask;
    TMP_INIT;

    FLINT_ASSERT(nvars + 1 <= FLINT_BITS);

    TMP_START;
    exps = (ulong *) TMP_ALLOC((nvars + 1)*sizeof(ulong));

    off = 1;
    for (j = 0; j <= nvars; j++)
    {
        FLINT_ASSERT(expect_deg[j] >= 0);
        off *= Adeg_bounds[j];
        exps[j] = expect_deg[j];
    }

    /* we are going to push back terms manually */
    A->length = 0;

    /* find exponent vector for all variables */
    pexps = TMP_ARRAY_ALLOC((N + 1)*(nvars + 1), ulong);
    for (k = 0; k < nvars; k++)
    {
        mpoly_gen_monomial_sp(pexps + (N + 1)*k, k, A->bits, ctx->minfo);
        (pexps + (N + 1)*k)[N] = 0;
    }
    mpoly_monomial_zero(pexps + (N + 1)*nvars, N);
    (pexps + (N + 1)*nvars)[N] = 1;

    /* get most significant exponent in exps and its vector in ptempexp */
    off--;
    pcurexp = TMP_ARRAY_ALLOC(N + 1, ulong);
    rangemask = TMP_ARRAY_ALLOC(nvars + 1, ulong);
    outrange = 0;
    mpoly_monomial_zero(pcurexp, N + 1);
    k = off;
    for (j = nvars; j >= 0; j--) 
    {
        exps[j] = k % Adeg_bounds[j];
        rangemask[j] = UWORD(1) << j;
        outrange ^= (outrange ^ FLINT_SIGN_EXT(expect_deg[j] - exps[j])) & rangemask[j];
        k = k / Adeg_bounds[j];
        mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + (N + 1)*j, N + 1);
    }

    /* scan down through the exponents */
    for (; off >= 0; off--)
    {
        if (off < D->length && !fmpz_is_zero(D->coeffs + off))
        {
            if (outrange)
            {
                ret = 0;
                goto cleanup;
            }

            if (A->length < 1 || !mpoly_monomial_equal(A->exps + N*(A->length - 1), pcurexp, N))
            {
                fmpz_mod_mpolyn_fit_length(A, A->length + 1, ctx);
                fmpz_mod_poly_zero(A->coeffs + A->length, ctx->ffinfo);
                fmpz_mod_poly_set_coeff_fmpz(A->coeffs + A->length, pcurexp[N], D->coeffs + off, ctx->ffinfo);
                mpoly_monomial_set(A->exps + N*A->length, pcurexp, N);
                A->length++;
            }
            else
            {
                fmpz_mod_poly_set_coeff_fmpz(A->coeffs + A->length - 1, pcurexp[N], D->coeffs + off, ctx->ffinfo);
            }
        }

        j = nvars;
        do {
            --exps[j];
            outrange ^= (outrange ^ FLINT_SIGN_EXT(expect_deg[j] - exps[j])) & rangemask[j];
            if (FLINT_SIGN_EXT(exps[j]) != 0)
            {
                FLINT_ASSERT(off == 0 || j > 0);
                FLINT_ASSERT(exps[j] == -UWORD(1));
                exps[j] = Adeg_bounds[j] - 1;
                outrange ^= (outrange ^ FLINT_SIGN_EXT(expect_deg[j] - exps[j])) & rangemask[j];
                mpoly_monomial_madd_inplace_mp(pcurexp, exps[j], pexps + (N + 1)*j, (N + 1));
            }
            else
            {
                mpoly_monomial_sub_mp(pcurexp, pcurexp, pexps + (N + 1)*j, N + 1);
                break;
            }
        } while (--j >= 0);
    }

    ret = 1;

cleanup:

    TMP_END;
    return ret;
}

/*
    return 1: yes it divides and Q=A/B
           0: no or dunno
*/
int fmpz_mod_mpolyn_divides(
    fmpz_mod_mpolyn_t Q,
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    flint_bitcnt_t bits = A->bits;
    slong i;
    slong nvars = ctx->minfo->nvars;
    fmpz_mod_poly_t Ad, Bd, Qd, Rd;
    slong * Abounds, * Bbounds, * Qbounds, * Edegs;
    slong prod_deg;
    TMP_INIT;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B->bits == bits);
    FLINT_ASSERT(Q->bits == bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (nvars + 1 > FLINT_BITS)
        return 0;

    TMP_START;

    /*
        for each variable v we need to pack to degree deg_v(A)
        except for the outermost variable
    */
    Abounds = TMP_ARRAY_ALLOC(4*(nvars + 1), slong);
    Bbounds = Abounds + nvars + 1;
    Qbounds = Bbounds + nvars + 1;
    Edegs   = Qbounds + nvars + 1;

    mpoly_degrees_si(Abounds, A->exps, A->length, bits, ctx->minfo);
    mpoly_degrees_si(Bbounds, B->exps, B->length, bits, ctx->minfo);
    Abounds[nvars] = fmpz_mod_mpolyn_lastdeg(A, ctx);
    Bbounds[nvars] = fmpz_mod_mpolyn_lastdeg(B, ctx);

    prod_deg = 1;
    for (i = 0; i <= nvars; i++)
    {
        /* if divides, expected degrees */
        Edegs[i] = Abounds[i] - Bbounds[i];

        if (Abounds[i] < Bbounds[i])
        {
            success = 0;
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
            success = 0;
            goto cleanup;
        }
    }

    _init_dense_mock(Ad, A, Abounds, ctx);
    _init_dense_mock(Bd, B, Bbounds, ctx);
    fmpz_mod_poly_init(Qd, ctx->ffinfo);
    fmpz_mod_poly_init(Rd, ctx->ffinfo);
    fmpz_mod_poly_divrem(Qd, Rd, Ad, Bd, ctx->ffinfo);
    success = fmpz_mod_poly_is_zero(Rd, ctx->ffinfo) &&
              _from_dense(Q, Qbounds, Edegs, Qd, ctx);
    fmpz_mod_poly_clear(Qd, ctx->ffinfo);
    fmpz_mod_poly_clear(Rd, ctx->ffinfo);
    _clear_dense_mock(Ad);
    _clear_dense_mock(Bd);

cleanup:

    TMP_END;
    return success;
}

