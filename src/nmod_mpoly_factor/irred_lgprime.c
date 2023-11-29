/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "fq_nmod_mpoly_factor.h"


static void _fq_nmod_mpoly_set_nmod_mpoly(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t Actx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t Bctx)
{
    slong d = fq_nmod_ctx_degree(Actx->fqctx);
    slong N = mpoly_words_per_exp(B->bits, Bctx->minfo);
    slong i;

    FLINT_ASSERT(Actx->minfo->ord == Bctx->minfo->ord);
    FLINT_ASSERT(Actx->minfo->nvars == Bctx->minfo->nvars);

    fq_nmod_mpoly_fit_length_reset_bits(A, B->length, B->bits, Actx);
    A->length = B->length;

    mpoly_copy_monomials(A->exps, B->exps, B->length, N);

    for (i = 0; i < B->length; i++)
        _n_fq_set_nmod(A->coeffs + d*i, B->coeffs[i], d);
}


static void _frob_combine(
    nmod_mpolyv_t Af,
    fq_nmod_mpolyv_t eAf,
    const nmod_mpoly_ctx_t ctx,
    const fq_nmod_mpoly_ctx_t ectx)
{
    slong d = fq_nmod_ctx_degree(ectx->fqctx);
    slong i, j, N;
    fq_nmod_mpolyv_t tfac;
    fq_nmod_mpoly_t t;
    nmod_mpoly_struct * s;

    FLINT_ASSERT(d > 1);

    fq_nmod_mpoly_init(t, ectx);
    fq_nmod_mpolyv_init(tfac, ectx);

    Af->length = 0;
    while (eAf->length > 0)
    {
        eAf->length--;
        fq_nmod_mpoly_swap(t, eAf->coeffs + eAf->length, ectx);

        fq_nmod_mpolyv_fit_length(tfac, 1, ectx);
        fq_nmod_mpoly_set(tfac->coeffs + 0, t, ectx);
        tfac->length = 1;

        for (i = 1; i < d; i++)
        {
            for (j = 0; j < t->length; j++)
            {
                n_fq_pow_ui(t->coeffs + d*j, t->coeffs + d*j,
                                     ectx->fqctx->modulus->mod.n, ectx->fqctx);
            }

            for (j = 0; j < eAf->length; j++)
            {
                if (fq_nmod_mpoly_equal(t, eAf->coeffs + j, ectx))
                    break;
            }

            if (j >= eAf->length)
                continue;   /* t should already be in tfac */

            fq_nmod_mpolyv_fit_length(tfac, tfac->length + 1, ectx);
            fq_nmod_mpoly_swap(tfac->coeffs + tfac->length, eAf->coeffs + j, ectx);
            tfac->length++;
            eAf->length--;
            fq_nmod_mpoly_swap(eAf->coeffs + j, eAf->coeffs + eAf->length, ectx);
        }

        fq_nmod_mpoly_swap(t, tfac->coeffs + 0, ectx);
        for (i = 1; i < tfac->length; i++)
            fq_nmod_mpoly_mul(t, t, tfac->coeffs + i, ectx);

        nmod_mpolyv_fit_length(Af, Af->length + 1, ctx);
        s = Af->coeffs + Af->length;
        Af->length++;

        nmod_mpoly_fit_length_reset_bits(s, t->length, t->bits, ctx);
        s->length = t->length;
        N = mpoly_words_per_exp(t->bits, ectx->minfo);
        mpoly_copy_monomials(s->exps, t->exps, t->length, N);
        for (i = 0; i < t->length; i++)
        {
            for (j = 1; j < d; j++)
            {
                if ((t->coeffs + d*i)[j] != 0)
                {
                    flint_throw(FLINT_ERROR, "fatal error in _frob_combine");
                }
            }
            s->coeffs[i] = (t->coeffs + d*i)[0];
        }
    }

    fq_nmod_mpolyv_clear(tfac, ectx);
    fq_nmod_mpoly_clear(t, ectx);
}

/*
    return:
        1: success
        0: failed, ran out of primes, don't try again
       -1: failed, don't try again
*/
int nmod_mpoly_factor_irred_lgprime_zassenhaus(
    nmod_mpolyv_t Af,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    fq_nmod_mpolyv_t eAf;
    fq_nmod_mpoly_t eA;
    fq_nmod_mpoly_ctx_t ectx;
    slong edeg;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    edeg = 2;
    fq_nmod_mpoly_ctx_init_deg(ectx, ctx->minfo->nvars, ORD_LEX, ctx->mod.n, edeg);
    fq_nmod_mpoly_init(eA, ectx);
    fq_nmod_mpolyv_init(eAf, ectx);

    goto have_prime;

choose_prime:

    edeg++;
    fq_nmod_mpoly_ctx_change_modulus(ectx, edeg);

have_prime:

    _fq_nmod_mpoly_set_nmod_mpoly(eA, ectx, A, ctx);

    success = fq_nmod_mpoly_factor_irred_smprime_zassenhaus(eAf, eA, ectx, state);
    if (success == 0)
        goto choose_prime;
    if (success < 0)
        goto cleanup;

    _frob_combine(Af, eAf, ctx, ectx);

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(eA, ectx);
    fq_nmod_mpolyv_clear(eAf, ectx);
    fq_nmod_mpoly_ctx_clear(ectx);

    return success;
}


/* the methods in the extended context want an irreducible factorization */
static int _map_fac(
    fq_nmod_mpoly_factor_t eAfac,
    const fq_nmod_mpoly_ctx_t ectx,
    const nmod_mpoly_factor_t Afac,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    fq_nmod_mpoly_t t;
    fq_nmod_mpoly_factor_t tfac;

    fq_nmod_mpoly_init(t, ectx);
    fq_nmod_mpoly_factor_init(tfac, ectx);

    fq_nmod_set_ui(eAfac->constant, Afac->constant, ectx->fqctx);
    eAfac->num = 0;
    for (i = 0; i < Afac->num; i++)
    {
        _fq_nmod_mpoly_set_nmod_mpoly(t, ectx, Afac->poly + i, ctx);
        success = fq_nmod_mpoly_factor(tfac, t, ectx);
        if (!success)
            goto cleanup;

        FLINT_ASSERT(fq_nmod_is_one(tfac->constant, ectx->fqctx));

        fq_nmod_mpoly_factor_fit_length(eAfac, eAfac->num + tfac->num, ectx);
        for (j = 0; j < tfac->num; j++)
        {
            fq_nmod_mpoly_swap(eAfac->poly + eAfac->num, tfac->poly + j, ectx);
            fmpz_mul(eAfac->exp + eAfac->num, tfac->exp + j, Afac->exp + i);
            eAfac->num++;
        }
    }

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(t, ectx);
    fq_nmod_mpoly_factor_clear(tfac, ectx);

    return success;
}

int nmod_mpoly_factor_irred_lgprime_wang(
    nmod_mpolyv_t Af,
    const nmod_mpoly_t A,
    const nmod_mpoly_factor_t lcAfac,
    const nmod_mpoly_t lcA,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    fq_nmod_mpoly_factor_t elcAfac;
    fq_nmod_mpolyv_t eAf;
    fq_nmod_mpoly_t eA, elcA;
    fq_nmod_mpoly_ctx_t ectx;
    slong edeg;
    const slong n = ctx->minfo->nvars - 1;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    edeg = 1 + n_clog(A->length + 1, ctx->mod.n)/2;
    edeg = FLINT_MAX(2, edeg);

    fq_nmod_mpoly_ctx_init_deg(ectx, n + 1, ORD_LEX, ctx->mod.n, edeg);
    fq_nmod_mpoly_init(eA, ectx);
    fq_nmod_mpolyv_init(eAf, ectx);
    fq_nmod_mpoly_init(elcA, ectx);
    fq_nmod_mpoly_factor_init(elcAfac, ectx);
    fq_nmod_mpoly_factor_fit_length(elcAfac, lcAfac->num, ectx);
    elcAfac->num = lcAfac->num;

    goto have_prime;

choose_prime:

    edeg++;
    fq_nmod_mpoly_ctx_change_modulus(ectx, edeg);

have_prime:

    _fq_nmod_mpoly_set_nmod_mpoly(eA, ectx, A, ctx);
    _fq_nmod_mpoly_set_nmod_mpoly(elcA, ectx, lcA, ctx);
    _map_fac(elcAfac, ectx, lcAfac, ctx);
    success = fq_nmod_mpoly_factor_irred_smprime_wang(eAf, eA, elcAfac, elcA, ectx, state);
    if (success == 0)
        goto choose_prime;
    if (success < 0)
        goto cleanup;

    _frob_combine(Af, eAf, ctx, ectx);

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(eA, ectx);
    fq_nmod_mpolyv_clear(eAf, ectx);
    fq_nmod_mpoly_clear(elcA, ectx);
    fq_nmod_mpoly_factor_clear(elcAfac, ectx);
    fq_nmod_mpoly_ctx_clear(ectx);

    return success;
}

int nmod_mpoly_factor_irred_lgprime_zippel(
    nmod_mpolyv_t Af,
    const nmod_mpoly_t A,
    const nmod_mpoly_factor_t lcAfac,
    const nmod_mpoly_t lcA,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    fq_nmod_mpoly_factor_t elcAfac;
    fq_nmod_mpolyv_t eAf;
    fq_nmod_mpoly_t eA, elcA;
    fq_nmod_mpoly_ctx_t ectx;
    slong edeg;
    const slong n = ctx->minfo->nvars - 1;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    edeg = 1 + n_clog(A->length + 1, ctx->mod.n);
    edeg = FLINT_MAX(2, edeg);
    fq_nmod_mpoly_ctx_init_deg(ectx, n + 1, ORD_LEX, ctx->mod.n, edeg);
    fq_nmod_mpoly_init(eA, ectx);
    fq_nmod_mpolyv_init(eAf, ectx);
    fq_nmod_mpoly_init(elcA, ectx);
    fq_nmod_mpoly_factor_init(elcAfac, ectx);
    fq_nmod_mpoly_factor_fit_length(elcAfac, lcAfac->num, ectx);
    elcAfac->num = lcAfac->num;

    goto have_prime;

choose_prime:

    edeg++;
    fq_nmod_mpoly_ctx_change_modulus(ectx, edeg);

have_prime:

    _fq_nmod_mpoly_set_nmod_mpoly(eA, ectx, A, ctx);
    _fq_nmod_mpoly_set_nmod_mpoly(elcA, ectx, lcA, ctx);
    fq_nmod_set_ui(elcAfac->constant, lcAfac->constant, ectx->fqctx);
    _map_fac(elcAfac, ectx, lcAfac, ctx);
    success = fq_nmod_mpoly_factor_irred_smprime_zippel(eAf, eA, elcAfac, elcA, ectx, state);
    if (success == 0)
        goto choose_prime;
    if (success < 0)
        goto cleanup;

    _frob_combine(Af, eAf, ctx, ectx);

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(eA, ectx);
    fq_nmod_mpolyv_clear(eAf, ectx);
    fq_nmod_mpoly_clear(elcA, ectx);
    fq_nmod_mpoly_factor_clear(elcAfac, ectx);
    fq_nmod_mpoly_ctx_clear(ectx);

    return success;
}

