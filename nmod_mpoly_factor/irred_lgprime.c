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
    slong i, N;

    FLINT_ASSERT(Actx->minfo->ord == Bctx->minfo->ord);
    FLINT_ASSERT(Actx->minfo->nvars == Bctx->minfo->nvars);

    fq_nmod_mpoly_fit_bits(A, B->bits, Actx);
    A->bits = B->bits;

    N = mpoly_words_per_exp(B->bits, Bctx->minfo);

    fq_nmod_mpoly_fit_length(A, B->length, Actx);
    A->length = B->length;

    mpoly_copy_monomials(A->exps, B->exps, B->length, N);

    for (i = 0; i < B->length; i++)
        fq_nmod_set_ui(A->coeffs + i, B->coeffs[i], Actx->fqctx);
}


static void _frob_combine(
    nmod_mpolyv_t Af,
    fq_nmod_mpolyv_t eAf,
    const nmod_mpoly_ctx_t ctx,
    const fq_nmod_mpoly_ctx_t ectx)
{
    slong i, j;
    fq_nmod_mpolyv_t tfac;
    fq_nmod_mpoly_t t;
    nmod_mpoly_struct * s;
    slong k = fq_nmod_ctx_degree(ectx->fqctx);

    FLINT_ASSERT(k > 1);

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

        for (i = 1; i < k; i++)
        {
            for (j = 0; j < t->length; j++)
            {
                fq_nmod_pow_ui(t->coeffs + j, t->coeffs + j,
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

        nmod_mpoly_fit_length_set_bits(s, t->length, t->bits, ctx);
        s->length = t->length;
        mpoly_copy_monomials(s->exps, t->exps,
                         mpoly_words_per_exp(t->bits, ectx->minfo), t->length);
        for (i = 0; i < t->length; i++)
        {
            if (t->coeffs[i].length != 1)
            {
                flint_printf("fatal error in _frob_combine");
                flint_abort();
            }
            s->coeffs[i] = t->coeffs[i].coeffs[0];
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
    fq_nmod_mpoly_ctx_init_deg(ectx, ctx->minfo->nvars, ORD_LEX, ctx->ffinfo->mod.n, edeg);
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
    slong i;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    edeg = 1 + n_clog(A->length + 1, ctx->ffinfo->mod.n)/2;
    edeg = FLINT_MAX(2, edeg);

    fq_nmod_mpoly_ctx_init_deg(ectx, n + 1, ORD_LEX, ctx->ffinfo->mod.n, edeg);
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
    for (i = 0; i < lcAfac->num; i++)
    {
        fmpz_set(elcAfac->exp + i, lcAfac->exp + i);
        _fq_nmod_mpoly_set_nmod_mpoly(elcAfac->poly + i, ectx, lcAfac->poly + i, ctx);
    }

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
    slong i;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    edeg = 1 + n_clog(A->length + 1, ctx->ffinfo->mod.n);
    edeg = FLINT_MAX(2, edeg);
    fq_nmod_mpoly_ctx_init_deg(ectx, n + 1, ORD_LEX, ctx->ffinfo->mod.n, edeg);
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
    for (i = 0; i < lcAfac->num; i++)
    {
        fmpz_set(elcAfac->exp + i, lcAfac->exp + i);
        _fq_nmod_mpoly_set_nmod_mpoly(elcAfac->poly + i, ectx, lcAfac->poly + i, ctx);
    }

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

