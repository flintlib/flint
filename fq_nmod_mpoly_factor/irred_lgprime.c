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


void bad_n_fq_embed_sm_to_lg(
    mp_limb_t * out_,            /* element of lgctx */
    const n_poly_t in_,  /* poly over smctx */
    const bad_fq_nmod_embed_t emb)
{
    fq_nmod_t out;
    fq_nmod_poly_t in;

    fq_nmod_init(out, emb->lgctx);
    fq_nmod_poly_init(in, emb->smctx);

    n_poly_fq_get_fq_nmod_poly(in, in_, emb->smctx);
    bad_fq_nmod_embed_sm_to_lg(out, in , emb);
    n_fq_set_fq_nmod(out_, out, emb->lgctx);

    fq_nmod_poly_clear(in, emb->smctx);
    fq_nmod_clear(out, emb->lgctx);
}

void bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(
    fq_nmod_t out,            /* element of lgctx */
    const n_poly_t in_,  /* poly over smctx */
    const bad_fq_nmod_embed_t emb)
{
    fq_nmod_poly_t in;

    fq_nmod_poly_init(in, emb->smctx);

    n_poly_fq_get_fq_nmod_poly(in, in_, emb->smctx);
    bad_fq_nmod_embed_sm_to_lg(out, in , emb);

    fq_nmod_poly_clear(in, emb->smctx);
}

void bad_n_fq_embed_lg_to_sm(
    n_poly_t out_,  /* poly over smctx */
    const mp_limb_t * in_,  /* element of lgctx */
    const bad_fq_nmod_embed_t emb)
{
    fq_nmod_poly_t out;
    fq_nmod_t in;

    fq_nmod_poly_init(out, emb->smctx);
    fq_nmod_init(in, emb->lgctx);

    n_fq_get_fq_nmod(in, in_, emb->lgctx);
    bad_fq_nmod_embed_lg_to_sm(out, in, emb);
    n_poly_fq_set_fq_nmod_poly(out_, out, emb->smctx);

    fq_nmod_clear(in, emb->lgctx);
    fq_nmod_poly_clear(out, emb->smctx);
}


static void _map_sm_to_lg(
    fq_nmod_mpoly_t eA,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ectx,
    const fq_nmod_mpoly_ctx_t ctx,
    const bad_fq_nmod_embed_t emb)
{
    slong i;
    fq_nmod_poly_t t;

    FLINT_ASSERT(eA != A);

    fq_nmod_poly_init(t, emb->smctx);

    fq_nmod_mpoly_fit_length_set_bits(eA, A->length, A->bits, ectx);
    mpoly_copy_monomials(eA->exps, A->exps,
                         mpoly_words_per_exp(A->bits, ectx->minfo), A->length);
    for (i = 0; i < A->length; i++)
    {
        fq_nmod_poly_set_fq_nmod(t, A->coeffs + i, emb->smctx);
        bad_fq_nmod_embed_sm_to_lg(eA->coeffs + i, t, emb);
    }
    eA->length = A->length;

    fq_nmod_poly_clear(t, emb->smctx);
}

static void _frob_combine(
    fq_nmod_mpolyv_t Af,
    fq_nmod_mpolyv_t eAf,
    const fq_nmod_mpoly_ctx_t ctx,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    slong i, j;
    fq_nmod_mpolyv_t tfac;
    fq_nmod_mpoly_t t;
    fq_nmod_mpoly_struct * s;
    fq_nmod_poly_t c;
    slong k = fq_nmod_ctx_degree(ectx->fqctx)/fq_nmod_ctx_degree(ctx->fqctx);
    fmpz_t q;

    FLINT_ASSERT(k > 1);

    fmpz_init(q);
    fq_nmod_mpoly_init(t, ectx);
    fq_nmod_mpolyv_init(tfac, ectx);
    fq_nmod_poly_init(c, ctx->fqctx);

    fmpz_pow_ui(q, fq_nmod_ctx_prime(ctx->fqctx), fq_nmod_ctx_degree(ctx->fqctx));

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
                fq_nmod_pow(t->coeffs + j, t->coeffs + j, q, ectx->fqctx);

            for (j = 0; j < eAf->length; j++)
            {
                if (fq_nmod_mpoly_equal(t, eAf->coeffs + j, ectx))
                    break;
            }

            if (j >= eAf->length)
                continue;   /* t, should already be in tfac */

            fq_nmod_mpolyv_fit_length(tfac, tfac->length + 1, ectx);
            fq_nmod_mpoly_swap(tfac->coeffs + tfac->length, eAf->coeffs + j, ectx);
            tfac->length++;
            eAf->length--;
            fq_nmod_mpoly_swap(eAf->coeffs + j, eAf->coeffs + eAf->length, ectx);
        }

        fq_nmod_mpoly_swap(t, tfac->coeffs + 0, ectx);
        for (i = 1; i < tfac->length; i++)
            fq_nmod_mpoly_mul(t, t, tfac->coeffs + i, ectx);

        fq_nmod_mpolyv_fit_length(Af, Af->length + 1, ctx);
        s = Af->coeffs + Af->length;
        Af->length++;

        fq_nmod_mpoly_fit_length_set_bits(s, t->length, t->bits, ctx);
        s->length = t->length;
        mpoly_copy_monomials(s->exps, t->exps,
                         mpoly_words_per_exp(t->bits, ectx->minfo), t->length);
        for (i = 0; i < t->length; i++)
        {
            bad_fq_nmod_embed_lg_to_sm(c, t->coeffs + i, emb);
            if (c->length != 1)
            {
                flint_printf("fatal error in _frob_combine");
                flint_abort();
            }
            fq_nmod_swap(s->coeffs + i, c->coeffs + 0, ctx->fqctx);
        }
    }

    fq_nmod_poly_clear(c, ctx->fqctx);
    fq_nmod_mpolyv_clear(tfac, ectx);
    fq_nmod_mpoly_clear(t, ectx);
    fmpz_clear(q);
}

/*
    return:
        1: success
        0: failed, ran out of primes, don't try again
       -1: failed, don't try again
*/
int fq_nmod_mpoly_factor_irred_lgprime_zassenhaus(
    fq_nmod_mpolyv_t Af,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    fq_nmod_mpolyv_t eAf;
    fq_nmod_mpoly_t eA;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx;    

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, ectx, ctx, state);

    fq_nmod_mpoly_init(eA, ectx);
    fq_nmod_mpolyv_init(eAf, ectx);

    goto have_prime;

choose_prime:

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, state);
    if (cur_emb == NULL)
    {
        success = 0;
        goto cleanup;
    }

have_prime:

    _map_sm_to_lg(eA, A, ectx, ctx, cur_emb);

    success = fq_nmod_mpoly_factor_irred_smprime_zassenhaus(eAf, eA, ectx, state);
    if (success == 0)
        goto choose_prime;
    if (success < 0)
        goto cleanup;

    _frob_combine(Af, eAf, ctx, ectx, cur_emb);

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(eA, ectx);
    fq_nmod_mpolyv_clear(eAf, ectx);
    bad_fq_nmod_mpoly_embed_chooser_clear(embc, ectx, ctx, state);

    return success;
}

int fq_nmod_mpoly_factor_irred_lgprime_wang(
    fq_nmod_mpolyv_t Af,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    fq_nmod_mpoly_factor_t elcAfac;
    fq_nmod_mpolyv_t eAf;
    fq_nmod_mpoly_t eA, elcA;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx;
    slong i;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(fq_nmod_is_one(A->coeffs + 0, ctx->fqctx));
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, ectx, ctx, state);

    fq_nmod_mpoly_init(eA, ectx);
    fq_nmod_mpolyv_init(eAf, ectx);
    fq_nmod_mpoly_init(elcA, ectx);
    fq_nmod_mpoly_factor_init(elcAfac, ectx);
    fq_nmod_mpoly_factor_fit_length(elcAfac, lcAfac->num, ectx);
    elcAfac->num = lcAfac->num;

    goto have_prime;

choose_prime:

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, state);
    if (cur_emb == NULL)
    {
        success = 0;
        goto cleanup;
    }

have_prime:

    _map_sm_to_lg(eA, A, ectx, ctx, cur_emb);
    _map_sm_to_lg(elcA, lcA, ectx, ctx, cur_emb);
    {
    fq_nmod_poly_t tt;
    fq_nmod_poly_init(tt, ctx->fqctx);
    fq_nmod_poly_set_fq_nmod(tt, lcAfac->constant, ctx->fqctx);
    bad_fq_nmod_embed_sm_to_lg(elcAfac->constant, tt, cur_emb);
    fq_nmod_poly_clear(tt, ctx->fqctx);
    }
    for (i = 0; i < lcAfac->num; i++)
    {
        fmpz_set(elcAfac->exp + i, lcAfac->exp + i);
        _map_sm_to_lg(elcAfac->poly + i, lcAfac->poly + i, ectx, ctx, cur_emb);
    }

    success = fq_nmod_mpoly_factor_irred_smprime_wang(eAf, eA, elcAfac, elcA, ectx, state);
    if (success == 0)
        goto choose_prime;
    if (success < 0)
        goto cleanup;

    _frob_combine(Af, eAf, ctx, ectx, cur_emb);

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(eA, ectx);
    fq_nmod_mpolyv_clear(eAf, ectx);
    fq_nmod_mpoly_clear(elcA, ectx);
    fq_nmod_mpoly_factor_clear(elcAfac, ectx);

    bad_fq_nmod_mpoly_embed_chooser_clear(embc, ectx, ctx, state);

    return success;
}

int fq_nmod_mpoly_factor_irred_lgprime_zippel(
    fq_nmod_mpolyv_t Af,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    fq_nmod_mpoly_factor_t elcAfac;
    fq_nmod_mpolyv_t eAf;
    fq_nmod_mpoly_t eA, elcA;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx;
    slong i;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(fq_nmod_is_one(A->coeffs + 0, ctx->fqctx));
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, ectx, ctx, state);

    fq_nmod_mpoly_init(eA, ectx);
    fq_nmod_mpolyv_init(eAf, ectx);
    fq_nmod_mpoly_init(elcA, ectx);
    fq_nmod_mpoly_factor_init(elcAfac, ectx);
    fq_nmod_mpoly_factor_fit_length(elcAfac, lcAfac->num, ectx);
    elcAfac->num = lcAfac->num;

    goto have_prime;

choose_prime:

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, state);
    if (cur_emb == NULL)
    {
        success = 0;
        goto cleanup;
    }

have_prime:

    _map_sm_to_lg(eA, A, ectx, ctx, cur_emb);
    _map_sm_to_lg(elcA, lcA, ectx, ctx, cur_emb);
    {
    fq_nmod_poly_t tt;
    fq_nmod_poly_init(tt, ctx->fqctx);
    fq_nmod_poly_set_fq_nmod(tt, lcAfac->constant, ctx->fqctx);
    bad_fq_nmod_embed_sm_to_lg(elcAfac->constant, tt, cur_emb);
    fq_nmod_poly_clear(tt, ctx->fqctx);
    }
    for (i = 0; i < lcAfac->num; i++)
    {
        fmpz_set(elcAfac->exp + i, lcAfac->exp + i);
        _map_sm_to_lg(elcAfac->poly + i, lcAfac->poly + i, ectx, ctx, cur_emb);
    }

    success = fq_nmod_mpoly_factor_irred_smprime_zippel(eAf, eA, elcAfac, elcA, ectx, state);
    if (success == 0)
        goto choose_prime;
    if (success < 0)
        goto cleanup;

    _frob_combine(Af, eAf, ctx, ectx, cur_emb);

    success = 1;

cleanup:

    fq_nmod_mpoly_clear(eA, ectx);
    fq_nmod_mpolyv_clear(eAf, ectx);
    fq_nmod_mpoly_clear(elcA, ectx);
    fq_nmod_mpoly_factor_clear(elcAfac, ectx);

    bad_fq_nmod_mpoly_embed_chooser_clear(embc, ectx, ctx, state);

    return success;
}

