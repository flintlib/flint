/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "nmod_mpoly_factor.h"

/* max = max abs coeffs(A) */
void fmpz_mpoly_height(
    fmpz_t max,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t t;

    fmpz_init(t);
    fmpz_zero(max);

    for (i = 0; i < A->length; i++)
    {
        fmpz_abs(t, A->coeffs + i);
        if (fmpz_cmp(max, t) < 0)
            fmpz_set(max, t);
    }

    fmpz_clear(t);
}

/* max = max abs coeffs(A), sum = sum abs coeffs(A) */
void fmpz_mpoly_heights(
    fmpz_t max,
    fmpz_t sum,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t t;

    fmpz_init(t);
    fmpz_zero(max);
    fmpz_zero(sum);

    for (i = 0; i < A->length; i++)
    {
        fmpz_abs(t, A->coeffs + i);
        fmpz_add(sum, sum, t);
        if (fmpz_cmp(max, t) < 0)
            fmpz_set(max, t);
    }

    fmpz_clear(t);
}


/* The inputs A and B are modified */
int fmpz_mpolyl_gcd_brown(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    fmpz_mpoly_t A,
    fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const mpoly_gcd_info_t I)
{
    int success;
    fmpz_t bound;
    slong offset, shift;
    mp_limb_t p, gammared;
    fmpz_t gamma, modulus;
    fmpz_t gnm, gns, anm, ans, bnm, bns;
    fmpz_t cA, cB, cG, cAbar, cBbar;
    fmpz_t temp;
    fmpz_mpoly_t T;
    nmod_mpolyn_t Gp, Abarp, Bbarp, Ap, Bp;
    nmod_mpoly_ctx_t pctx;
    flint_bitcnt_t bits = G->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    nmod_poly_stack_t Sp;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    mpoly_gen_offset_shift_sp(&offset, &shift, ctx->minfo->nvars - 1, bits, ctx->minfo);

    fmpz_init(gamma);
    fmpz_init(gnm);
    fmpz_init(gns);
    fmpz_init(anm);
    fmpz_init(ans);
    fmpz_init(bnm);
    fmpz_init(bns);
    fmpz_init(bound);
    fmpz_init(temp);
    fmpz_init_set_si(modulus, 1);

    /* compute contents of G, Abar, Bbar, A, B */
    fmpz_init(cA);
    fmpz_init(cB);
    fmpz_init(cG);
    fmpz_init(cAbar);
    fmpz_init(cBbar);
    _fmpz_vec_content(cA, A->coeffs, A->length);
    _fmpz_vec_content(cB, B->coeffs, B->length);
    fmpz_gcd(cG, cA, cB);
    fmpz_divexact(cAbar, cA, cG);
    fmpz_divexact(cBbar, cB, cG);

    /* remove content from inputs */
    fmpz_mpoly_scalar_divexact_fmpz(A, A, cA, ctx);
    fmpz_mpoly_scalar_divexact_fmpz(B, B, cB, ctx);

    fmpz_gcd(gamma, A->coeffs + 0, B->coeffs + 0);

    fmpz_mpoly_height(bound, A, ctx);
    fmpz_mpoly_height(temp, B, ctx);
    if (fmpz_cmp(bound, temp) < 0)
        fmpz_swap(bound, temp);
    fmpz_mul(bound, bound, gamma);
    fmpz_add(bound, bound, bound);

    fmpz_mpoly_init3(T, 0, bits, ctx);

    nmod_mpoly_ctx_init(pctx, ctx->minfo->nvars, ORD_LEX, 2);
    nmod_poly_stack_init(Sp, bits, pctx);
    nmod_mpolyn_init(Ap, bits, pctx);
    nmod_mpolyn_init(Bp, bits, pctx);
    nmod_mpolyn_init(Gp, bits, pctx);
    nmod_mpolyn_init(Abarp, bits, pctx);
    nmod_mpolyn_init(Bbarp, bits, pctx);

    p = UWORD(1) << (FLINT_BITS - 1);

choose_prime:

    if (p >= UWORD_MAX_PRIME)
    {
        /* ran out of primes - absolute failure */
        success = 0;
        goto cleanup;
    }
    p = n_nextprime(p, 1);

    /* make sure reduction does not kill both lc(A) and lc(B) */
    gammared = fmpz_fdiv_ui(gamma, p);
    if (gammared == 0)
    {
        goto choose_prime;
    }

    nmod_mpoly_ctx_set_modulus(pctx, p);
    /* the unfortunate nmod poly's store their own context :( */
    nmod_poly_stack_set_ctx(Sp, pctx);
    nmod_mpolyn_set_mod(Ap, pctx->mod);
    nmod_mpolyn_set_mod(Bp, pctx->mod);
    nmod_mpolyn_set_mod(Gp, pctx->mod);
    nmod_mpolyn_set_mod(Abarp, pctx->mod);
    nmod_mpolyn_set_mod(Bbarp, pctx->mod);

    /* reduction should kill neither A nor B */
    fmpz_mpoly_interp_reduce_p_mpolyn(Ap, pctx, A, ctx);
    fmpz_mpoly_interp_reduce_p_mpolyn(Bp, pctx, B, ctx);
    FLINT_ASSERT(Ap->length > 0);
    FLINT_ASSERT(Bp->length > 0);

    success = nmod_mpolyn_gcd_brown_smprime(Gp, Abarp, Bbarp,
                                   Ap, Bp, ctx->minfo->nvars - 1, pctx, I, Sp);
    if (!success)
    {
        goto choose_prime;
    }

    if (nmod_mpolyn_is_nonzero_nmod(Gp, pctx))
    {
        fmpz_mpoly_one(G, ctx);
        fmpz_mpoly_swap(Abar, A, ctx);
        fmpz_mpoly_swap(Bbar, B, ctx);
        goto successful_put_content;
    }

    if (!fmpz_is_one(modulus))
    {
        int cmp;
        slong k;
        FLINT_ASSERT(G->length > 0);

        k = n_poly_degree(Gp->coeffs + 0);
        cmp = mpoly_monomial_cmp_nomask_extra(G->exps + N*0,
                                        Gp->exps + N*0, N, offset, k << shift);
        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            fmpz_one(modulus);
        }
    }

    FLINT_ASSERT(1 == nmod_mpolyn_leadcoeff(Gp, pctx));
    nmod_mpolyn_scalar_mul_nmod(Gp, gammared, pctx);

    if (!fmpz_is_one(modulus))
    {
        fmpz_mpoly_interp_crt_p_mpolyn(G, T, ctx, modulus, Gp, pctx);
        fmpz_mpoly_interp_crt_p_mpolyn(Abar, T, ctx, modulus, Abarp, pctx);
        fmpz_mpoly_interp_crt_p_mpolyn(Bbar, T, ctx, modulus, Bbarp, pctx);
    }
    else
    {
        fmpz_mpoly_interp_lift_p_mpolyn(G, ctx, Gp, pctx);
        fmpz_mpoly_interp_lift_p_mpolyn(Abar, ctx, Abarp, pctx);
        fmpz_mpoly_interp_lift_p_mpolyn(Bbar, ctx, Bbarp, pctx);
    }

    fmpz_mul_ui(modulus, modulus, p);

    if (fmpz_cmp(modulus, bound) <= 0)
    {
        goto choose_prime;
    }

    fmpz_mpoly_heights(gnm, gns, G, ctx);
    fmpz_mpoly_heights(anm, ans, Abar, ctx);
    fmpz_mpoly_heights(bnm, bns, Bbar, ctx);
    fmpz_mul(ans, ans, gnm);
    fmpz_mul(anm, anm, gns);
    fmpz_mul(bns, bns, gnm);
    fmpz_mul(bnm, bnm, gns);

    if (fmpz_cmp(ans, anm) > 0)
        fmpz_swap(ans, anm);
    if (fmpz_cmp(bns, bnm) > 0)
        fmpz_swap(bns, bnm);
    fmpz_add(ans, ans, ans);
    fmpz_add(bns, bns, bns);
    if (fmpz_cmp(ans, modulus) < 0 && fmpz_cmp(bns, modulus) < 0)
    {
        goto successful;
    }

    /* do not reset modulus to 1 */
    goto choose_prime;

successful:

    FLINT_ASSERT(fmpz_equal(gamma, G->coeffs + 0));

    _fmpz_vec_content(temp, G->coeffs, G->length);
    fmpz_mpoly_scalar_divexact_fmpz(G, G, temp, ctx);
    fmpz_mpoly_scalar_divexact_fmpz(Abar, Abar, G->coeffs + 0, ctx);
    fmpz_mpoly_scalar_divexact_fmpz(Bbar, Bbar, G->coeffs + 0, ctx);

successful_put_content:

    fmpz_mpoly_scalar_mul_fmpz(G, G, cG, ctx);
    fmpz_mpoly_scalar_mul_fmpz(Abar, Abar, cAbar, ctx);
    fmpz_mpoly_scalar_mul_fmpz(Bbar, Bbar, cBbar, ctx);

    success = 1;

cleanup:

    fmpz_clear(cA);
    fmpz_clear(cB);
    fmpz_clear(cG);
    fmpz_clear(cAbar);
    fmpz_clear(cBbar);

    fmpz_clear(gamma);
    fmpz_clear(gnm);
    fmpz_clear(gns);
    fmpz_clear(anm);
    fmpz_clear(ans);
    fmpz_clear(bnm);
    fmpz_clear(bns);
    fmpz_clear(bound);
    fmpz_clear(temp);
    fmpz_clear(modulus);

    nmod_mpolyn_clear(Gp, pctx);
    nmod_mpolyn_clear(Abarp, pctx);
    nmod_mpolyn_clear(Bbarp, pctx);
    nmod_mpolyn_clear(Ap, pctx);
    nmod_mpolyn_clear(Bp, pctx);
    nmod_poly_stack_clear(Sp);
    nmod_mpoly_ctx_clear(pctx);

    fmpz_mpoly_clear(T, ctx);

    return success;
}

