/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"


static void _fmpz_mod_poly_gcd_cofactors(
    fmpz_mod_poly_t G,
    fmpz_mod_poly_t Abar,
    fmpz_mod_poly_t Bbar,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_t t)
{
    fmpz_mod_poly_gcd(G, A, B, ctx);
    fmpz_mod_poly_divrem(Abar, t, A, G, ctx);
    fmpz_mod_poly_divrem(Bbar, t, B, G, ctx);
}


int fmpz_mod_polyu1n_gcd_brown_smprime(
    fmpz_mod_polyun_t G,
    fmpz_mod_polyun_t Abar,
    fmpz_mod_polyun_t Bbar,
    fmpz_mod_polyun_t A,
    fmpz_mod_polyun_t B,
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_stack_t St_poly,
    fmpz_mod_polyun_stack_t St_polyun)
{
    int success;
    slong bound;
    fmpz_t alpha, temp, gammaevalp, gammaevalm;
    fmpz_mod_poly_struct * Aevalp, * Bevalp, * Gevalp, * Abarevalp, * Bbarevalp;
    fmpz_mod_poly_struct * Aevalm, * Bevalm, * Gevalm, * Abarevalm, * Bbarevalm;
    fmpz_mod_polyun_struct * T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fmpz_mod_poly_struct * cA, * cB, * cG, * cAbar, * cBbar, * gamma;
    fmpz_mod_poly_struct * modulus, * alphapow, * r;
    int gstab, astab, bstab, use_stab;
#if FLINT_WANT_ASSERT
    fmpz_mod_poly_t leadA, leadB;
#endif

    fmpz_init(alpha);
    fmpz_init(temp);
    fmpz_init(gammaevalp);
    fmpz_init(gammaevalm);

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_init(leadA, ctx);
    fmpz_mod_poly_init(leadB, ctx);
    fmpz_mod_poly_set(leadA, A->coeffs + 0, ctx);
    fmpz_mod_poly_set(leadB, B->coeffs + 0, ctx);
#endif

    fmpz_mod_poly_stack_fit_request(St_poly, 19);
    cA          = fmpz_mod_poly_stack_take_top(St_poly);
    cB          = fmpz_mod_poly_stack_take_top(St_poly);
    cG          = fmpz_mod_poly_stack_take_top(St_poly);
    cAbar       = fmpz_mod_poly_stack_take_top(St_poly);
    cBbar       = fmpz_mod_poly_stack_take_top(St_poly);
    gamma       = fmpz_mod_poly_stack_take_top(St_poly);
    Aevalp      = fmpz_mod_poly_stack_take_top(St_poly);
    Bevalp      = fmpz_mod_poly_stack_take_top(St_poly);
    Gevalp      = fmpz_mod_poly_stack_take_top(St_poly);
    Abarevalp   = fmpz_mod_poly_stack_take_top(St_poly);
    Bbarevalp   = fmpz_mod_poly_stack_take_top(St_poly);
    Aevalm      = fmpz_mod_poly_stack_take_top(St_poly);
    Bevalm      = fmpz_mod_poly_stack_take_top(St_poly);
    Gevalm      = fmpz_mod_poly_stack_take_top(St_poly);
    Abarevalm   = fmpz_mod_poly_stack_take_top(St_poly);
    Bbarevalm   = fmpz_mod_poly_stack_take_top(St_poly);
    r           = fmpz_mod_poly_stack_take_top(St_poly);
    alphapow    = fmpz_mod_poly_stack_take_top(St_poly);
    modulus     = fmpz_mod_poly_stack_take_top(St_poly);

    fmpz_mod_polyun_stack_fit_request(St_polyun, 1);
    T = fmpz_mod_polyun_stack_take_top(St_polyun);

    _fmpz_mod_poly_vec_remove_content(cA, A->coeffs, A->length, ctx);
    _fmpz_mod_poly_vec_remove_content(cB, B->coeffs, B->length, ctx);

    _fmpz_mod_poly_gcd_cofactors(cG, cAbar, cBbar, cA, cB, ctx, r);

    fmpz_mod_poly_gcd(gamma, A->coeffs + 0, B->coeffs + 0, ctx);

    ldegA = _fmpz_mod_poly_vec_max_degree(A->coeffs, A->length, ctx);
    ldegB = _fmpz_mod_poly_vec_max_degree(B->coeffs, B->length, ctx);
    deggamma = fmpz_mod_poly_degree(gamma, ctx);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fmpz_mod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1), ctx);
    fmpz_mod_poly_one(modulus, ctx);

    use_stab = 1;
    gstab = bstab = astab = 0;

    fmpz_fdiv_q_2exp(alpha, fmpz_mod_ctx_modulus(ctx), 1);

choose_prime:

    fmpz_sub_ui(alpha, alpha, 1);
    if (fmpz_sgn(alpha) <= 0)
    {
        success = 0;
        goto cleanup;
    }

    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    fmpz_one(alphapow->coeffs + 0);
    fmpz_set(alphapow->coeffs + 1, alpha);

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_eval2_pow(gammaevalp, gammaevalm, gamma, alphapow, ctx);
    if (fmpz_is_zero(gammaevalp) || fmpz_is_zero(gammaevalm))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fmpz_mod_polyu1n_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, ctx);
    fmpz_mod_polyu1n_interp_reduce_2sm_poly(Bevalp, Bevalm, B, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    if (use_stab && gstab)
    {
        slong Gdeg;
        fmpz_mod_polyu1n_interp_reduce_2sm_poly(Gevalp, Gevalm, G, alphapow, ctx);
        Gdeg = G->exps[0];
        success = 1;
        success = success && _fmpz_mod_poly_degree(Gevalp) == Gdeg;
        success = success && _fmpz_mod_poly_degree(Gevalm) == Gdeg;
        success = success && fmpz_equal(Gevalp->coeffs + Gdeg, gammaevalp);
        success = success && fmpz_equal(Gevalm->coeffs + Gdeg, gammaevalm);
        fmpz_mod_poly_divrem(Abarevalp, r, Aevalp, Gevalp, ctx);
        success = success && (r->length == 0);
        fmpz_mod_poly_divrem(Abarevalm, r, Aevalm, Gevalm, ctx);
        success = success && (r->length == 0);
        fmpz_mod_poly_divrem(Bbarevalp, r, Bevalp, Gevalp, ctx);
        success = success && (r->length == 0);
        fmpz_mod_poly_divrem(Bbarevalm, r, Bevalm, Gevalm, ctx);
        success = success && (r->length == 0);

        if (!success)
        {
            use_stab = 0;
            fmpz_mod_poly_one(modulus, ctx);
            goto choose_prime;
        }

        fmpz_mod_poly_scalar_mul_fmpz(Abarevalp, Abarevalp, gammaevalp, ctx);
        fmpz_mod_poly_scalar_mul_fmpz(Abarevalm, Abarevalm, gammaevalm, ctx);
        fmpz_mod_poly_scalar_mul_fmpz(Bbarevalp, Bbarevalp, gammaevalp, ctx);
        fmpz_mod_poly_scalar_mul_fmpz(Bbarevalm, Bbarevalm, gammaevalm, ctx);
    }
    else
    {
        _fmpz_mod_poly_gcd_cofactors(Gevalp, Abarevalp, Bbarevalp,
                                                       Aevalp, Bevalp, ctx, r);
        _fmpz_mod_poly_gcd_cofactors(Gevalm, Abarevalm, Bbarevalm,
                                                       Aevalm, Bevalm, ctx, r);
    }

    FLINT_ASSERT(Gevalp->length > 0);
    FLINT_ASSERT(Abarevalp->length > 0);
    FLINT_ASSERT(Bbarevalp->length > 0);
    FLINT_ASSERT(Gevalm->length > 0);
    FLINT_ASSERT(Abarevalm->length > 0);
    FLINT_ASSERT(Bbarevalm->length > 0);

    if (_fmpz_mod_poly_degree(Gevalp) == 0 || _fmpz_mod_poly_degree(Gevalm) == 0)
    {
        fmpz_mod_polyun_one(G, ctx);
        fmpz_mod_polyun_swap(Abar, A);
        fmpz_mod_polyun_swap(Bbar, B);
        goto successful_put_content;
    }

    if (_fmpz_mod_poly_degree(Gevalp) != _fmpz_mod_poly_degree(Gevalm))
    {
        goto choose_prime;
    }

    if (_fmpz_mod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (_fmpz_mod_poly_degree(Gevalp) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (_fmpz_mod_poly_degree(Gevalp) < G->exps[0])
        {
            fmpz_mod_poly_one(modulus, ctx);
        }
    }

    fmpz_mod_poly_scalar_mul_fmpz(Gevalp, Gevalp, gammaevalp, ctx);
    fmpz_mod_poly_scalar_mul_fmpz(Gevalm, Gevalm, gammaevalm, ctx);
    if (fmpz_mod_poly_degree(modulus, ctx) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(temp, modulus, alpha, ctx);
        fmpz_mod_mul(temp, temp, alpha, ctx);
        fmpz_mod_add(temp, temp, temp, ctx);
        fmpz_mod_poly_scalar_div_fmpz(modulus, modulus, temp, ctx);
        gstab = gstab || !fmpz_mod_polyu1n_interp_crt_2sm_poly(&ldegG, G, T, Gevalp, Gevalm, modulus, alphapow, ctx);
        fmpz_mod_polyu1n_interp_crt_2sm_poly(&ldegAbar, Abar, T, Abarevalp, Abarevalm, modulus, alphapow, ctx);
        fmpz_mod_polyu1n_interp_crt_2sm_poly(&ldegBbar, Bbar, T, Bbarevalp, Bbarevalm, modulus, alphapow, ctx);
    }
    else
    {
        fmpz_mod_polyu1n_interp_lift_2sm_poly(&ldegG, G, Gevalp, Gevalm, alpha, ctx);
        fmpz_mod_polyu1n_interp_lift_2sm_poly(&ldegAbar, Abar, Abarevalp, Abarevalm, alpha, ctx);
        fmpz_mod_polyu1n_interp_lift_2sm_poly(&ldegBbar, Bbar, Bbarevalp, Bbarevalm, alpha, ctx);
        gstab = astab = bstab = 0;
    }

    fmpz_mod_mul(temp, alpha, alpha, ctx);
    fmpz_mod_neg(temp, temp, ctx);
    fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 2, temp, ctx);

    if (fmpz_mod_poly_degree(modulus, ctx) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (deggamma + ldegA == ldegG + ldegAbar &&
        deggamma + ldegB == ldegG + ldegBbar)
    {
        goto successful;
    }

    fmpz_mod_poly_one(modulus, ctx);
    goto choose_prime;

successful:

    _fmpz_mod_poly_vec_remove_content(r, G->coeffs, G->length, ctx);
    _fmpz_mod_poly_vec_divexact_poly(Abar->coeffs, Abar->length, G->coeffs + 0, ctx);
    _fmpz_mod_poly_vec_divexact_poly(Bbar->coeffs, Bbar->length, G->coeffs + 0, ctx);

successful_put_content:

    _fmpz_mod_poly_vec_mul_poly(G->coeffs, G->length, cG, ctx);
    _fmpz_mod_poly_vec_mul_poly(Abar->coeffs, Abar->length, cAbar, ctx);
    _fmpz_mod_poly_vec_mul_poly(Bbar->coeffs, Bbar->length, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_poly_lead(G->coeffs + 0, ctx)));
        fmpz_mod_poly_mul(modulus, G->coeffs + 0, Abar->coeffs + 0, ctx);
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadA, ctx));
        fmpz_mod_poly_mul(modulus, G->coeffs + 0, Bbar->coeffs + 0, ctx);
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadB, ctx));
    }
    fmpz_mod_poly_clear(leadA, ctx);
    fmpz_mod_poly_clear(leadB, ctx);
#endif

    fmpz_mod_poly_stack_give_back(St_poly, 19);
    fmpz_mod_polyun_stack_give_back(St_polyun, 1);

    fmpz_clear(alpha);
    fmpz_clear(temp);
    fmpz_clear(gammaevalp);
    fmpz_clear(gammaevalm);

    return success;
}


static void fmpz_mod_mpolyn_set_polyun_swap(
    fmpz_mod_mpolyn_t A,
    fmpz_mod_polyun_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N, off, shift;
    flint_bitcnt_t bits = A->bits;
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, bits, ctx->minfo);

    fmpz_mod_mpolyn_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        mpoly_monomial_zero(A->exps + N*i, N);
        (A->exps + N*i)[off] = B->exps[i] << shift;
        fmpz_mod_poly_swap(A->coeffs + i, B->coeffs + i, ctx->ffinfo);
    }

    A->length = B->length;
}

static void fmpz_mod_mpolyn_get_polyun_swap(
    fmpz_mod_polyun_t B,
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N, off, shift;
    flint_bitcnt_t bits = A->bits;
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, bits, ctx->minfo);

    fmpz_mod_polyun_fit_length(B, A->length, ctx->ffinfo);
    for (i = 0; i < A->length; i++)
    {
        B->exps[i] = (A->exps + N*i)[off] >> shift;
        fmpz_mod_poly_swap(B->coeffs + i, A->coeffs + i, ctx->ffinfo);
    }

    B->length = A->length;
}


int fmpz_mod_mpolyn_gcd_brown_smprime(
    fmpz_mod_mpolyn_t G,
    fmpz_mod_mpolyn_t Abar,
    fmpz_mod_mpolyn_t Bbar,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx,
    const mpoly_gcd_info_t I,
    fmpz_mod_poly_polyun_mpolyn_stack_t St)
{
    int success, changed;
    slong bound, upper_limit;
    slong offset, shift;
    fmpz_t alpha, temp, gammaevalp, gammaevalm;
    fmpz_mod_mpolyn_struct * Aevalp, * Bevalp, * Gevalp, * Abarevalp, * Bbarevalp;
    fmpz_mod_mpolyn_struct * Aevalm, * Bevalm, * Gevalm, * Abarevalm, * Bbarevalm;
    fmpz_mod_mpolyn_struct * T1, * T2;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fmpz_mod_poly_struct * cA, * cB, * cG, * cAbar, * cBbar, * gamma;
    fmpz_mod_poly_struct * modulus, * alphapow, * t1;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
#if FLINT_WANT_ASSERT
    fmpz_mod_mpolyn_t Aorg, Borg;
    fmpz_mod_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(bits == St->mpolyn_stack->bits);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);

    FLINT_ASSERT(var > 0);

    if (var == 1)
    {
        fmpz_mod_polyun_struct * mA, * mB, * mG, * mAbar, * mBbar;

        fmpz_mod_polyun_stack_fit_request(St->polyun_stack, 5);
        mA = fmpz_mod_polyun_stack_take_top(St->polyun_stack);
        mB = fmpz_mod_polyun_stack_take_top(St->polyun_stack);
        mG = fmpz_mod_polyun_stack_take_top(St->polyun_stack);
        mAbar = fmpz_mod_polyun_stack_take_top(St->polyun_stack);
        mBbar = fmpz_mod_polyun_stack_take_top(St->polyun_stack);

        fmpz_mod_mpolyn_get_polyun_swap(mA, A, ctx);
        fmpz_mod_mpolyn_get_polyun_swap(mB, B, ctx);

        success = fmpz_mod_polyu1n_gcd_brown_smprime(mG, mAbar, mBbar, mA, mB,
                                ctx->ffinfo, St->poly_stack, St->polyun_stack);
        fmpz_mod_mpolyn_set_polyun_swap(G, mG, ctx);
        fmpz_mod_mpolyn_set_polyun_swap(Abar, mAbar, ctx);
        fmpz_mod_mpolyn_set_polyun_swap(Bbar, mBbar, ctx);

        fmpz_mod_polyun_stack_give_back(St->polyun_stack, 5);

        return success;
    }

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_init(leadA, ctx->ffinfo);
    fmpz_mod_poly_init(leadB, ctx->ffinfo);
    fmpz_mod_poly_set(leadA, A->coeffs + 0, ctx->ffinfo);
    fmpz_mod_poly_set(leadB, B->coeffs + 0, ctx->ffinfo);
    fmpz_mod_mpolyn_init(Aorg, bits, ctx);
    fmpz_mod_mpolyn_init(Borg, bits, ctx);
    fmpz_mod_mpolyn_set(Aorg, A, ctx);
    fmpz_mod_mpolyn_set(Borg, B, ctx);
#endif

    fmpz_init(alpha);
    fmpz_init(temp);
    fmpz_init(gammaevalp);
    fmpz_init(gammaevalm);

    fmpz_mod_poly_stack_fit_request(St->poly_stack, 9);
    cA          = fmpz_mod_poly_stack_take_top(St->poly_stack);
    cB          = fmpz_mod_poly_stack_take_top(St->poly_stack);
    cG          = fmpz_mod_poly_stack_take_top(St->poly_stack);
    cAbar       = fmpz_mod_poly_stack_take_top(St->poly_stack);
    cBbar       = fmpz_mod_poly_stack_take_top(St->poly_stack);
    gamma       = fmpz_mod_poly_stack_take_top(St->poly_stack);
    alphapow    = fmpz_mod_poly_stack_take_top(St->poly_stack);
    modulus     = fmpz_mod_poly_stack_take_top(St->poly_stack);
    t1          = fmpz_mod_poly_stack_take_top(St->poly_stack);

    fmpz_mod_mpolyn_stack_fit_request(St->mpolyn_stack, 12, ctx);
    Aevalp      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Bevalp      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Gevalp      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Abarevalp   = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Bbarevalp   = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Aevalm      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Bevalm      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Gevalm      = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Abarevalm   = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    Bbarevalm   = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    T1          = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);
    T2          = fmpz_mod_mpolyn_stack_take_top(St->mpolyn_stack);

    _fmpz_mod_poly_vec_remove_content(cA, A->coeffs, A->length, ctx->ffinfo);
    _fmpz_mod_poly_vec_remove_content(cB, B->coeffs, B->length, ctx->ffinfo);
    _fmpz_mod_poly_gcd_cofactors(cG, cAbar, cBbar, cA, cB, ctx->ffinfo, t1);
    fmpz_mod_poly_gcd(gamma, A->coeffs + 0, B->coeffs + 0, ctx->ffinfo);

    ldegA = fmpz_mod_mpolyn_lastdeg(A, ctx);
    ldegB = fmpz_mod_mpolyn_lastdeg(B, ctx);
    deggamma = fmpz_mod_poly_degree(gamma, ctx->ffinfo);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    upper_limit = mpoly_gcd_info_get_brown_upper_limit(I, var, bound);

    fmpz_mod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1), ctx->ffinfo);

    fmpz_mod_poly_one(modulus, ctx->ffinfo);

    fmpz_fdiv_q_2exp(alpha, fmpz_mod_ctx_modulus(ctx->ffinfo), 1);

choose_prime:

    fmpz_sub_ui(alpha, alpha, 1);
    if (fmpz_sgn(alpha) <= 0)
    {
        success = 0;
        goto cleanup;
    }

    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    fmpz_one(alphapow->coeffs + 0);
    fmpz_set(alphapow->coeffs + 1, alpha);

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_eval2_pow(gammaevalp, gammaevalm, gamma, alphapow, ctx->ffinfo);
    if (fmpz_is_zero(gammaevalp) || fmpz_is_zero(gammaevalm))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fmpz_mod_mpolyn_interp_reduce_2sm_mpolyn(Aevalp, Aevalm, A, var, alphapow, ctx);
    fmpz_mod_mpolyn_interp_reduce_2sm_mpolyn(Bevalp, Bevalm, B, var, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    success = fmpz_mod_mpolyn_gcd_brown_smprime(Gevalp, Abarevalp, Bbarevalp,
                                          Aevalp, Bevalp, var - 1, ctx, I, St);
    success = success &&
              fmpz_mod_mpolyn_gcd_brown_smprime(Gevalm, Abarevalm, Bbarevalm,
                                          Aevalm, Bevalm, var - 1, ctx, I, St);
    if (success == 0)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(Gevalp->length > 0);
    FLINT_ASSERT(Abarevalp->length > 0);
    FLINT_ASSERT(Bbarevalp->length > 0);
    FLINT_ASSERT(Gevalm->length > 0);
    FLINT_ASSERT(Abarevalm->length > 0);
    FLINT_ASSERT(Bbarevalm->length > 0);

    if (fmpz_mod_mpolyn_is_nonzero_fmpz(Gevalp, ctx) ||
        fmpz_mod_mpolyn_is_nonzero_fmpz(Gevalm, ctx))
    {
        fmpz_mod_mpolyn_one(G, ctx);
        fmpz_mod_mpolyn_swap(Abar, A, ctx);
        fmpz_mod_mpolyn_swap(Bbar, B, ctx);
        goto successful_put_content;
    }

    if (Gevalp->coeffs[0].length != Gevalm->coeffs[0].length ||
        !mpoly_monomial_equal(Gevalp->exps + N*0, Gevalm->exps + N*0, N))
    {
        goto choose_prime;
    }

    /* the Geval have matching degrees */
    if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) > 0)
    {
        slong k = fmpz_mod_poly_degree(Gevalp->coeffs + 0, ctx->ffinfo);
        int cmp = mpoly_monomial_cmp_nomask_extra(G->exps + N*0,
                                    Gevalp->exps + N*0, N, offset, k << shift);
        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            fmpz_mod_poly_one(modulus, ctx->ffinfo);
        }
    }

    /* update interpolants */
    fmpz_mod_mpolyn_scalar_mul_fmpz_mod(Gevalp, gammaevalp, ctx);
    fmpz_mod_mpolyn_scalar_mul_fmpz_mod(Gevalm, gammaevalm, ctx);
    if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(temp, modulus, alpha, ctx->ffinfo);
        fmpz_mod_mul(temp, temp, alpha, ctx->ffinfo);
        fmpz_mod_add(temp, temp, temp, ctx->ffinfo);
        fmpz_mod_poly_scalar_div_fmpz(modulus, modulus, temp, ctx->ffinfo);

        changed = fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(&ldegG, G, T1,
                                  Gevalp, Gevalm, var, modulus, alphapow, ctx);
        if (!changed && _fmpz_mod_poly_degree(modulus) < upper_limit)
        {
            _fmpz_mod_poly_vec_remove_content(t1, G->coeffs, G->length, ctx->ffinfo);
            if (fmpz_mod_mpolyn_divides(T1, A, G, ctx) &&
                fmpz_mod_mpolyn_divides(T2, B, G, ctx))
            {
                fmpz_mod_mpolyn_swap(Abar, T1, ctx);
                fmpz_mod_mpolyn_swap(Bbar, T2, ctx);
                goto successful_fix_lc;
            }
            _fmpz_mod_poly_vec_mul_poly(G->coeffs, G->length, t1, ctx->ffinfo);
        }

        changed = fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(&ldegAbar, Abar, T1,
                            Abarevalp, Abarevalm, var, modulus, alphapow, ctx);
        if (!changed && _fmpz_mod_poly_degree(modulus) < upper_limit)
        {
            _fmpz_mod_poly_vec_remove_content(t1, Abar->coeffs, Abar->length, ctx->ffinfo);
            if (fmpz_mod_mpolyn_divides(T1, A, Abar, ctx) &&
                fmpz_mod_mpolyn_divides(T2, B, T1, ctx))
            {
                fmpz_mod_mpolyn_swap(G, T1, ctx);
                fmpz_mod_mpolyn_swap(Bbar, T2, ctx);
                goto successful_fix_lc;
            }
            _fmpz_mod_poly_vec_mul_poly(Abar->coeffs, Abar->length, t1, ctx->ffinfo);
        }

        changed = fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(&ldegBbar, Bbar, T1,
                            Bbarevalp, Bbarevalm, var, modulus, alphapow, ctx);
        if (!changed && _fmpz_mod_poly_degree(modulus) < upper_limit)
        {
            _fmpz_mod_poly_vec_remove_content(t1, Bbar->coeffs, Bbar->length, ctx->ffinfo);
            if (fmpz_mod_mpolyn_divides(T1, B, Bbar, ctx) &&
                fmpz_mod_mpolyn_divides(T2, A, T1, ctx))
            {
                fmpz_mod_mpolyn_swap(G, T1, ctx);
                fmpz_mod_mpolyn_swap(Abar, T2, ctx);
                goto successful_fix_lc;
            }
            _fmpz_mod_poly_vec_mul_poly(Bbar->coeffs, Bbar->length, t1, ctx->ffinfo);
        }
    }
    else
    {
        fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(&ldegG, G,
                                              Gevalp, Gevalm, var, alpha, ctx);
        fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(&ldegAbar, Abar,
                                        Abarevalp, Abarevalm, var, alpha, ctx);
        fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(&ldegBbar, Bbar,
                                        Bbarevalp, Bbarevalm, var, alpha, ctx);
    }

    fmpz_mod_mul(temp, alpha, alpha, ctx->ffinfo);
    fmpz_mod_neg(temp, temp, ctx->ffinfo);
    fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 2, temp, ctx->ffinfo);

    if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (deggamma + ldegA == ldegG + ldegAbar &&
        deggamma + ldegB == ldegG + ldegBbar)
    {
        goto successful;
    }

    fmpz_mod_poly_one(modulus, ctx->ffinfo);
    goto choose_prime;

successful_fix_lc:
    
    fmpz_set(temp, fmpz_mod_poly_lead(G->coeffs + 0, ctx->ffinfo));
    if (fmpz_is_one(temp))
        goto successful_put_content;

    _fmpz_mod_poly_vec_mul_fmpz_mod(Abar->coeffs, Abar->length, temp, ctx->ffinfo);
    _fmpz_mod_poly_vec_mul_fmpz_mod(Bbar->coeffs, Bbar->length, temp, ctx->ffinfo);
    fmpz_mod_inv(temp, temp, ctx->ffinfo);
    _fmpz_mod_poly_vec_mul_fmpz_mod(G->coeffs, G->length, temp, ctx->ffinfo);

    goto successful_put_content;

successful:

    _fmpz_mod_poly_vec_remove_content(t1, G->coeffs, G->length, ctx->ffinfo);
    _fmpz_mod_poly_vec_divexact_poly(Abar->coeffs, Abar->length, G->coeffs + 0, ctx->ffinfo);
    _fmpz_mod_poly_vec_divexact_poly(Bbar->coeffs, Bbar->length, G->coeffs + 0, ctx->ffinfo);

successful_put_content:

    _fmpz_mod_poly_vec_mul_poly(G->coeffs, G->length, cG, ctx->ffinfo);
    _fmpz_mod_poly_vec_mul_poly(Abar->coeffs, Abar->length, cAbar, ctx->ffinfo);
    _fmpz_mod_poly_vec_mul_poly(Bbar->coeffs, Bbar->length, cBbar, ctx->ffinfo);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_mpolyn_leadcoeff(G)));
        fmpz_mod_poly_mul(modulus, G->coeffs + 0, Abar->coeffs + 0, ctx->ffinfo);
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadA, ctx->ffinfo));
        fmpz_mod_poly_mul(modulus, G->coeffs + 0, Bbar->coeffs + 0, ctx->ffinfo);
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadB, ctx->ffinfo));
    }
    fmpz_mod_poly_clear(leadA, ctx->ffinfo);
    fmpz_mod_poly_clear(leadB, ctx->ffinfo);
    fmpz_mod_mpolyn_clear(Aorg, ctx);
    fmpz_mod_mpolyn_clear(Borg, ctx);
#endif

    fmpz_clear(alpha);
    fmpz_clear(temp);
    fmpz_clear(gammaevalp);
    fmpz_clear(gammaevalm);

    fmpz_mod_poly_stack_give_back(St->poly_stack, 9);
    fmpz_mod_mpolyn_stack_give_back(St->mpolyn_stack, 12);

    FLINT_ASSERT(success || fmpz_bits(fmpz_mod_mpoly_ctx_modulus(ctx)) < 60);

    return success;
}

