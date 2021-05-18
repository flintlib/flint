/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"


int nmod_mpolyn_gcd_brown_smprime_bivar(
    nmod_mpolyn_t G,
    nmod_mpolyn_t Abar,
    nmod_mpolyn_t Bbar,
    nmod_mpolyn_t A,
    nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx,
    nmod_poly_stack_t Sp)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaevalp, gammaevalm;
    n_poly_struct * Aevalp, * Bevalp, * Gevalp, * Abarevalp, * Bbarevalp;
    n_poly_struct * Aevalm, * Bevalm, * Gevalm, * Abarevalm, * Bbarevalm;
    nmod_mpolyn_struct * T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    n_poly_struct * cA, * cB, * cG, * cAbar, * cBbar, * gamma;
    n_poly_struct * modulus, * modulus2, * alphapow, * r;
    int gstab, astab, bstab, use_stab;
    slong N, off, shift;
    flint_bitcnt_t bits = A->bits;
#if FLINT_WANT_ASSERT
    n_poly_t leadA, leadB;
    const slong Sp_size_poly = nmod_poly_stack_size_poly(Sp);
    const slong Sp_size_mpolyn = nmod_poly_stack_size_mpolyn(Sp);
#endif

    FLINT_ASSERT(Sp->ctx->mod.n == ctx->mod.n);
    FLINT_ASSERT(Sp->ctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(Sp->bits == bits);
    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B->bits == bits);
    FLINT_ASSERT(G->bits == bits);
    FLINT_ASSERT(Abar->bits == bits);
    FLINT_ASSERT(Bbar->bits == bits);

#if FLINT_WANT_ASSERT
    n_poly_init(leadA);
    n_poly_init(leadB);
    n_poly_set(leadA, nmod_mpolyn_leadcoeff_poly(A, ctx));
    n_poly_set(leadB, nmod_mpolyn_leadcoeff_poly(B, ctx));
#endif

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, bits, ctx->minfo);

    nmod_poly_stack_fit_request_poly(Sp, 20);
    cA          = nmod_poly_stack_take_top_poly(Sp);
    cB          = nmod_poly_stack_take_top_poly(Sp);
    cG          = nmod_poly_stack_take_top_poly(Sp);
    cAbar       = nmod_poly_stack_take_top_poly(Sp);
    cBbar       = nmod_poly_stack_take_top_poly(Sp);
    gamma       = nmod_poly_stack_take_top_poly(Sp);
    Aevalp      = nmod_poly_stack_take_top_poly(Sp);
    Bevalp      = nmod_poly_stack_take_top_poly(Sp);
    Gevalp      = nmod_poly_stack_take_top_poly(Sp);
    Abarevalp   = nmod_poly_stack_take_top_poly(Sp);
    Bbarevalp   = nmod_poly_stack_take_top_poly(Sp);
    Aevalm      = nmod_poly_stack_take_top_poly(Sp);
    Bevalm      = nmod_poly_stack_take_top_poly(Sp);
    Gevalm      = nmod_poly_stack_take_top_poly(Sp);
    Abarevalm   = nmod_poly_stack_take_top_poly(Sp);
    Bbarevalm   = nmod_poly_stack_take_top_poly(Sp);
    r           = nmod_poly_stack_take_top_poly(Sp);
    alphapow    = nmod_poly_stack_take_top_poly(Sp);
    modulus     = nmod_poly_stack_take_top_poly(Sp);
    modulus2    = nmod_poly_stack_take_top_poly(Sp);

    nmod_poly_stack_fit_request_mpolyn(Sp, 1);
    T           = nmod_poly_stack_take_top_mpolyn(Sp);

    nmod_mpolyn_content_last(cA, A, ctx);
    nmod_mpolyn_content_last(cB, B, ctx);
    nmod_mpolyn_divexact_last(A, cA, ctx);
    nmod_mpolyn_divexact_last(B, cB, ctx);

    n_poly_mod_gcd(cG, cA, cB, ctx->mod);

    n_poly_mod_div(cAbar, cA, cG, ctx->mod);
    n_poly_mod_div(cBbar, cB, cG, ctx->mod);

    n_poly_mod_gcd(gamma, nmod_mpolyn_leadcoeff_poly(A, ctx),
                          nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->mod);

    ldegA = nmod_mpolyn_lastdeg(A, ctx);
    ldegB = nmod_mpolyn_lastdeg(B, ctx);
    deggamma = n_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    n_poly_one(modulus);

    if ((ctx->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    use_stab = 1;
    gstab = bstab = astab = 0;

    alpha = (ctx->mod.n - UWORD(1))/UWORD(2);

choose_prime: /* primes are v - alpha, v + alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm, gamma, alphapow, ctx->mod);
    if (gammaevalp == 0 || gammaevalm == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    nmod_mpolyn_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, ctx);
    nmod_mpolyn_interp_reduce_2sm_poly(Bevalp, Bevalm, B, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    if (use_stab && gstab)
    {
        slong Gdeg;
        nmod_mpolyn_interp_reduce_2sm_poly(Gevalp, Gevalm, G, alphapow, ctx);
        Gdeg = ((G->exps + N*0)[off]>>shift);
        success = 1;
        success = success && n_poly_degree(Gevalp) == Gdeg;
        success = success && n_poly_degree(Gevalm) == Gdeg;
        success = success && Gevalp->coeffs[Gdeg] == gammaevalp;
        success = success && Gevalm->coeffs[Gdeg] == gammaevalm;
        n_poly_mod_divrem(Abarevalp, r, Aevalp, Gevalp, ctx->mod);
        success = success && (r->length == 0);
        n_poly_mod_divrem(Abarevalm, r, Aevalm, Gevalm, ctx->mod);
        success = success && (r->length == 0);
        n_poly_mod_divrem(Bbarevalp, r, Bevalp, Gevalp, ctx->mod);
        success = success && (r->length == 0);
        n_poly_mod_divrem(Bbarevalm, r, Bevalm, Gevalm, ctx->mod);
        success = success && (r->length == 0);

        if (!success)
        {
            use_stab = 0;
            n_poly_one(modulus);
            alpha = (ctx->mod.n - UWORD(1))/UWORD(2);
            goto choose_prime;
        }

        _n_poly_mod_scalar_mul_nmod(Abarevalp, Abarevalp, gammaevalp, ctx->mod);
        _n_poly_mod_scalar_mul_nmod(Abarevalm, Abarevalm, gammaevalm, ctx->mod);
        _n_poly_mod_scalar_mul_nmod(Bbarevalp, Bbarevalp, gammaevalp, ctx->mod);
        _n_poly_mod_scalar_mul_nmod(Bbarevalm, Bbarevalm, gammaevalm, ctx->mod);
    }
    else
    {
        n_poly_mod_gcd(Gevalp, Aevalp, Bevalp, ctx->mod);
        n_poly_mod_div(Abarevalp, Aevalp, Gevalp, ctx->mod);
        n_poly_mod_div(Bbarevalp, Bevalp, Gevalp, ctx->mod);
        n_poly_mod_gcd(Gevalm, Aevalm, Bevalm, ctx->mod);
        n_poly_mod_div(Abarevalm, Aevalm, Gevalm, ctx->mod);
        n_poly_mod_div(Bbarevalm, Bevalm, Gevalm, ctx->mod);
        gstab = astab = bstab = 0;
    }

    FLINT_ASSERT(Gevalp->length > 0);
    FLINT_ASSERT(Abarevalp->length > 0);
    FLINT_ASSERT(Bbarevalp->length > 0);
    FLINT_ASSERT(Gevalm->length > 0);
    FLINT_ASSERT(Abarevalm->length > 0);
    FLINT_ASSERT(Bbarevalm->length > 0);

    if (n_poly_degree(Gevalp) == 0 || n_poly_degree(Gevalm) == 0)
    {
        nmod_mpolyn_one(G, ctx);
        nmod_mpolyn_swap(Abar, A);
        nmod_mpolyn_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (n_poly_degree(Gevalp) != n_poly_degree(Gevalm))
    {
        goto choose_prime;
    }

    /* the Geval have matching degrees */
    if (n_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (n_poly_degree(Gevalp) > ((G->exps + N*0)[off]>>shift))
        {
            goto choose_prime;
        }
        else if (n_poly_degree(Gevalp) < ((G->exps + N*0)[off]>>shift))
        {
            n_poly_one(modulus);
        }
    }

    /* update interpolants */
    _n_poly_mod_scalar_mul_nmod(Gevalp, Gevalp, gammaevalp, ctx->mod);
    _n_poly_mod_scalar_mul_nmod(Gevalm, Gevalm, gammaevalm, ctx->mod);
    if (n_poly_degree(modulus) > 0)
    {
        temp = n_poly_mod_evaluate_nmod(modulus, alpha, ctx->mod);
        FLINT_ASSERT(temp == n_poly_mod_evaluate_nmod(modulus, ctx->mod.n - alpha, ctx->mod));
        temp = nmod_mul(temp, alpha, ctx->mod);
        temp = nmod_add(temp, temp, ctx->mod);
        temp = n_invmod(temp, ctx->mod.n);
        _n_poly_mod_scalar_mul_nmod(modulus, modulus, temp, ctx->mod);
        if (!gstab)
        {
            gstab = !nmod_mpolyn_interp_crt_2sm_poly(&ldegG, G, T,
                                       Gevalp, Gevalm, modulus, alphapow, ctx);
        }
        nmod_mpolyn_interp_crt_2sm_poly(&ldegAbar, Abar, T,
                                 Abarevalp, Abarevalm, modulus, alphapow, ctx);
        nmod_mpolyn_interp_crt_2sm_poly(&ldegBbar, Bbar, T,
                                 Bbarevalp, Bbarevalm, modulus, alphapow, ctx);
    }
    else
    {
        nmod_mpolyn_interp_lift_2sm_poly(&ldegG, G, Gevalp, Gevalm, alpha, ctx);
        nmod_mpolyn_interp_lift_2sm_poly(&ldegAbar, Abar,
                                             Abarevalp, Abarevalm, alpha, ctx);
        nmod_mpolyn_interp_lift_2sm_poly(&ldegBbar, Bbar,
                                             Bbarevalp, Bbarevalm, alpha, ctx);
        gstab = astab = bstab = 0;
    }
    temp = nmod_mul(alpha, alpha, ctx->mod);
    _n_poly_mod_scalar_mul_nmod(modulus2, modulus, temp, ctx->mod);
    n_poly_shift_left(modulus, modulus, 2);
    n_poly_mod_sub(modulus, modulus, modulus2, ctx->mod);

    if (n_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    n_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyn_content_last(modulus, G, ctx);
    nmod_mpolyn_divexact_last(G, modulus, ctx);
    nmod_mpolyn_divexact_last(Abar, nmod_mpolyn_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyn_divexact_last(Bbar, nmod_mpolyn_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyn_mul_last(G, cG, ctx);
    nmod_mpolyn_mul_last(Abar, cAbar, ctx);
    nmod_mpolyn_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyn_leadcoeff(G, ctx));
        n_poly_mod_mul(modulus, nmod_mpolyn_leadcoeff_poly(G, ctx),
                               nmod_mpolyn_leadcoeff_poly(Abar, ctx), ctx->mod);
        FLINT_ASSERT(n_poly_equal(modulus, leadA));
        n_poly_mod_mul(modulus, nmod_mpolyn_leadcoeff_poly(G, ctx),
                                nmod_mpolyn_leadcoeff_poly(Bbar, ctx), ctx->mod);
        FLINT_ASSERT(n_poly_equal(modulus, leadB));
    }
    n_poly_clear(leadA);
    n_poly_clear(leadB);
#endif

    nmod_poly_stack_give_back_poly(Sp, 20);
    nmod_poly_stack_give_back_mpolyn(Sp, 1);
    FLINT_ASSERT(Sp_size_poly == nmod_poly_stack_size_poly(Sp));
    FLINT_ASSERT(Sp_size_mpolyn == nmod_poly_stack_size_mpolyn(Sp));

    return success;
}


int nmod_mpolyn_gcd_brown_smprime(
    nmod_mpolyn_t G,
    nmod_mpolyn_t Abar,
    nmod_mpolyn_t Bbar,
    nmod_mpolyn_t A,
    nmod_mpolyn_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx,
    const mpoly_gcd_info_t I,
    nmod_poly_stack_t Sp)
{
    int success;
    int changed;
    slong bound;
    slong upper_limit;
    slong offset, shift;
    mp_limb_t alpha, temp, gammaevalp, gammaevalm;
    nmod_mpolyn_struct * Aevalp, * Bevalp, * Gevalp, * Abarevalp, * Bbarevalp;
    nmod_mpolyn_struct * Aevalm, * Bevalm, * Gevalm, * Abarevalm, * Bbarevalm;
    nmod_mpolyn_struct * T1, * T2;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    n_poly_struct * cA, * cB, * cG, * cAbar, * cBbar, * gamma;
    n_poly_struct * modulus, * modulus2, * alphapow;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
#if FLINT_WANT_ASSERT
    nmod_mpolyn_t Aorg, Borg;
    n_poly_t leadA, leadB;
    slong Sp_size_poly = nmod_poly_stack_size_poly(Sp);
    slong Sp_size_mpolyn = nmod_poly_stack_size_mpolyn(Sp);
#endif

    FLINT_ASSERT(Sp->ctx->mod.n == ctx->mod.n);
    FLINT_ASSERT(Sp->ctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(Sp->bits == A->bits);
    FLINT_ASSERT(Sp->bits == B->bits);
    FLINT_ASSERT(Sp->bits == G->bits);
    FLINT_ASSERT(Sp->bits == Abar->bits);
    FLINT_ASSERT(Sp->bits == Bbar->bits);

    FLINT_ASSERT(var > 0);

    if (var == 1)
        return nmod_mpolyn_gcd_brown_smprime_bivar(G, Abar, Bbar, A, B, ctx, Sp);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if FLINT_WANT_ASSERT
    n_poly_init(leadA);
    n_poly_init(leadB);
    n_poly_set(leadA, nmod_mpolyn_leadcoeff_poly(A, ctx));
    n_poly_set(leadB, nmod_mpolyn_leadcoeff_poly(B, ctx));
    nmod_mpolyn_init(Aorg, A->bits, ctx);
    nmod_mpolyn_init(Borg, B->bits, ctx);
    nmod_mpolyn_set(Aorg, A, ctx);
    nmod_mpolyn_set(Borg, B, ctx);
#endif

    nmod_poly_stack_fit_request_poly(Sp, 9);
    cA          = nmod_poly_stack_take_top_poly(Sp);
    cB          = nmod_poly_stack_take_top_poly(Sp);
    cG          = nmod_poly_stack_take_top_poly(Sp);
    cAbar       = nmod_poly_stack_take_top_poly(Sp);
    cBbar       = nmod_poly_stack_take_top_poly(Sp);
    gamma       = nmod_poly_stack_take_top_poly(Sp);
    alphapow    = nmod_poly_stack_take_top_poly(Sp);
    modulus     = nmod_poly_stack_take_top_poly(Sp);
    modulus2    = nmod_poly_stack_take_top_poly(Sp);

    nmod_poly_stack_fit_request_mpolyn(Sp, 12);
    Aevalp      = nmod_poly_stack_take_top_mpolyn(Sp);
    Bevalp      = nmod_poly_stack_take_top_mpolyn(Sp);
    Gevalp      = nmod_poly_stack_take_top_mpolyn(Sp);
    Abarevalp   = nmod_poly_stack_take_top_mpolyn(Sp);
    Bbarevalp   = nmod_poly_stack_take_top_mpolyn(Sp);
    Aevalm      = nmod_poly_stack_take_top_mpolyn(Sp);
    Bevalm      = nmod_poly_stack_take_top_mpolyn(Sp);
    Gevalm      = nmod_poly_stack_take_top_mpolyn(Sp);
    Abarevalm   = nmod_poly_stack_take_top_mpolyn(Sp);
    Bbarevalm   = nmod_poly_stack_take_top_mpolyn(Sp);
    T1          = nmod_poly_stack_take_top_mpolyn(Sp);
    T2          = nmod_poly_stack_take_top_mpolyn(Sp);

    nmod_mpolyn_content_last(cA, A, ctx);
    nmod_mpolyn_content_last(cB, B, ctx);
    nmod_mpolyn_divexact_last(A, cA, ctx);
    nmod_mpolyn_divexact_last(B, cB, ctx);

    n_poly_mod_gcd(cG, cA, cB, ctx->mod);

    n_poly_mod_div(cAbar, cA, cG, ctx->mod);
    n_poly_mod_div(cBbar, cB, cG, ctx->mod);

    n_poly_mod_gcd(gamma, nmod_mpolyn_leadcoeff_poly(A, ctx),
                          nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->mod);

    ldegA = nmod_mpolyn_lastdeg(A, ctx);
    ldegB = nmod_mpolyn_lastdeg(B, ctx);
    deggamma = n_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    upper_limit = mpoly_gcd_info_get_brown_upper_limit(I, var, bound);

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));

    n_poly_one(modulus);

    if ((ctx->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx->mod.n - UWORD(1))/UWORD(2);

choose_prime:

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx->mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    _nmod_poly_eval2_pow(&gammaevalp, &gammaevalm, gamma, alphapow, ctx->mod);
    if (gammaevalp == 0 || gammaevalm == 0)
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    nmod_mpolyn_interp_reduce_2sm_mpolyn(Aevalp, Aevalm, A, var, alphapow, ctx);
    nmod_mpolyn_interp_reduce_2sm_mpolyn(Bevalp, Bevalm, B, var, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    success = nmod_mpolyn_gcd_brown_smprime(Gevalp, Abarevalp, Bbarevalp,
                                          Aevalp, Bevalp, var - 1, ctx, I, Sp);
    success = success
           && nmod_mpolyn_gcd_brown_smprime(Gevalm, Abarevalm, Bbarevalm,
                                          Aevalm, Bevalm, var - 1, ctx, I, Sp);
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

    if (   nmod_mpolyn_is_nonzero_nmod(Gevalp, ctx)
        || nmod_mpolyn_is_nonzero_nmod(Gevalm, ctx))
    {
        nmod_mpolyn_one(G, ctx);
        nmod_mpolyn_swap(Abar, A);
        nmod_mpolyn_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (n_poly_degree(Gevalp->coeffs + 0) !=
        n_poly_degree(Gevalm->coeffs + 0))
    {
        goto choose_prime;
    }
    if (!mpoly_monomial_equal(Gevalp->exps + N*0, Gevalm->exps + N*0, N))
    {
        goto choose_prime;
    }

    /* the Geval have matching degrees */
    if (n_poly_degree(modulus) > 0)
    {
        int cmp;
        slong k;

        FLINT_ASSERT(G->length > 0);

        k = n_poly_degree(Gevalp->coeffs + 0);
        cmp = mpoly_monomial_cmp_nomask_extra(G->exps + N*0,
                                    Gevalp->exps + N*0, N, offset, k << shift);
        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            n_poly_one(modulus);
        }
    }

    /* update interpolants */
    temp = nmod_mpolyn_leadcoeff(Gevalp, ctx);
    temp = n_invmod(temp, ctx->mod.n);
    temp = nmod_mul(gammaevalp, temp, ctx->mod);
    nmod_mpolyn_scalar_mul_nmod(Gevalp, temp, ctx);
    temp = nmod_mpolyn_leadcoeff(Gevalm, ctx);
    temp = n_invmod(temp, ctx->mod.n);
    temp = nmod_mul(gammaevalm, temp, ctx->mod);
    nmod_mpolyn_scalar_mul_nmod(Gevalm, temp, ctx);
    if (n_poly_degree(modulus) > 0)
    {
        temp = n_poly_mod_evaluate_nmod(modulus, alpha, ctx->mod);
        FLINT_ASSERT(temp == n_poly_mod_evaluate_nmod(modulus, ctx->mod.n - alpha, ctx->mod));
        temp = nmod_mul(temp, alpha, ctx->mod);
        temp = nmod_add(temp, temp, ctx->mod);
        temp = n_invmod(temp, ctx->mod.n);
        _n_poly_mod_scalar_mul_nmod_inplace(modulus, temp, ctx->mod);

        changed = nmod_mpolyn_interp_crt_2sm_mpolyn(&ldegG, G, T1,
                                  Gevalp, Gevalm, var, modulus, alphapow, ctx);
        if (!changed && n_poly_degree(modulus) < upper_limit)
        {
            nmod_mpolyn_content_last(modulus2, G, ctx);
            nmod_mpolyn_divexact_last(G, modulus2, ctx);
            success =            nmod_mpolyn_divides(T1, A, G, ctx);
            success = success && nmod_mpolyn_divides(T2, B, G, ctx);
            if (success)
            {
                nmod_mpolyn_swap(T1, Abar);
                nmod_mpolyn_swap(T2, Bbar);
successful_fix_lc:
                temp = nmod_mpolyn_leadcoeff(G, ctx);
                nmod_mpolyn_scalar_mul_nmod(Abar, temp, ctx);
                nmod_mpolyn_scalar_mul_nmod(Bbar, temp, ctx);
                temp = n_invmod(temp, ctx->mod.n);
                nmod_mpolyn_scalar_mul_nmod(G, temp, ctx);
                goto successful_put_content;
            }
            else
            {
                nmod_mpolyn_mul_last(G, modulus2, ctx); /* restore G */
            }
        }

        changed = nmod_mpolyn_interp_crt_2sm_mpolyn(&ldegAbar, Abar, T1,
                            Abarevalp, Abarevalm, var, modulus, alphapow, ctx);
        if (!changed && n_poly_degree(modulus) < upper_limit)
        {
            nmod_mpolyn_content_last(modulus2, Abar, ctx);
            nmod_mpolyn_divexact_last(Abar, modulus2, ctx);
            success =            nmod_mpolyn_divides(T1, A, Abar, ctx);
            success = success && nmod_mpolyn_divides(T2, B, T1, ctx);
            if (success)
            {
                nmod_mpolyn_swap(T1, G);
                nmod_mpolyn_swap(T2, Bbar);
                goto successful_fix_lc;
            }
            else
            {
                nmod_mpolyn_mul_last(Abar, modulus2, ctx); /* restore Abar */
            }
        }

        changed = nmod_mpolyn_interp_crt_2sm_mpolyn(&ldegBbar, Bbar, T1,
                            Bbarevalp, Bbarevalm, var, modulus, alphapow, ctx);
        if (!changed && n_poly_degree(modulus) < upper_limit)
        {
            nmod_mpolyn_content_last(modulus2, Bbar, ctx);
            nmod_mpolyn_divexact_last(Bbar, modulus2, ctx);
            success =            nmod_mpolyn_divides(T1, B, Bbar, ctx);
            success = success && nmod_mpolyn_divides(T2, A, T1, ctx);
            if (success)
            {
                nmod_mpolyn_swap(T1, G);
                nmod_mpolyn_swap(T2, Abar);
                goto successful_fix_lc;
            }
            else
            {
                nmod_mpolyn_mul_last(Bbar, modulus2, ctx); /* restore Bbar */
            }
        }
    }
    else
    {
        nmod_mpolyn_interp_lift_2sm_mpolyn(&ldegG, G,
                                              Gevalp, Gevalm, var, alpha, ctx);
        nmod_mpolyn_interp_lift_2sm_mpolyn(&ldegAbar, Abar,
                                        Abarevalp, Abarevalm, var, alpha, ctx);
        nmod_mpolyn_interp_lift_2sm_mpolyn(&ldegBbar, Bbar,
                                        Bbarevalp, Bbarevalm, var, alpha, ctx);
    }

    temp = nmod_mul(alpha, alpha, ctx->mod);
    _n_poly_mod_scalar_mul_nmod(modulus2, modulus, temp, ctx->mod);
    n_poly_shift_left(modulus, modulus, 2);
    n_poly_mod_sub(modulus, modulus, modulus2, ctx->mod);

    if (n_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    n_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyn_content_last(modulus2, G, ctx);
    nmod_mpolyn_divexact_last(G, modulus2, ctx);
    nmod_mpolyn_divexact_last(Abar, nmod_mpolyn_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyn_divexact_last(Bbar, nmod_mpolyn_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyn_mul_last(G, cG, ctx);
    nmod_mpolyn_mul_last(Abar, cAbar, ctx);
    nmod_mpolyn_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyn_leadcoeff(G, ctx));
        n_poly_mod_mul(modulus, nmod_mpolyn_leadcoeff_poly(G, ctx),
                                nmod_mpolyn_leadcoeff_poly(Abar, ctx), ctx->mod);
        FLINT_ASSERT(n_poly_equal(modulus, leadA));
        n_poly_mod_mul(modulus, nmod_mpolyn_leadcoeff_poly(G, ctx),
                               nmod_mpolyn_leadcoeff_poly(Bbar, ctx), ctx->mod);
        FLINT_ASSERT(n_poly_equal(modulus, leadB));

        success =            nmod_mpolyn_divides(T1, Aorg, G, ctx);
        success = success && nmod_mpolyn_divides(T2, Borg, G, ctx);
        FLINT_ASSERT(success);
        FLINT_ASSERT(nmod_mpolyn_equal(T1, Abar, ctx));
        FLINT_ASSERT(nmod_mpolyn_equal(T2, Bbar, ctx));

        success = 1;
    }
    n_poly_clear(leadA);
    n_poly_clear(leadB);
    nmod_mpolyn_clear(Aorg, ctx);
    nmod_mpolyn_clear(Borg, ctx);
#endif

    nmod_poly_stack_give_back_poly(Sp, 9);
    nmod_poly_stack_give_back_mpolyn(Sp, 12);
    FLINT_ASSERT(Sp_size_poly == nmod_poly_stack_size_poly(Sp));
    FLINT_ASSERT(Sp_size_mpolyn == nmod_poly_stack_size_mpolyn(Sp));

    return success;
}



int nmod_mpolyn_gcd_brown_lgprime_bivar(
    nmod_mpolyn_t G,
    nmod_mpolyn_t Abar,
    nmod_mpolyn_t Bbar,
    nmod_mpolyn_t A,
    nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t temp, gammaeval;
    fq_nmod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyn_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    n_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    n_poly_t modulus;
    slong deg;
    fq_nmod_mpoly_ctx_t ectx;
    slong N, off, shift;
#if FLINT_WANT_ASSERT
    n_poly_t leadA, leadB;
#endif

#if FLINT_WANT_ASSERT
    n_poly_init(leadA);
    n_poly_init(leadB);
    n_poly_set(leadA, nmod_mpolyn_leadcoeff_poly(A, ctx));
    n_poly_set(leadB, nmod_mpolyn_leadcoeff_poly(B, ctx));
#endif

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    n_poly_init(cA);
    n_poly_init(cB);
    nmod_mpolyn_content_last(cA, A, ctx);
    nmod_mpolyn_content_last(cB, B, ctx);
    nmod_mpolyn_divexact_last(A, cA, ctx);
    nmod_mpolyn_divexact_last(B, cB, ctx);

    n_poly_init(cG);
    n_poly_mod_gcd(cG, cA, cB, ctx->mod);

    n_poly_init(cAbar);
    n_poly_init(cBbar);
    n_poly_mod_div(cAbar, cA, cG, ctx->mod);
    n_poly_mod_div(cBbar, cB, cG, ctx->mod);

    n_poly_init(gamma);
    n_poly_mod_gcd(gamma, nmod_mpolyn_leadcoeff_poly(A, ctx),
                          nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->mod);

    ldegA = nmod_mpolyn_lastdeg(A, ctx);
    ldegB = nmod_mpolyn_lastdeg(B, ctx);
    deggamma = n_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_mpolyn_init(T, A->bits, ctx);

    n_poly_init(modulus);
    n_poly_one(modulus);

    deg = WORD(20)/(FLINT_BIT_COUNT(ctx->mod.n));
    deg = FLINT_MAX(WORD(2), deg);
    fq_nmod_mpoly_ctx_init_deg(ectx, ctx->minfo->nvars, ORD_LEX,
                                                      ctx->mod.n, deg);

    fq_nmod_poly_init(Aeval, ectx->fqctx);
    fq_nmod_poly_init(Beval, ectx->fqctx);
    fq_nmod_poly_init(Geval, ectx->fqctx);
    fq_nmod_poly_init(Abareval, ectx->fqctx);
    fq_nmod_poly_init(Bbareval, ectx->fqctx);
    fq_nmod_init(gammaeval, ectx->fqctx);
    fq_nmod_init(temp, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime: /* prime will be irreducible element of Fp[v] */

    /* same TODO */
    deg++;
    if (deg > 10000)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

    fq_nmod_mpoly_ctx_change_modulus(ectx, deg);

have_prime:

    /* make sure reduction does not kill both lc(A) and lc(B) */
    n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(gammaeval), gamma,
          evil_const_cast_nmod_poly_to_n_poly(ectx->fqctx->modulus), ctx->mod);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
    {
        goto choose_prime;
    }

    /* reduction should kill neither A nor B */
    nmod_mpolyn_interp_reduce_lg_poly(Aeval, ectx->fqctx, A, ctx);
    nmod_mpolyn_interp_reduce_lg_poly(Beval, ectx->fqctx, B, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    fq_nmod_poly_gcd(Geval, Aeval, Beval, ectx->fqctx);
    success = fq_nmod_poly_divides(Abareval, Aeval, Geval, ectx->fqctx);
    FLINT_ASSERT(success);
    success = fq_nmod_poly_divides(Bbareval, Beval, Geval, ectx->fqctx);
    FLINT_ASSERT(success);

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fq_nmod_poly_degree(Geval, ectx->fqctx) == 0)
    {
        nmod_mpolyn_one(G, ctx);
        nmod_mpolyn_swap(Abar, A);
        nmod_mpolyn_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (n_poly_degree(modulus) > 0)
    {
        slong Gdeg;
        FLINT_ASSERT(G->length > 0);
        Gdeg = (G->exps + N*0)[off]>>shift;
        if (fq_nmod_poly_degree(Geval, ectx->fqctx) > Gdeg)
        {
            goto choose_prime;
        }
        else if (fq_nmod_poly_degree(Geval, ectx->fqctx) < Gdeg)
        {
            n_poly_one(modulus);
        }
    }

    FLINT_ASSERT(fq_nmod_is_one(Geval->coeffs + Geval->length - 1, ectx->fqctx));
    fq_nmod_poly_scalar_mul_fq_nmod(Geval, Geval, gammaeval, ectx->fqctx);

    if (n_poly_degree(modulus) > 0)
    {
        nmod_mpolyn_interp_crt_lg_poly(&ldegG, G, T, modulus, ctx,
                                                           Geval, ectx->fqctx);
        nmod_mpolyn_interp_crt_lg_poly(&ldegAbar, Abar, T, modulus, ctx,
                                                        Abareval, ectx->fqctx);
        nmod_mpolyn_interp_crt_lg_poly(&ldegBbar, Bbar, T, modulus, ctx,
                                                        Bbareval, ectx->fqctx);
    }
    else
    {
        nmod_mpolyn_interp_lift_lg_poly(&ldegG, G, ctx, Geval, ectx->fqctx);
        nmod_mpolyn_interp_lift_lg_poly(&ldegAbar, Abar, ctx,
                                                        Abareval, ectx->fqctx);
        nmod_mpolyn_interp_lift_lg_poly(&ldegBbar, Bbar, ctx,
                                                        Bbareval, ectx->fqctx);
    }

    n_poly_mod_mul(modulus, modulus,
          evil_const_cast_nmod_poly_to_n_poly(ectx->fqctx->modulus), ctx->mod);

    if (n_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    n_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyn_content_last(modulus, G, ctx);
    nmod_mpolyn_divexact_last(G, modulus, ctx);
    nmod_mpolyn_divexact_last(Abar, nmod_mpolyn_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyn_divexact_last(Bbar, nmod_mpolyn_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyn_mul_last(G, cG, ctx);
    nmod_mpolyn_mul_last(Abar, cAbar, ctx);
    nmod_mpolyn_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyn_leadcoeff(G, ctx));
        n_poly_mod_mul(modulus, nmod_mpolyn_leadcoeff_poly(G, ctx),
                              nmod_mpolyn_leadcoeff_poly(Abar, ctx), ctx->mod);
        FLINT_ASSERT(n_poly_equal(modulus, leadA));
        n_poly_mod_mul(modulus, nmod_mpolyn_leadcoeff_poly(G, ctx),
                              nmod_mpolyn_leadcoeff_poly(Bbar, ctx), ctx->mod);
        FLINT_ASSERT(n_poly_equal(modulus, leadB));
    }
    n_poly_clear(leadA);
    n_poly_clear(leadB);
#endif

    n_poly_clear(cA);
    n_poly_clear(cB);
    n_poly_clear(cG);
    n_poly_clear(cAbar);
    n_poly_clear(cBbar);
    n_poly_clear(gamma);
    n_poly_clear(modulus);

    nmod_mpolyn_clear(T, ctx);

    fq_nmod_poly_clear(Aeval, ectx->fqctx);
    fq_nmod_poly_clear(Beval, ectx->fqctx);
    fq_nmod_poly_clear(Geval, ectx->fqctx);
    fq_nmod_poly_clear(Abareval, ectx->fqctx);
    fq_nmod_poly_clear(Bbareval, ectx->fqctx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(temp, ectx->fqctx);

    fq_nmod_mpoly_ctx_clear(ectx);

    return success;
}

int nmod_mpolyn_gcd_brown_lgprime(
    nmod_mpolyn_t G,
    nmod_mpolyn_t Abar,
    nmod_mpolyn_t Bbar,
    nmod_mpolyn_t A,
    nmod_mpolyn_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    slong offset, shift;
    fq_nmod_t temp, gammaeval;
    fq_nmod_mpolyn_t Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyn_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    n_poly_t cA, cB, cG, cAbar, cBbar, gamma;
    n_poly_t modulus;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong deg;
    fq_nmod_mpoly_ctx_t ectx;
#if FLINT_WANT_ASSERT
    n_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(var > 0);

    if (var == 1)
        return nmod_mpolyn_gcd_brown_lgprime_bivar(G, Abar, Bbar, A, B, ctx);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if FLINT_WANT_ASSERT
    n_poly_init(leadA);
    n_poly_init(leadB);
    n_poly_set(leadA, nmod_mpolyn_leadcoeff_poly(A, ctx));
    n_poly_set(leadB, nmod_mpolyn_leadcoeff_poly(B, ctx));
#endif

    n_poly_init(cA);
    n_poly_init(cB);
    nmod_mpolyn_content_last(cA, A, ctx);
    nmod_mpolyn_content_last(cB, B, ctx);
    nmod_mpolyn_divexact_last(A, cA, ctx);
    nmod_mpolyn_divexact_last(B, cB, ctx);

    n_poly_init(cG);
    n_poly_mod_gcd(cG, cA, cB, ctx->mod);

    n_poly_init(cAbar);
    n_poly_init(cBbar);
    n_poly_mod_div(cAbar, cA, cG, ctx->mod);
    n_poly_mod_div(cBbar, cB, cG, ctx->mod);

    n_poly_init(gamma);
    n_poly_mod_gcd(gamma, nmod_mpolyn_leadcoeff_poly(A, ctx),
                          nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->mod);

    ldegA = nmod_mpolyn_lastdeg(A, ctx);
    ldegB = nmod_mpolyn_lastdeg(B, ctx);
    deggamma = n_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    nmod_mpolyn_init(T, bits, ctx);

    n_poly_init(modulus);
    n_poly_one(modulus);

    deg = WORD(20)/(FLINT_BIT_COUNT(ctx->mod.n));
    deg = FLINT_MAX(WORD(2), deg);

    fq_nmod_mpoly_ctx_init_deg(ectx, ctx->minfo->nvars, ORD_LEX,
                                                      ctx->mod.n, deg);

    fq_nmod_mpolyn_init(Aeval, bits, ectx);
    fq_nmod_mpolyn_init(Beval, bits, ectx);
    fq_nmod_mpolyn_init(Geval, bits, ectx);
    fq_nmod_mpolyn_init(Abareval, bits, ectx);
    fq_nmod_mpolyn_init(Bbareval, bits, ectx);
    fq_nmod_init(gammaeval, ectx->fqctx);
    fq_nmod_init(temp, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime: /* prime will be irreducible element of Fp[v] */

    /* same TODO */
    deg++;
    if (deg > 10000)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

    fq_nmod_mpoly_ctx_change_modulus(ectx, deg);

have_prime:

    /* make sure reduction does not kill both lc(A) and lc(B) */
    n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(gammaeval), gamma,
          evil_const_cast_nmod_poly_to_n_poly(ectx->fqctx->modulus), ctx->mod);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
    {
        goto choose_prime;
    }

    /* make sure reduction does not kill either A or B */
    nmod_mpolyn_interp_reduce_lg_mpolyn(Aeval, ectx, A, var, ctx);
    nmod_mpolyn_interp_reduce_lg_mpolyn(Beval, ectx, B, var, ctx);
    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(Aeval, ectx));
    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(Beval, ectx));
    if (Aeval->length == 0 || Beval->length == 0)
    {
        goto choose_prime;
    }

    success = fq_nmod_mpolyn_gcd_brown_smprime(Geval, Abareval, Bbareval,
                                                  Aeval, Beval, var - 1, ectx);
    if (success == 0)
    {
        goto choose_prime;
    }

    if (fq_nmod_mpolyn_is_nonzero_fq_nmod(Geval, ectx))
    {
        nmod_mpolyn_one(G, ctx);
        nmod_mpolyn_swap(Abar, A);
        nmod_mpolyn_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (n_poly_degree(modulus) > 0)
    {
        /* compare leading monomials of Geval and G */
        int cmp = 0;
        slong k;

        FLINT_ASSERT(G->length > 0);

        k = n_poly_degree(Geval->coeffs + 0);
        cmp = mpoly_monomial_cmp_nomask_extra(G->exps + N*0,
                                     Geval->exps + N*0, N, offset, k << shift);
        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            n_poly_one(modulus);
        }
    }

    n_fq_get_fq_nmod(temp, fq_nmod_mpolyn_leadcoeff(Geval, ectx), ectx->fqctx);
    fq_nmod_inv(temp, temp, ectx->fqctx);
    fq_nmod_mul(temp, temp, gammaeval, ectx->fqctx);
    fq_nmod_mpolyn_scalar_mul_fq_nmod(Geval, temp, ectx);

    if (n_poly_degree(modulus) > 0)
    {
        nmod_mpolyn_interp_crt_lg_mpolyn(&ldegG, G, T, modulus,
                                                        var, ctx, Geval, ectx);
        nmod_mpolyn_interp_crt_lg_mpolyn(&ldegAbar, Abar, T, modulus,
                                                     var, ctx, Abareval, ectx);
        nmod_mpolyn_interp_crt_lg_mpolyn(&ldegBbar, Bbar, T, modulus,
                                                     var, ctx, Bbareval, ectx);
    }
    else
    {
        nmod_mpolyn_interp_lift_lg_mpolyn(&ldegG, G, var, ctx, Geval, ectx);
        nmod_mpolyn_interp_lift_lg_mpolyn(&ldegAbar, Abar,
                                                     var, ctx, Abareval, ectx);
        nmod_mpolyn_interp_lift_lg_mpolyn(&ldegBbar, Bbar,
                                                     var, ctx, Bbareval, ectx);
    }
    n_poly_mod_mul(modulus, modulus,
          evil_const_cast_nmod_poly_to_n_poly(ectx->fqctx->modulus), ctx->mod);

    if (n_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    n_poly_one(modulus);
    goto choose_prime;

successful:

    nmod_mpolyn_content_last(modulus, G, ctx);
    nmod_mpolyn_divexact_last(G, modulus, ctx);
    nmod_mpolyn_divexact_last(Abar, nmod_mpolyn_leadcoeff_poly(G, ctx), ctx);
    nmod_mpolyn_divexact_last(Bbar, nmod_mpolyn_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    nmod_mpolyn_mul_last(G, cG, ctx);
    nmod_mpolyn_mul_last(Abar, cAbar, ctx);
    nmod_mpolyn_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == nmod_mpolyn_leadcoeff(G, ctx));
        n_poly_mod_mul(modulus, nmod_mpolyn_leadcoeff_poly(G, ctx),
                              nmod_mpolyn_leadcoeff_poly(Abar, ctx), ctx->mod);
        FLINT_ASSERT(n_poly_equal(modulus, leadA));
        n_poly_mod_mul(modulus, nmod_mpolyn_leadcoeff_poly(G, ctx),
                              nmod_mpolyn_leadcoeff_poly(Bbar, ctx), ctx->mod);
        FLINT_ASSERT(n_poly_equal(modulus, leadB));
    }
    n_poly_clear(leadA);
    n_poly_clear(leadB);
#endif

    n_poly_clear(cA);
    n_poly_clear(cB);
    n_poly_clear(cG);
    n_poly_clear(cAbar);
    n_poly_clear(cBbar);
    n_poly_clear(gamma);
    n_poly_clear(modulus);

    nmod_mpolyn_clear(T, ctx);

    fq_nmod_mpolyn_clear(Aeval, ectx);
    fq_nmod_mpolyn_clear(Beval, ectx);
    fq_nmod_mpolyn_clear(Geval, ectx);
    fq_nmod_mpolyn_clear(Abareval, ectx);
    fq_nmod_mpolyn_clear(Bbareval, ectx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(temp, ectx->fqctx);

    fq_nmod_mpoly_ctx_clear(ectx);

    return success;
}
