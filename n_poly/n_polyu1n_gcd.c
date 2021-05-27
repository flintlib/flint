/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "nmod_mpoly_factor.h"

static void n_polyu1n_mod_interp_reduce_2sm_poly(
    n_poly_t E,
    n_poly_t F,
    const n_polyun_t A,
    n_poly_t alphapow,
    nmod_t ctx)
{
    slong i;
    mp_limb_t u, v;

    n_poly_zero(E);
    n_poly_zero(F);
    for (i = 0; i < A->length; i++)
    {
        n_poly_mod_eval2_pow(&u, &v, A->coeffs + i, alphapow, ctx);
        n_poly_set_coeff(E, A->exps[i], u);
        n_poly_set_coeff(F, A->exps[i], v);
    }
}

static void n_polyu1n_mod_interp_lift_2sm_poly(
    slong * lastdeg,
    n_polyun_t F,
    const n_poly_t A,
    const n_poly_t B,
    mp_limb_t alpha,
    nmod_t ctx)
{
    slong lastlen = 0;
    slong Fi, Aexp, Bexp;
    const mp_limb_t * Acoeffs = A->coeffs;
    const mp_limb_t * Bcoeffs = B->coeffs;
    slong e;
    mp_limb_t d0 = (1 + ctx.n)/2;
    mp_limb_t d1 = nmod_inv(nmod_add(alpha, alpha, ctx), ctx);
    mp_limb_t Avalue, Bvalue, u, v;

    Aexp = n_poly_degree(A);
    Bexp = n_poly_degree(B);

    n_polyun_fit_length(F, FLINT_MAX(Aexp, Bexp) + 1);

    Fi = 0;
    while (Aexp >= 0 || Bexp >= 0)
    {
        e = Aexp;
        Avalue = 0;
        Bvalue = 0;
        if (Aexp == Bexp)
        {
            Avalue = Acoeffs[Aexp];
            Bvalue = Bcoeffs[Bexp];
        }
        else if (Aexp > Bexp)
        {
            Avalue = Acoeffs[Aexp];
        }
        else
        {
            FLINT_ASSERT(Bexp > Aexp);
            e = Bexp;
            Bvalue = Bcoeffs[Bexp];
        }
        FLINT_ASSERT(Avalue != 0 || Bvalue != 0);
        u = nmod_add(Avalue, Bvalue, ctx);
        v = nmod_sub(Avalue, Bvalue, ctx);
        u = nmod_mul(u, d0, ctx);
        v = nmod_mul(v, d1, ctx);

        FLINT_ASSERT(Fi < F->alloc);

        F->exps[Fi] = e;

        FLINT_ASSERT(u != 0 || v != 0);
        n_poly_fit_length(F->coeffs + Fi, 2);
        F->coeffs[Fi].coeffs[0] = u;
        F->coeffs[Fi].coeffs[1] = v;
        F->coeffs[Fi].length = 1 + (v != 0);
        lastlen = FLINT_MAX(lastlen, F->coeffs[Fi].length);
        Fi++;

        if (e == Aexp)
        {
            do {
                Aexp--;
            } while (Aexp >= 0 && Acoeffs[Aexp] == 0);
        }
        if (e == Bexp)
        {
            do {
                Bexp--;
            } while (Bexp >= 0 && Bcoeffs[Bexp] == 0);
        }
    }
    F->length = Fi;

    *lastdeg = lastlen - 1;
    return;
}

static int n_polyu1n_mod_interp_crt_2sm_poly(
    slong * lastdeg,
    n_polyun_t F,
    n_polyun_t T,
    const n_poly_t A,
    const n_poly_t B,
    const n_poly_t modulus,
    n_poly_t alphapow,
    nmod_t ctx)
{
    int changed = 0, Finc;
    slong lastlen = 0;
    n_poly_struct * Fvalue;
    mp_limb_t u, v, FvalueA, FvalueB;
    slong Fi, Ti, Aexp, Bexp, e, fexp;
    const mp_limb_t * Acoeff = A->coeffs;
    const mp_limb_t * Bcoeff = B->coeffs;
    slong Flen = F->length;
    n_poly_t zero;

    zero->alloc = 0;
    zero->length = 0;
    zero->coeffs = NULL;

    Fi = 0;
    Aexp = n_poly_degree(A);
    Bexp = n_poly_degree(B);

    n_polyun_fit_length(T, Flen + FLINT_MAX(Aexp, Bexp) + 1);
    Ti = 0;

#if FLINT_WANT_ASSERT
    u = n_poly_mod_evaluate_nmod(modulus, alphapow->coeffs[1], ctx);
    u = nmod_mul(u, alphapow->coeffs[1], ctx);
    u = nmod_add(u, u, ctx);
    FLINT_ASSERT(u == 1);

    v = nmod_neg(alphapow->coeffs[1], ctx);
    u = n_poly_mod_evaluate_nmod(modulus, v, ctx);
    u = nmod_mul(u, alphapow->coeffs[1], ctx);
    u = nmod_add(u, u, ctx);
    FLINT_ASSERT(u == 1);
#endif

    while (Fi < Flen || Aexp >= 0 || Bexp >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        fexp = e = -WORD(1);
        if (Fi < Flen)
        {
            fexp = e = F->exps[Fi];
            FLINT_ASSERT(!n_poly_is_zero(F->coeffs + Fi));
            FLINT_ASSERT(n_poly_degree(F->coeffs + Fi) < n_poly_degree(modulus));
        }
        if (Aexp >= 0)
        {
            e = FLINT_MAX(e, Aexp);
            FLINT_ASSERT(Acoeff[Aexp] != 0);
        }
        if (Bexp >= 0)
        {
            e = FLINT_MAX(e, Bexp);
            FLINT_ASSERT(Bcoeff[Bexp] != 0);
        }

        FLINT_ASSERT(e >= 0);

        T->exps[Ti] = e;

        Fvalue = zero;
        FvalueA = 0;
        FvalueB = 0;
        Finc = 0;
        if (Fi < Flen && e == fexp)
        {
            Finc = 1;
            Fvalue = F->coeffs + Fi;
            n_poly_mod_eval2_pow(&FvalueA, &FvalueB, Fvalue, alphapow, ctx);
        }

        if (e == Aexp)
        {
            FvalueA = nmod_sub(FvalueA, Acoeff[Aexp], ctx);
        }
        if (e == Bexp)
        {
            FvalueB = nmod_sub(FvalueB, Bcoeff[Bexp], ctx);
        }

        u = nmod_sub(FvalueB, FvalueA, ctx);
        v = nmod_add(FvalueB, FvalueA, ctx);
        v = nmod_mul(v, alphapow->coeffs[1], ctx);
        v = nmod_neg(v, ctx);
        changed |= u != 0 || v != 0;
        n_poly_mod_addmul_linear(T->coeffs + Ti, Fvalue, modulus, u, v, ctx);

        FLINT_ASSERT(T->coeffs[Ti].length > 0);
        lastlen = FLINT_MAX(lastlen, T->coeffs[Ti].length);
        Ti++;

        Fi += Finc;
        if (e == Aexp)
        {
            do {
                Aexp--;
            } while (Aexp >= 0 && Acoeff[Aexp] == 0);
        }
        if (e == Bexp)
        {
            do {
                Bexp--;
            } while (Bexp >= 0 && Bcoeff[Bexp] == 0);
        }
    }
    T->length = Ti;

    *lastdeg = lastlen - 1;

    if (changed)
        n_polyun_swap(T, F);

    return changed;
}


int n_polyu1n_mod_gcd_brown_smprime(
    n_polyun_t G,
    n_polyun_t Abar,
    n_polyun_t Bbar,
    n_polyun_t A,
    n_polyun_t B,
    nmod_t ctx,
    n_poly_polyun_stack_t St)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaevalp, gammaevalm;
    n_poly_struct * Aevalp, * Bevalp, * Gevalp, * Abarevalp, * Bbarevalp;
    n_poly_struct * Aevalm, * Bevalm, * Gevalm, * Abarevalm, * Bbarevalm;
    n_polyun_struct * T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    n_poly_struct * cA, * cB, * cG, * cAbar, * cBbar, * gamma;
    n_poly_struct * modulus, * alphapow, * r;
    int gstab, astab, bstab, use_stab;
#if FLINT_WANT_ASSERT
    n_poly_t leadA, leadB;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    n_poly_init(leadA);
    n_poly_init(leadB);
    n_poly_set(leadA, A->coeffs + 0);
    n_poly_set(leadB, B->coeffs + 0);
#endif

    n_poly_stack_fit_request(St->poly_stack, 19);
    cA          = n_poly_stack_take_top(St->poly_stack);
    cB          = n_poly_stack_take_top(St->poly_stack);
    cG          = n_poly_stack_take_top(St->poly_stack);
    cAbar       = n_poly_stack_take_top(St->poly_stack);
    cBbar       = n_poly_stack_take_top(St->poly_stack);
    gamma       = n_poly_stack_take_top(St->poly_stack);
    Aevalp      = n_poly_stack_take_top(St->poly_stack);
    Bevalp      = n_poly_stack_take_top(St->poly_stack);
    Gevalp      = n_poly_stack_take_top(St->poly_stack);
    Abarevalp   = n_poly_stack_take_top(St->poly_stack);
    Bbarevalp   = n_poly_stack_take_top(St->poly_stack);
    Aevalm      = n_poly_stack_take_top(St->poly_stack);
    Bevalm      = n_poly_stack_take_top(St->poly_stack);
    Gevalm      = n_poly_stack_take_top(St->poly_stack);
    Abarevalm   = n_poly_stack_take_top(St->poly_stack);
    Bbarevalm   = n_poly_stack_take_top(St->poly_stack);
    r           = n_poly_stack_take_top(St->poly_stack);
    alphapow    = n_poly_stack_take_top(St->poly_stack);
    modulus     = n_poly_stack_take_top(St->poly_stack);

    n_polyun_stack_fit_request(St->polyun_stack, 1);
    T           = n_polyun_stack_take_top(St->polyun_stack);

    _n_poly_vec_mod_remove_content(cA, A->coeffs, A->length, ctx);
    _n_poly_vec_mod_remove_content(cB, B->coeffs, B->length, ctx);

    n_poly_mod_gcd(cG, cA, cB, ctx);
    n_poly_mod_div(cAbar, cA, cG, ctx);
    n_poly_mod_div(cBbar, cB, cG, ctx);

    n_poly_mod_gcd(gamma, A->coeffs + 0, B->coeffs + 0, ctx);

    ldegA = _n_poly_vec_max_degree(A->coeffs, A->length);
    ldegB = _n_poly_vec_max_degree(B->coeffs, B->length);
    deggamma = n_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    n_poly_one(modulus);

    use_stab = 1;
    gstab = bstab = astab = 0;

    alpha = (ctx.n - 1)/2;

choose_prime: /* primes are v - alpha, v + alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha < ctx.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_poly_mod_eval2_pow(&gammaevalp, &gammaevalm, gamma, alphapow, ctx);
    if (gammaevalp == 0 || gammaevalm == 0)
        goto choose_prime;

    /* evaluation point should kill neither A nor B */
    n_polyu1n_mod_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, ctx);
    n_polyu1n_mod_interp_reduce_2sm_poly(Bevalp, Bevalm, B, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    if (use_stab && gstab)
    {
        slong Gdeg;
        n_polyu1n_mod_interp_reduce_2sm_poly(Gevalp, Gevalm, G, alphapow, ctx);
        Gdeg = G->exps[0];
        success = 1;
        success = success && n_poly_degree(Gevalp) == Gdeg;
        success = success && n_poly_degree(Gevalm) == Gdeg;
        success = success && Gevalp->coeffs[Gdeg] == gammaevalp;
        success = success && Gevalm->coeffs[Gdeg] == gammaevalm;
        success = success && (n_poly_mod_divrem(Abarevalp, r,
                                         Aevalp, Gevalp, ctx), r->length == 0);
        success = success && (n_poly_mod_divrem(Abarevalm, r,
                                         Aevalm, Gevalm, ctx), r->length == 0);
        success = success && (n_poly_mod_divrem(Bbarevalp, r,
                                         Bevalp, Gevalp, ctx), r->length == 0);
        success = success && (n_poly_mod_divrem(Bbarevalm, r,
                                         Bevalm, Gevalm, ctx), r->length == 0);
        if (!success)
        {
            use_stab = 0;
            n_poly_one(modulus);
            goto choose_prime;
        }

        _n_poly_mod_scalar_mul_nmod_inplace(Abarevalp, gammaevalp, ctx);
        _n_poly_mod_scalar_mul_nmod_inplace(Abarevalm, gammaevalm, ctx);
        _n_poly_mod_scalar_mul_nmod_inplace(Bbarevalp, gammaevalp, ctx);
        _n_poly_mod_scalar_mul_nmod_inplace(Bbarevalm, gammaevalm, ctx);
    }
    else
    {
        n_poly_mod_gcd(Gevalp, Aevalp, Bevalp, ctx);
        n_poly_mod_div(Abarevalp, Aevalp, Gevalp, ctx);
        n_poly_mod_div(Bbarevalp, Bevalp, Gevalp, ctx);
        n_poly_mod_gcd(Gevalm, Aevalm, Bevalm, ctx);
        n_poly_mod_div(Abarevalm, Aevalm, Gevalm, ctx);
        n_poly_mod_div(Bbarevalm, Bevalm, Gevalm, ctx);
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
        n_polyun_one(G);
        n_polyun_swap(Abar, A);
        n_polyun_swap(Bbar, B);
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
        if (n_poly_degree(Gevalp) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (n_poly_degree(Gevalp) < G->exps[0])
        {
            n_poly_one(modulus);
        }
    }

    /* update interpolants */
    _n_poly_mod_scalar_mul_nmod_inplace(Gevalp, gammaevalp, ctx);
    _n_poly_mod_scalar_mul_nmod_inplace(Gevalm, gammaevalm, ctx);
    if (n_poly_degree(modulus) > 0)
    {
        temp = n_poly_mod_evaluate_nmod(modulus, alpha, ctx);
        FLINT_ASSERT(temp == n_poly_mod_evaluate_nmod(modulus, ctx.n - alpha, ctx));
        temp = nmod_mul(temp, alpha, ctx);
        temp = nmod_add(temp, temp, ctx);
        temp = nmod_inv(temp, ctx);
        _n_poly_mod_scalar_mul_nmod_inplace(modulus, temp, ctx);
        gstab = gstab || !n_polyu1n_mod_interp_crt_2sm_poly(&ldegG, G, T, Gevalp, Gevalm, modulus, alphapow, ctx);
        n_polyu1n_mod_interp_crt_2sm_poly(&ldegAbar, Abar, T, Abarevalp, Abarevalm, modulus, alphapow, ctx);
        n_polyu1n_mod_interp_crt_2sm_poly(&ldegBbar, Bbar, T, Bbarevalp, Bbarevalm, modulus, alphapow, ctx);
    }
    else
    {
        n_polyu1n_mod_interp_lift_2sm_poly(&ldegG, G, Gevalp, Gevalm, alpha, ctx);
        n_polyu1n_mod_interp_lift_2sm_poly(&ldegAbar, Abar, Abarevalp, Abarevalm, alpha, ctx);
        n_polyu1n_mod_interp_lift_2sm_poly(&ldegBbar, Bbar, Bbarevalp, Bbarevalm, alpha, ctx);
        gstab = astab = bstab = 0;
    }

    temp = ctx.n - nmod_mul(alpha, alpha, ctx);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, temp, ctx);

    if (n_poly_degree(modulus) < bound)
        goto choose_prime;

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (deggamma + ldegA == ldegG + ldegAbar &&
        deggamma + ldegB == ldegG + ldegBbar)
    {
        goto successful;
    }

    n_poly_one(modulus);
    goto choose_prime;

successful:

    _n_poly_vec_mod_remove_content(modulus, G->coeffs, G->length, ctx);
    _n_poly_vec_mod_divexact_poly(Abar->coeffs, Abar->length, G->coeffs + 0, ctx);
    _n_poly_vec_mod_divexact_poly(Bbar->coeffs, Bbar->length, G->coeffs + 0, ctx);

successful_put_content:

    _n_poly_vec_mod_mul_poly(G->coeffs, G->length, cG, ctx);
    _n_poly_vec_mod_mul_poly(Abar->coeffs, Abar->length, cAbar, ctx);
    _n_poly_vec_mod_mul_poly(Bbar->coeffs, Bbar->length, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(1 == n_poly_lead(G->coeffs + 0));
        n_poly_mod_mul(modulus, G->coeffs + 0, Abar->coeffs + 0, ctx);
        FLINT_ASSERT(n_poly_equal(modulus, leadA));
        n_poly_mod_mul(modulus, G->coeffs + 0, Bbar->coeffs + 0, ctx);
        FLINT_ASSERT(n_poly_equal(modulus, leadB));
    }
    n_poly_clear(leadA);
    n_poly_clear(leadB);
#endif

    n_poly_stack_give_back(St->poly_stack, 19);
    n_polyun_stack_give_back(St->polyun_stack, 1);

    return success;
}

