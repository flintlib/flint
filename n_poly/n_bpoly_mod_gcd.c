/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


void n_bpoly_mod_interp_reduce_2sm_poly(
    n_poly_t Ap,
    n_poly_t Am,
    const n_bpoly_t A,
    n_poly_t alphapow,
    nmod_t mod)
{
    slong i, Alen = A->length;
    const n_poly_struct * Ac = A->coeffs;
    mp_limb_t * Apc, * Amc;

    n_poly_fit_length(Ap, Alen);
    n_poly_fit_length(Am, Alen);

    Apc = Ap->coeffs;
    Amc = Am->coeffs;

    for (i = 0; i < Alen; i++)
        n_poly_mod_eval2_pow(Apc + i, Amc + i, Ac + i, alphapow, mod);

    Ap->length = Alen;
    _n_poly_normalise(Ap);
    Am->length = Alen;
    _n_poly_normalise(Am);
}

void n_bpoly_mod_interp_lift_2sm_poly(
    slong * deg1,
    n_bpoly_t T,
    const n_poly_t A,
    const n_poly_t B,    
    mp_limb_t alpha,
    nmod_t mod)
{
    slong i;
    slong lastlength = 0;
    const mp_limb_t * Acoeffs = A->coeffs;
    const mp_limb_t * Bcoeffs = B->coeffs;
    n_poly_struct * Tcoeffs;
    slong Alen = A->length;
    slong Blen = B->length;
    slong Tlen = FLINT_MAX(Alen, Blen);
    mp_limb_t d0 = (1 + mod.n)/2;
    mp_limb_t d1 = nmod_inv(nmod_add(alpha, alpha, mod), mod);
    mp_limb_t Avalue, Bvalue, u, v;

    n_bpoly_fit_length(T, Tlen);

    Tcoeffs = T->coeffs;

    for (i = 0; i < Tlen; i++)
    {
        Avalue = (i < Alen) ? Acoeffs[i] : 0;
        Bvalue = (i < Blen) ? Bcoeffs[i] : 0;
        u = nmod_sub(Avalue, Bvalue, mod);
        v = nmod_add(Avalue, Bvalue, mod);
        u = nmod_mul(u, d1, mod);
        v = nmod_mul(v, d0, mod);
        if ((u | v) == 0)
        {
            n_poly_zero(Tcoeffs + i);
        }
        else
        {
            n_poly_fit_length(Tcoeffs + i, 2);
            Tcoeffs[i].coeffs[0] = v;
            Tcoeffs[i].coeffs[1] = u;
            Tcoeffs[i].length = 1 + (u != 0);
            lastlength = FLINT_MAX(lastlength, Tcoeffs[i].length);
        }
    }

    *deg1 = lastlength - 1;

    FLINT_ASSERT(Tlen <= 0 || !n_poly_is_zero(Tcoeffs + Tlen - 1));
    T->length = Tlen;
}

int n_bpoly_mod_interp_crt_2sm_poly(
    slong * deg1,
    n_bpoly_t F,
    n_bpoly_t T,
    n_poly_t A,
    n_poly_t B,
    const n_poly_t modulus,
    n_poly_t alphapow,
    nmod_t mod)
{
    int changed = 0;
    slong i, lastlength = 0;
    slong Alen = A->length;
    slong Blen = B->length;
    slong Flen = F->length;
    slong Tlen = FLINT_MAX(FLINT_MAX(Alen, Blen), Flen);
    n_poly_struct * Tcoeffs, * Fcoeffs;
    mp_limb_t * Acoeffs, * Bcoeffs;
    n_poly_t zero;
    mp_limb_t Avalue, Bvalue, FvalueA, FvalueB, u, v;
    n_poly_struct * Fvalue;
    mp_limb_t alpha = alphapow->coeffs[1];

    zero->alloc = 0;
    zero->length = 0;
    zero->coeffs = NULL;

    n_bpoly_fit_length(T, Tlen);
    Tcoeffs = T->coeffs;
    Acoeffs = A->coeffs;
    Bcoeffs = B->coeffs;
    Fcoeffs = F->coeffs;

    for (i = 0; i < Tlen; i++)
    {
        Fvalue = (i < Flen) ? Fcoeffs + i : zero;
        n_poly_mod_eval2_pow(&FvalueA, &FvalueB, Fvalue, alphapow, mod);
        Avalue = (i < Alen) ? Acoeffs[i] : 0;
        Bvalue = (i < Blen) ? Bcoeffs[i] : 0;
        FvalueA = nmod_sub(FvalueA, Avalue, mod);
        FvalueB = nmod_sub(FvalueB, Bvalue, mod);
        u = nmod_sub(FvalueB, FvalueA, mod);
        v = nmod_mul(mod.n - alpha, nmod_add(FvalueB, FvalueA, mod), mod);
        if ((u | v) != 0)
        {
            changed = 1;
            n_poly_mod_addmul_linear(Tcoeffs + i, Fvalue, modulus, u, v, mod);
        }
        else
        {
            n_poly_set(Tcoeffs + i, Fvalue);
        }

        lastlength = FLINT_MAX(lastlength, Tcoeffs[i].length);
    }

    FLINT_ASSERT(Tlen <= 0 || !n_poly_is_zero(Tcoeffs + Tlen - 1));
    T->length = Tlen;

    if (changed)
        n_bpoly_swap(T, F);

    FLINT_ASSERT(n_bpoly_mod_is_canonical(F, mod));

    *deg1 = lastlength - 1;
    return changed;
}


int n_bpoly_mod_gcd_brown_smprime(
    n_bpoly_t G,
    n_bpoly_t Abar,
    n_bpoly_t Bbar,
    n_bpoly_t A,
    n_bpoly_t B,
    nmod_t ctx,
    n_poly_bpoly_stack_t Sp)
{
    int success;
    slong bound;
    mp_limb_t alpha, temp, gammaevalp, gammaevalm;
    n_poly_struct * Aevalp, * Bevalp, * Gevalp, * Abarevalp, * Bbarevalp;
    n_poly_struct * Aevalm, * Bevalm, * Gevalm, * Abarevalm, * Bbarevalm;
    n_bpoly_struct * T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    n_poly_struct * cA, * cB, * cG, * cAbar, * cBbar, * gamma;
    n_poly_struct * modulus, * alphapow, * r;
    int gstab, astab, bstab, use_stab;
#if FLINT_WANT_ASSERT
    n_poly_t leadA, leadB;
    const slong Sp_size_poly = n_poly_stack_size(Sp->poly_stack);
    const slong Sp_size_bpoly = n_bpoly_stack_size(Sp->bpoly_stack);
#endif

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

#if FLINT_WANT_ASSERT
    n_poly_init(leadA);
    n_poly_init(leadB);
    n_poly_set(leadA, A->coeffs + A->length - 1);
    n_poly_set(leadB, B->coeffs + B->length - 1);
#endif

    n_poly_stack_fit_request(Sp->poly_stack, 19);
    cA          = n_poly_stack_take_top(Sp->poly_stack);
    cB          = n_poly_stack_take_top(Sp->poly_stack);
    cG          = n_poly_stack_take_top(Sp->poly_stack);
    cAbar       = n_poly_stack_take_top(Sp->poly_stack);
    cBbar       = n_poly_stack_take_top(Sp->poly_stack);
    gamma       = n_poly_stack_take_top(Sp->poly_stack);
    Aevalp      = n_poly_stack_take_top(Sp->poly_stack);
    Bevalp      = n_poly_stack_take_top(Sp->poly_stack);
    Gevalp      = n_poly_stack_take_top(Sp->poly_stack);
    Abarevalp   = n_poly_stack_take_top(Sp->poly_stack);
    Bbarevalp   = n_poly_stack_take_top(Sp->poly_stack);
    Aevalm      = n_poly_stack_take_top(Sp->poly_stack);
    Bevalm      = n_poly_stack_take_top(Sp->poly_stack);
    Gevalm      = n_poly_stack_take_top(Sp->poly_stack);
    Abarevalm   = n_poly_stack_take_top(Sp->poly_stack);
    Bbarevalm   = n_poly_stack_take_top(Sp->poly_stack);
    r           = n_poly_stack_take_top(Sp->poly_stack);
    alphapow    = n_poly_stack_take_top(Sp->poly_stack);
    modulus     = n_poly_stack_take_top(Sp->poly_stack);

    n_bpoly_stack_fit_request(Sp->bpoly_stack, 1);
    T           = n_bpoly_stack_take_top(Sp->bpoly_stack);

    n_bpoly_mod_content_last(cA, A, ctx);
    n_bpoly_mod_content_last(cB, B, ctx);
    n_bpoly_mod_divexact_last(A, cA, ctx);
    n_bpoly_mod_divexact_last(B, cB, ctx);

    n_poly_mod_gcd(cG, cA, cB, ctx);

    n_poly_mod_div(cAbar, cA, cG, ctx);
    n_poly_mod_div(cBbar, cB, cG, ctx);

    n_poly_mod_gcd(gamma, A->coeffs + A->length - 1, B->coeffs + B->length - 1, ctx);

    ldegA = n_bpoly_degree1(A);
    ldegB = n_bpoly_degree1(B);
    deggamma = n_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    n_poly_one(modulus);

    use_stab = 1;
    gstab = bstab = astab = 0;

    if ((ctx.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (ctx.n - UWORD(1))/UWORD(2);

choose_prime: /* primes are v - alpha, v + alpha */

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_poly_mod_eval2_pow(&gammaevalp, &gammaevalm, gamma, alphapow, ctx);
    if (gammaevalp == 0 || gammaevalm == 0)
        goto choose_prime;

    n_bpoly_mod_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, ctx);
    n_bpoly_mod_interp_reduce_2sm_poly(Bevalp, Bevalm, B, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    if (use_stab && gstab)
    {
        slong Gdeg;
        n_bpoly_mod_interp_reduce_2sm_poly(Gevalp, Gevalm, G, alphapow, ctx);
        Gdeg = n_bpoly_degree0(G);
        success = 1;
        success = success && n_poly_degree(Gevalp) == Gdeg;
        success = success && n_poly_degree(Gevalm) == Gdeg;
        success = success && Gevalp->coeffs[Gdeg] == gammaevalp;
        success = success && Gevalm->coeffs[Gdeg] == gammaevalm;
        n_poly_mod_divrem(Abarevalp, r, Aevalp, Gevalp, ctx);
        success = success && (r->length == 0);
        n_poly_mod_divrem(Abarevalm, r, Aevalm, Gevalm, ctx);
        success = success && (r->length == 0);
        n_poly_mod_divrem(Bbarevalp, r, Bevalp, Gevalp, ctx);
        success = success && (r->length == 0);
        n_poly_mod_divrem(Bbarevalm, r, Bevalm, Gevalm, ctx);
        success = success && (r->length == 0);

        if (!success)
        {
            use_stab = 0;
            n_poly_one(modulus);
            alpha = (ctx.n - UWORD(1))/UWORD(2);
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
        n_bpoly_one(G);
        n_bpoly_swap(Abar, A);
        n_bpoly_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (n_poly_degree(Gevalp) != n_poly_degree(Gevalm))
    {
        goto choose_prime;
    }

    if (n_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (n_poly_degree(Gevalp) > n_bpoly_degree0(G))
        {
            goto choose_prime;
        }
        else if (n_poly_degree(Gevalp) < n_bpoly_degree0(G))
        {
            n_poly_one(modulus);
        }
    }

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
        if (!gstab)
        {
            gstab = !n_bpoly_mod_interp_crt_2sm_poly(&ldegG, G, T,
                                       Gevalp, Gevalm, modulus, alphapow, ctx);
        }
        n_bpoly_mod_interp_crt_2sm_poly(&ldegAbar, Abar, T,
                                 Abarevalp, Abarevalm, modulus, alphapow, ctx);
        n_bpoly_mod_interp_crt_2sm_poly(&ldegBbar, Bbar, T,
                                 Bbarevalp, Bbarevalm, modulus, alphapow, ctx);
    }
    else
    {
        n_bpoly_mod_interp_lift_2sm_poly(&ldegG, G, Gevalp, Gevalm, alpha, ctx);
        n_bpoly_mod_interp_lift_2sm_poly(&ldegAbar, Abar,
                                             Abarevalp, Abarevalm, alpha, ctx);
        n_bpoly_mod_interp_lift_2sm_poly(&ldegBbar, Bbar,
                                             Bbarevalp, Bbarevalm, alpha, ctx);
        gstab = astab = bstab = 0;
    }

    temp = ctx.n - nmod_mul(alpha, alpha, ctx);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, temp, ctx);

    if (n_poly_degree(modulus) < bound)
        goto choose_prime;

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

    n_bpoly_mod_content_last(modulus, G, ctx);
    n_bpoly_mod_divexact_last(G, modulus, ctx);
    n_bpoly_mod_divexact_last(Abar, G->coeffs + G->length - 1, ctx);
    n_bpoly_mod_divexact_last(Bbar, G->coeffs + G->length - 1, ctx);

successful_put_content:

    n_bpoly_mod_mul_last(G, cG, ctx);
    n_bpoly_mod_mul_last(Abar, cAbar, ctx);
    n_bpoly_mod_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        n_poly_struct * Glead = G->coeffs + G->length - 1;
        FLINT_ASSERT(1 == Glead->coeffs[Glead->length - 1]);
        n_poly_mod_mul(modulus, G->coeffs + G->length - 1, Abar->coeffs + Abar->length - 1, ctx);
        FLINT_ASSERT(n_poly_equal(modulus, leadA));
        n_poly_mod_mul(modulus, G->coeffs + G->length - 1, Bbar->coeffs + Bbar->length - 1, ctx);
        FLINT_ASSERT(n_poly_equal(modulus, leadB));
    }
    n_poly_clear(leadA);
    n_poly_clear(leadB);
#endif

    n_poly_stack_give_back(Sp->poly_stack, 19);
    n_bpoly_stack_give_back(Sp->bpoly_stack, 1);
    FLINT_ASSERT(Sp_size_poly == n_poly_stack_size(Sp->poly_stack));
    FLINT_ASSERT(Sp_size_bpoly == n_bpoly_stack_size(Sp->bpoly_stack));

    return success;
}

