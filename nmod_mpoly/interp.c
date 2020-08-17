/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/*
    interp_reduce: map from Fp[x] to Fp[x]/poly(x)
    interp_lift:   map from Fp[x]/poly(x) to Fp[x]
    interp_crt:    update element of Fp[x] with a new image in Fp[x]/poly(x)
    interp_mcrt:   same as interp_crt, but monomial match, thus easier
*/


void _nmod_poly_eval2_pow(
    mp_limb_t * vp,
    mp_limb_t * vm,
    nmod_poly_t P, 
    nmod_poly_t alphapow,
    const nmodf_ctx_t fctx)
{
    mp_limb_t * Pcoeffs = P->coeffs;
    slong Plen = P->length;
    mp_limb_t * alpha_powers = alphapow->coeffs;
    mp_limb_t p1, p0, a0, a1, a2, q1, q0, b0, b1, b2;
    slong k;

    a0 = a1 = a2 = UWORD(0);
    b0 = b1 = b2 = UWORD(0);

    if (Plen > alphapow->length)
    {
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        nmod_poly_fit_length(alphapow, Plen);
        for (k = oldlength; k < Plen; k++)
        {
            alphapow->coeffs[k] = nmod_mul(alphapow->coeffs[k - 1],
                                               alphapow->coeffs[1], fctx->mod);
        }
        alphapow->length = Plen;
        alpha_powers = alphapow->coeffs;
    }

    for (k = 0; k + 2 <= Plen; k += 2)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        umul_ppmm(q1, q0, Pcoeffs[k + 1], alpha_powers[k + 1]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, WORD(0), p1, p0);
        add_sssaaaaaa(b2, b1, b0, b2, b1, b0, WORD(0), q1, q0);
    }

    if (k < Plen)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, WORD(0), p1, p0);
        k++;
    }

    FLINT_ASSERT(k == Plen);

    NMOD_RED3(p0, a2, a1, a0, fctx->mod);
    NMOD_RED3(q0, b2, b1, b0, fctx->mod);

    *vp = nmod_add(p0, q0, fctx->mod);
    *vm = nmod_sub(p0, q0, fctx->mod);
}

/*****************************************************************************/

void nmod_mpolyn_interp_reduce_2sm_poly(
    nmod_poly_t E,
    nmod_poly_t F,
    const nmod_mpolyn_t A,
    nmod_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    mp_limb_t u, v;
    slong Ai, Alen, k;
    nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    nmod_poly_zero(E);
    nmod_poly_zero(F);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        _nmod_poly_eval2_pow(&u, &v, Acoeff + Ai, alphapow, ctx->ffinfo);
        k = (Aexp + N*Ai)[off] >> shift;
        nmod_poly_set_coeff_ui(E, k, u);
        nmod_poly_set_coeff_ui(F, k, v);
    }
}

void nmod_mpolyn_interp_lift_2sm_poly(
    slong * lastdeg_,
    nmod_mpolyn_t F,
    const nmod_poly_t A,
    const nmod_poly_t B,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong lastdeg = -WORD(1);
    mp_limb_t u, v, d0, d1, Avalue, Bvalue;
    slong Fi, Aexp, Bexp;
    mp_limb_t * Acoeff = A->coeffs;
    mp_limb_t * Bcoeff = B->coeffs;
    nmod_poly_struct * Fcoeff;
    ulong * Fexp;
    slong e;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    Aexp = nmod_poly_degree(A);
    Bexp = nmod_poly_degree(B);

    nmod_mpolyn_fit_length(F, FLINT_MAX(Aexp, Bexp) + 1, ctx);
    Fcoeff = F->coeffs;
    Fexp = F->exps;

    d0 = n_invmod(UWORD(2), ctx->ffinfo->mod.n);
    d1 = n_invmod(nmod_add(alpha, alpha, ctx->ffinfo->mod), ctx->ffinfo->mod.n);

    Fi = 0;
    while (Aexp >= 0 || Bexp >= 0)
    {
        e = Aexp;
        Avalue = 0;
        Bvalue = 0;
        if (Aexp == Bexp)
        {
            Avalue = Acoeff[Aexp];
            Bvalue = Bcoeff[Bexp];
        }
        else if (Aexp > Bexp)
        {
            Avalue = Acoeff[Aexp];
        }
        else
        {
            FLINT_ASSERT(Bexp > Aexp);
            e = Bexp;
            Bvalue = Bcoeff[Bexp];
        }
        FLINT_ASSERT(Avalue != 0 || Bvalue != 0);
        u = nmod_add(Avalue, Bvalue, ctx->ffinfo->mod);
        v = nmod_sub(Avalue, Bvalue, ctx->ffinfo->mod);
        u = nmod_mul(u, d0, ctx->ffinfo->mod);
        v = nmod_mul(v, d1, ctx->ffinfo->mod);

        FLINT_ASSERT(Fi < F->alloc);
        mpoly_monomial_zero(Fexp + N*Fi, N);
        (Fexp + N*Fi)[off] = e << shift;

        FLINT_ASSERT(u != 0 || v != 0);
        nmod_poly_fit_length(Fcoeff + Fi, 2);
        (Fcoeff + Fi)->coeffs[0] = u;
        (Fcoeff + Fi)->coeffs[1] = v;
        (Fcoeff + Fi)->length = 1 + (v != 0);
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Fcoeff + Fi));
        Fi++;

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
    F->length = Fi;

    *lastdeg_ = lastdeg;
    return;
}

int nmod_mpolyn_interp_crt_2sm_poly(
    slong * lastdeg_,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    const nmod_poly_t A,
    const nmod_poly_t B,
    const nmod_poly_t modulus,
    nmod_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0, Finc;
    mp_limb_t alpha = nmod_poly_get_coeff_ui(alphapow, 1);
    slong lastdeg = -WORD(1);
    mp_limb_t u, v, FvalueA, FvalueB;
    slong Fi, Toff, Aexp, Bexp, e, fexp;
    mp_limb_t * Acoeff = A->coeffs;
    mp_limb_t * Bcoeff = B->coeffs;
    slong Flen = F->length;
    nmod_poly_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong N, off, shift;

    FLINT_ASSERT(T->bits == F->bits);

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    Fi = 0;
    Aexp = nmod_poly_degree(A);
    Bexp = nmod_poly_degree(B);

    nmod_mpolyn_fit_length(T, Flen + FLINT_MAX(Aexp, Bexp) + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Toff = 0;

#if WANT_ASSERT
    u = nmod_poly_evaluate_nmod(modulus, alpha);
    u = nmod_mul(u, alpha, ctx->ffinfo->mod);
    u = nmod_mul(u, 2, ctx->ffinfo->mod);
    FLINT_ASSERT(u == 1);
    u = nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha);
    u = nmod_mul(u, alpha, ctx->ffinfo->mod);
    u = nmod_mul(u, 2, ctx->ffinfo->mod);
    FLINT_ASSERT(u == 1);
#endif

    while (Fi < Flen || Aexp >= 0 || Bexp >= 0)
    {
        FLINT_ASSERT(Toff < T->alloc);

        fexp = e = -WORD(1);
        if (Fi < Flen)
        {
            fexp = e = (Fexp + N*Fi)[off]>>shift;
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeff + Fi));
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + Fi) < nmod_poly_degree(modulus));
        }
        if (Aexp >= 0)
        {
            e = FLINT_MAX(e, Aexp);
            FLINT_ASSERT(Acoeff[Aexp] != UWORD(0));
        }
        if (Bexp >= 0)
        {
            e = FLINT_MAX(e, Bexp);
            FLINT_ASSERT(Bcoeff[Bexp] != UWORD(0));
        }

        FLINT_ASSERT(e >= 0);

        mpoly_monomial_zero(Texp + N*Toff, N);
        (Texp + N*Toff)[off] = e << shift;

        FvalueA = FvalueB = 0;
        Finc = 0;
        if (Fi < Flen && e == fexp)
        {
            Finc = 1;
            _nmod_poly_eval2_pow(&FvalueA, &FvalueB, Fcoeff + Fi, alphapow, ctx->ffinfo);
        }

        if (e == Aexp)
        {
            FvalueA = nmod_sub(FvalueA, Acoeff[Aexp], ctx->ffinfo->mod);
        }
        if (e == Bexp)
        {
            FvalueB = nmod_sub(FvalueB, Bcoeff[Bexp], ctx->ffinfo->mod);
        }

        u = nmod_sub(FvalueB, FvalueA, ctx->ffinfo->mod);
        v = nmod_mul(ctx->ffinfo->mod.n - alpha, nmod_add(FvalueB, FvalueA, ctx->ffinfo->mod), ctx->ffinfo->mod);

        if (u != 0 || v != 0)
        {
            changed = 1;
            if (u != 0)
            {
                nmod_poly_scalar_mul_nmod(Tcoeff + Toff, modulus, u);
                nmod_poly_shift_left(Tcoeff + Toff, Tcoeff + Toff, 1);
                _nmod_vec_scalar_addmul_nmod((Tcoeff + Toff)->coeffs,
                            modulus->coeffs, modulus->length, v, ctx->ffinfo->mod);
            }
            else
            {
                nmod_poly_scalar_mul_nmod(Tcoeff + Toff, modulus, v);
            }

            if (Finc)
            {
                nmod_poly_add(Tcoeff + Toff, Tcoeff + Toff, Fcoeff + Fi);
            }
        }
        else
        {
            FLINT_ASSERT(Finc == 1);
            nmod_poly_set(Tcoeff + Toff, Fcoeff + Fi);
        }

        FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + Toff));
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Tcoeff + Toff));
        Toff++;

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
    T->length = Toff;

    if (changed)
    {
        nmod_mpolyn_swap(T, F);
    }

    *lastdeg_ = lastdeg;
    return changed;
}


/*****************************************************************************/

/*
    E = A(x_var = alpha)
    A is in Fp[x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyn_interp_reduce_sm_mpolyn(
    nmod_mpolyn_t E,
    nmod_mpolyn_t A,
    slong var,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    mp_limb_t v;
    nmod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        v = nmod_poly_evaluate_nmod(Acoeff + Ai, alpha);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (v == 0)
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            nmod_poly_set_coeff_ui(Ecoeff + Ei - 1, k, v);
        }
        else
        {
            FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

            /* create new */
            if (Ei >= E->alloc)
            {
                nmod_mpolyn_fit_length(E, Ei + 1, ctx);
                Ecoeff = E->coeffs;
                Eexp = E->exps;
            }
            mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
            nmod_poly_zero(Ecoeff + Ei);
            nmod_poly_set_coeff_ui(Ecoeff + Ei, k, v);
            Ei++;
        }
    }
    E->length = Ei;
}


/*
    A = B
    A is in Fp[x_0, ..., x_(var-1), x_(var-1)][x_var]
    B is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyn_interp_lift_sm_mpolyn(
    nmod_mpolyn_t A,
    nmod_mpolyn_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;

    nmod_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;

    nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;

    nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Ai = 0;
    for (Bi = 0; Bi < Blen; Bi++)
    {
        if (Ai + (Bcoeff + Bi)->length >= A->alloc)
        {
            nmod_mpolyn_fit_length(A, Ai + (Bcoeff + Bi)->length, ctx);
            Acoeff = A->coeffs;
            Aexp = A->exps;
        }
        for (vi = (Bcoeff + Bi)->length - 1; vi >= 0; vi--)
        {
            if ((Bcoeff + Bi)->coeffs[vi] != 0)
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                nmod_poly_zero(Acoeff + Ai);
                nmod_poly_set_coeff_ui(Acoeff + Ai, 0, (Bcoeff + Bi)->coeffs[vi]);
                Ai++;
            }
        }
    }
    A->length = Ai;
}


/*
    T = F + modulus*(A - F(x_var = alpha))
    no assumptions about matching monomials
    F is in Fp[x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
    in order to fxn correctly, modulus(alpha) should be 1
*/
int nmod_mpolyn_interp_crt_sm_mpolyn(
    slong * lastdeg_,
    nmod_mpolyn_t T,
    nmod_mpolyn_t F,
    nmod_mpolyn_t A,
    slong var,
    nmod_poly_t modulus,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    slong vi;
    mp_limb_t v;
    nmod_poly_t tp;

    nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    nmod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    nmod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    Fi = Ai = vi = 0;
    if (Ai < Alen)
    {
        vi = nmod_poly_degree(A->coeffs + Ai);
    }
    while (Fi < Flen || Ai < Alen)
    {
        if (Ti >= T->alloc)
        {
            nmod_mpolyn_fit_length(T, Ti + FLINT_MAX(Flen - Fi, Alen - Ai), ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        if (Fi < Flen && Ai < Alen && mpoly_monomial_equal_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift))
        {
            /* F term ok, A term ok */
            v = nmod_poly_evaluate_nmod(Fcoeff + Fi, alpha);
            v = nmod_sub((Acoeff + Ai)->coeffs[vi], v, ctx->ffinfo->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_add(Tcoeff + Ti, Fcoeff + Fi, tp);
            }
            else
            {
                nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && (Acoeff + Ai)->coeffs[vi] == 0);
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = nmod_poly_degree(A->coeffs + Ai);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            v = nmod_poly_evaluate_nmod(Fcoeff + Fi, alpha);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, tp);
            }
            else
            {
                nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            changed = 1;
            nmod_poly_scalar_mul_nmod(Tcoeff + Ti, modulus, (Acoeff + Ai)->coeffs[vi]);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && (Acoeff + Ai)->coeffs[vi] == 0);
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = nmod_poly_degree(A->coeffs + Ai);
                }
            }
        }

        FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    nmod_poly_clear(tp);

    *lastdeg_ = FLINT_MAX(*lastdeg_, lastdeg);
    return changed;
}


/****************************************************************************/

void nmod_mpolyn_interp_reduce_2sm_mpolyn(
    nmod_mpolyn_t E,
    nmod_mpolyn_t F,
    nmod_mpolyn_t A,
    slong var,
    nmod_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    mp_limb_t e, f;
    nmod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;
    nmod_poly_struct * Fcoeff;
    ulong * Fexp;
    slong Fi;

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    Fcoeff = F->coeffs;
    Fexp = F->exps;
    Fi = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        _nmod_poly_eval2_pow(&e, &f, Acoeff + Ai, alphapow, ctx->ffinfo);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;

        if (e != 0)
        {
            if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
            {
                /* append to previous */
                nmod_poly_set_coeff_ui(Ecoeff + Ei - 1, k, e);
            }
            else
            {
                FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

                /* create new */
                if (Ei >= E->alloc)
                {
                    nmod_mpolyn_fit_length(E, Ei + 1, ctx);
                    Ecoeff = E->coeffs;
                    Eexp = E->exps;
                }
                mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
                nmod_poly_zero(Ecoeff + Ei);
                nmod_poly_set_coeff_ui(Ecoeff + Ei, k, e);
                Ei++;
            }
        }

        if (f != 0)
        {
            if (Fi > 0 && mpoly_monomial_equal_extra(Fexp + N*(Fi - 1), Aexp + N*Ai, N, offset, -(k << shift)))
            {
                /* append to previous */
                nmod_poly_set_coeff_ui(Fcoeff + Fi - 1, k, f);
            }
            else
            {
                FLINT_ASSERT(Fi == 0 || mpoly_monomial_gt_nomask_extra(Fexp + N*(Fi - 1), Aexp + N*Ai, N, offset, -(k << shift)));

                /* create new */
                if (Fi >= F->alloc)
                {
                    nmod_mpolyn_fit_length(F, Fi + 1, ctx);
                    Fcoeff = F->coeffs;
                    Fexp = F->exps;
                }
                mpoly_monomial_set_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, -(k << shift));
                nmod_poly_zero(Fcoeff + Fi);
                nmod_poly_set_coeff_ui(Fcoeff + Fi, k, f);
                Fi++;
            }
        }
    }
    E->length = Ei;
    F->length = Fi;
}


void nmod_mpolyn_interp_lift_2sm_mpolyn(
    slong * lastdeg_,
    nmod_mpolyn_t T,
    nmod_mpolyn_t A,
    nmod_mpolyn_t B,
    slong var,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    nmod_poly_t tp;
    nmod_poly_t zero;

    nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai, ai;

    nmod_poly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    slong Bi, bi;

    mp_limb_t u, v, Avalue, Bvalue, FvalueA, FvalueB;
    int cmp;
    mp_limb_t d0 = n_invmod(alpha + alpha, ctx->ffinfo->mod.n);

    nmod_poly_init(tp, ctx->ffinfo->mod.n);
    nmod_poly_init(zero, ctx->ffinfo->mod.n);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    nmod_mpolyn_fit_length(T, FLINT_MAX(Alen, Blen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;

    Ti = Ai = Bi = 0;
    ai = (Ai >= Alen) ? 0 : nmod_poly_degree(A->coeffs + Ai);
    bi = (Bi >= Blen) ? 0 : nmod_poly_degree(B->coeffs + Bi);

    while (Ai < Alen || Bi < Blen)
    {
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Alen - Ai, Blen - Bi);
            nmod_mpolyn_fit_length(T, Ti + extra, ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        FLINT_ASSERT(Ai >= Alen || (Acoeff + Ai)->coeffs[ai] != 0);
        FLINT_ASSERT(Bi >= Blen || (Bcoeff + Bi)->coeffs[bi] != 0);

        Avalue = 0;
        if (Ai < Alen)
        {
            Avalue = (Acoeff + Ai)->coeffs[ai];
            mpoly_monomial_set_extra(Texp + N*Ti,
                                     Aexp + N*Ai, N, offset, ai << shift);
        }

        Bvalue = 0;
        if (Bi < Blen)
        {
            cmp = (Avalue == 0) ? -1 : mpoly_monomial_cmp_nomask_extra(
                             Texp + N*Ti, Bexp + N*Bi, N, offset, bi << shift);
            if (cmp <= 0)
            {
                Bvalue = (Bcoeff + Bi)->coeffs[bi];
            }
            if (cmp < 0)
            {
                Avalue = 0;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Bexp + N*Bi, N, offset, bi << shift);
            }
        }

        FvalueA = nmod_neg(Avalue, ctx->ffinfo->mod);
        FvalueB = nmod_neg(Bvalue, ctx->ffinfo->mod);
        u = nmod_sub(FvalueB, FvalueA, ctx->ffinfo->mod);
        v = nmod_mul(ctx->ffinfo->mod.n - alpha, nmod_add(FvalueB, FvalueA, ctx->ffinfo->mod), ctx->ffinfo->mod);

        FLINT_ASSERT(u != 0 || v != 0);
        nmod_poly_zero(Tcoeff + Ti);
        u = nmod_mul(u, d0, ctx->ffinfo->mod);
        v = nmod_mul(v, d0, ctx->ffinfo->mod);
        nmod_poly_set_coeff_ui(Tcoeff + Ti, 0, v);
        nmod_poly_set_coeff_ui(Tcoeff + Ti, 1, u);

        if (Avalue != 0)
        {
            do {
                ai--;
            } while (ai >= 0 && (Acoeff + Ai)->coeffs[ai] == 0);
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree(A->coeffs + Ai);
                }
            }
        }
        if (Bvalue != 0)
        {
            do {
                bi--;
            } while (bi >= 0 && (Bcoeff + Bi)->coeffs[bi] == 0);
            if (bi < 0)
            {
                Bi++;
                if (Bi < Blen)
                {
                    bi = nmod_poly_degree(B->coeffs + Bi);
                }
            }
        }

        FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    FLINT_ASSERT(nmod_mpolyn_is_canonical(T, ctx));

    *lastdeg_ = lastdeg;
    return;
}


int nmod_mpolyn_interp_crt_2sm_mpolyn(
    slong * lastdeg_,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    nmod_mpolyn_t A,
    nmod_mpolyn_t B,
    slong var,
    nmod_poly_t modulus,
    nmod_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    nmod_poly_t tp;
    nmod_poly_t zero;

    nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    nmod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai, ai;

    nmod_poly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    slong Bi, bi;

    nmod_poly_struct * Fvalue;
    mp_limb_t u, v, Avalue, Bvalue, FvalueA, FvalueB;
    int texp_set, cmp;
    mp_limb_t alpha = nmod_poly_get_coeff_ui(alphapow, 1);

#if WANT_ASSERT
    u = nmod_poly_evaluate_nmod(modulus, alpha);
    u = nmod_mul(u, alpha, ctx->ffinfo->mod);
    u = nmod_mul(u, 2, ctx->ffinfo->mod);
    FLINT_ASSERT(u == 1);
    u = nmod_poly_evaluate_nmod(modulus, ctx->ffinfo->mod.n - alpha);
    u = nmod_mul(u, alpha, ctx->ffinfo->mod);
    u = nmod_mul(u, 2, ctx->ffinfo->mod);
    FLINT_ASSERT(u == 1);
#endif

    FLINT_ASSERT(nmod_mpolyn_is_canonical(A, ctx));
    FLINT_ASSERT(nmod_mpolyn_is_canonical(B, ctx));
    FLINT_ASSERT(nmod_mpolyn_is_canonical(F, ctx));

    nmod_poly_init_mod(tp, ctx->ffinfo->mod);
    nmod_poly_init_mod(zero, ctx->ffinfo->mod);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    nmod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;

    Ti = Fi = Ai = Bi = 0;
    ai = (Ai >= Alen) ? 0 : nmod_poly_degree(A->coeffs + Ai);
    bi = (Bi >= Blen) ? 0 : nmod_poly_degree(B->coeffs + Bi);

    while (Fi < Flen || Ai < Alen || Bi < Blen)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Alen - Ai);
            extra = FLINT_MAX(extra, Blen - Bi);
            nmod_mpolyn_fit_length(T, Ti + extra, ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        FLINT_ASSERT(Fi >= Flen || (Fcoeff + Fi)->length != 0);
        FLINT_ASSERT(Ai >= Alen || (Acoeff + Ai)->coeffs[ai] != 0);
        FLINT_ASSERT(Bi >= Blen || (Bcoeff + Bi)->coeffs[bi] != 0);

        Fvalue = zero;
        texp_set = 0;
        if (Fi < Flen)
        {
            Fvalue = Fcoeff + Fi;
            texp_set = 1;
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);
        }

        Avalue = 0;
        if (Ai < Alen)
        {
            cmp = (!texp_set) ? -1
                     : mpoly_monomial_cmp_nomask_extra(Texp + N*Ti,
                                      Aexp + N*Ai, N, offset, ai << shift);

            if (cmp <= 0)
            {
                Avalue = (Acoeff + Ai)->coeffs[ai];
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                texp_set = 1;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Aexp + N*Ai, N, offset, ai << shift);
            }
        }

        Bvalue = 0;
        if (Bi < Blen)
        {
            cmp = (!texp_set) ? -1 : mpoly_monomial_cmp_nomask_extra(
                             Texp + N*Ti, Bexp + N*Bi, N, offset, bi << shift);

            if (cmp <= 0)
            {
                Bvalue = (Bcoeff + Bi)->coeffs[bi];
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                Avalue = 0;
                texp_set = 1;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Bexp + N*Bi, N, offset, bi << shift);
            }
        }

        FLINT_ASSERT(texp_set);

        _nmod_poly_eval2_pow(&FvalueA, &FvalueB, Fvalue, alphapow, ctx->ffinfo);
        FvalueA = nmod_sub(FvalueA, Avalue, ctx->ffinfo->mod);
        FvalueB = nmod_sub(FvalueB, Bvalue, ctx->ffinfo->mod);
        u = nmod_sub(FvalueB, FvalueA, ctx->ffinfo->mod);
        v = nmod_mul(ctx->ffinfo->mod.n - alpha, nmod_add(FvalueB, FvalueA, ctx->ffinfo->mod), ctx->ffinfo->mod);
        if (u != 0 || v != 0)
        {
            changed = 1;
            nmod_poly_set_coeff_ui(tp, 0, v);
            nmod_poly_set_coeff_ui(tp, 1, u);
            nmod_poly_mul_classical(Tcoeff + Ti, modulus, tp);
            nmod_poly_add(Tcoeff + Ti, Tcoeff + Ti, Fvalue);

        }
        else
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fvalue));
            nmod_poly_set(Tcoeff + Ti, Fvalue);
        }

        Fi += (Fvalue != zero);
        if (Avalue != 0)
        {
            do {
                ai--;
            } while (ai >= 0 && (Acoeff + Ai)->coeffs[ai] == 0);
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = nmod_poly_degree(A->coeffs + Ai);
                }
            }
        }
        if (Bvalue != 0)
        {
            do {
                bi--;
            } while (bi >= 0 && (Bcoeff + Bi)->coeffs[bi] == 0);
            if (bi < 0)
            {
                Bi++;
                if (Bi < Blen)
                {
                    bi = nmod_poly_degree(B->coeffs + Bi);
                }
            }
        }

        FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    if (changed)
    {
        nmod_mpolyn_swap(T, F);
    }

    nmod_poly_clear(tp);
    nmod_poly_clear(zero);

    FLINT_ASSERT(nmod_mpolyn_is_canonical(F, ctx));

    *lastdeg_ = lastdeg;
    return changed;
}


/*****************************************************************************/

/* evaluate A at lastvar = alpha */
void nmod_mpolyn_interp_reduce_sm_mpoly(
    nmod_mpoly_t B,
    nmod_mpolyn_t A,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, N, k;

    FLINT_ASSERT(B->bits == A->bits);

    nmod_mpoly_fit_length(B, A->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    k = 0;
    for (i = 0; i < A->length; i++)
    {
        mpoly_monomial_set(B->exps + N*k, A->exps + N*i, N);
        B->coeffs[k] = nmod_poly_evaluate_nmod(A->coeffs + i, alpha);
        if (B->coeffs[k] != UWORD(0))
        {
            k++;
        }
    }
    B->length = k;
}

void nmod_mpolyun_interp_reduce_sm_mpolyu(
    nmod_mpolyu_t B,
    nmod_mpolyun_t A,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, k;

    FLINT_ASSERT(B->bits == A->bits);

    nmod_mpolyu_fit_length(B, A->length, ctx);
    k = 0;
    for (i = 0; i < A->length; i++)
    {
        B->exps[k] = A->exps[i];
        nmod_mpolyn_interp_reduce_sm_mpoly(B->coeffs + k, A->coeffs + i, alpha, ctx);
        k += !nmod_mpoly_is_zero(B->coeffs + k, ctx);
    }
    B->length = k;
}


void nmod_mpolyn_interp_lift_sm_mpoly(
    nmod_mpolyn_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, N;
    nmod_poly_struct * Acoeff;
    mp_limb_t * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Blen;

    FLINT_ASSERT(A->bits == B->bits);

    Blen = B->length;
    nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < Blen; i++)
    {
        nmod_poly_zero(Acoeff + i);
        nmod_poly_set_coeff_ui(Acoeff + i, 0, Bcoeff[i]);
        mpoly_monomial_set(Aexp + N*i, Bexp + N*i, N);
    }
    A->length = Blen;
}

void nmod_mpolyun_interp_lift_sm_mpolyu(
    nmod_mpolyun_t A,
    const nmod_mpolyu_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(A->bits == B->bits);
    nmod_mpolyun_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        nmod_mpolyn_interp_lift_sm_mpoly(A->coeffs + i, B->coeffs + i, ctx);
        FLINT_ASSERT((A->coeffs + i)->bits == B->bits);
    }
    A->length = B->length;
}

/*
    F = F + modulus*(A - F(alpha))
    no assumptions about matching monomials
*/
int nmod_mpolyn_interp_crt_sm_mpoly(
    slong * lastdeg,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    nmod_mpoly_t A,
    nmod_poly_t modulus,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    slong N;
    mp_limb_t v;
    flint_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    mp_limb_t * Acoeff = A->coeffs;
    nmod_poly_struct * Fcoeff = F->coeffs;
    nmod_poly_struct * Tcoeff;
    nmod_poly_t tp;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    nmod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen
                        || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + i) < nmod_poly_degree(modulus));

            /* F term ok, A term missing */
            v = nmod_poly_evaluate_nmod(Fcoeff + i, alpha);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_sub(Tcoeff + k, Fcoeff + i, tp);
            } else {
                nmod_poly_set(Tcoeff + k, Fcoeff + i);                
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + k));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (Acoeff[j] != UWORD(0))
            {
                changed = 1;
                nmod_poly_zero(Tcoeff + k);
                nmod_poly_scalar_mul_nmod(Tcoeff + k, modulus, Acoeff[j]);
                lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + i) < nmod_poly_degree(modulus));

            /* F term ok, A term ok */
            v = nmod_poly_evaluate_nmod(Fcoeff + i, alpha);
            v = nmod_sub(Acoeff[j], v, ctx->ffinfo->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                nmod_poly_scalar_mul_nmod(tp, modulus, v);
                nmod_poly_add(Tcoeff + k, Fcoeff + i, tp);
            } else {
                nmod_poly_set(Tcoeff + k, Fcoeff + i);                
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], nmod_poly_degree(Tcoeff + k));
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + k));
            k++;
            i++;
            j++;
        }
        else
        {
            FLINT_ASSERT(0);
        }
    }

    nmod_mpolyn_set_length(T, k, ctx);

    if (changed)
    {
        nmod_mpolyn_swap(T, F);
    }

    nmod_poly_clear(tp);
    return changed;
}

int nmod_mpolyun_interp_crt_sm_mpolyu(
    slong * lastdeg,
    nmod_mpolyun_t F,
    nmod_mpolyun_t T,
    nmod_mpolyu_t A,
    nmod_poly_t modulus,
    mp_limb_t alpha,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    nmod_mpolyn_t S;
    nmod_mpolyn_struct * Tcoeff;
    nmod_mpolyn_struct * Fcoeff;
    nmod_mpoly_struct  * Acoeff;
    nmod_mpoly_t zero;

    lastdeg[0] = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    nmod_mpolyn_init(S, F->bits, ctx);

    Flen = F->length;
    Alen = A->length;
    nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    nmod_mpoly_init(zero, ctx);
    nmod_mpoly_fit_bits(zero, A->bits, ctx);
    zero->bits = A->bits;

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= nmod_mpolyn_interp_crt_sm_mpoly(lastdeg, Tcoeff + k,
                                                 S, zero, modulus, alpha, ctx);
            Texp[k] = Fexp[i];
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || Aexp[j] > Fexp[i]))
        {
            /* F term missing, A term ok */
            nmod_mpolyn_zero(Tcoeff + k, ctx);
            changed |= nmod_mpolyn_interp_crt_sm_mpoly(lastdeg, Tcoeff + k,
                                           S, Acoeff + j, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            k++;
            j++;
        }
        else if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= nmod_mpolyn_interp_crt_sm_mpoly(lastdeg, Tcoeff + k,
                                           S, Acoeff + j, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            FLINT_ASSERT(!nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            i++;
            j++;
        }
        else
        {
            FLINT_ASSERT(0);
        }
    }
    T->length = k;

    if (changed)
    {
        nmod_mpolyun_swap(T, F);
    }

    nmod_mpolyn_clear(S, ctx);
    nmod_mpoly_clear(zero, ctx);
    return changed;    
}

