/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "n_poly.h"
#include "mpoly.h"
#include "nmod_mpoly.h"

/*
    interp_reduce: map from Fp[x] to Fp[x]/poly(x)
    interp_lift:   map from Fp[x]/poly(x) to Fp[x]
    interp_crt:    update element of Fp[x] with a new image in Fp[x]/poly(x)
    interp_mcrt:   same as interp_crt, but monomial match, thus easier
*/


void _nmod_poly_eval2_pow(
    ulong * vp,
    ulong * vm,
    n_poly_t P,
    n_poly_t alphapow,
    nmod_t fctx)
{
    ulong * Pcoeffs = P->coeffs;
    slong Plen = P->length;
    ulong * alpha_powers = alphapow->coeffs;
    ulong p1, p0, a0, a1, a2, q1, q0, b0, b1, b2;
    slong k;

    a0 = a1 = a2 = UWORD(0);
    b0 = b1 = b2 = UWORD(0);

    if (Plen > alphapow->length)
    {
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        n_poly_fit_length(alphapow, Plen);
        for (k = oldlength; k < Plen; k++)
        {
            alphapow->coeffs[k] = nmod_mul(alphapow->coeffs[k - 1],
                                           alphapow->coeffs[1], fctx);
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

    NMOD_RED3(p0, a2, a1, a0, fctx);
    NMOD_RED3(q0, b2, b1, b0, fctx);

    *vp = nmod_add(p0, q0, fctx);
    *vm = nmod_sub(p0, q0, fctx);
}

/*****************************************************************************/

void nmod_mpolyn_interp_reduce_2sm_poly(
    n_poly_t E,
    n_poly_t F,
    const nmod_mpolyn_t A,
    n_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    ulong u, v;
    slong Ai, Alen, k;
    n_poly_struct * Acoeff;
    ulong * Aexp;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    n_poly_zero(E);
    n_poly_zero(F);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        _nmod_poly_eval2_pow(&u, &v, Acoeff + Ai, alphapow, ctx->mod);
        k = (Aexp + N*Ai)[off] >> shift;
        n_poly_set_coeff(E, k, u);
        n_poly_set_coeff(F, k, v);
    }
}

void nmod_mpolyn_interp_lift_2sm_poly(
    slong * lastdeg_,
    nmod_mpolyn_t F,
    const n_poly_t A,
    const n_poly_t B,
    ulong alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong lastdeg = -WORD(1);
    ulong u, v, d0, d1, Avalue, Bvalue;
    slong Fi, Aexp, Bexp;
    ulong * Acoeff = A->coeffs;
    ulong * Bcoeff = B->coeffs;
    n_poly_struct * Fcoeff;
    ulong * Fexp;
    slong e;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    Aexp = n_poly_degree(A);
    Bexp = n_poly_degree(B);

    nmod_mpolyn_fit_length(F, FLINT_MAX(Aexp, Bexp) + 1, ctx);
    Fcoeff = F->coeffs;
    Fexp = F->exps;

    d0 = n_invmod(UWORD(2), ctx->mod.n);
    d1 = n_invmod(nmod_add(alpha, alpha, ctx->mod), ctx->mod.n);

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
        u = nmod_add(Avalue, Bvalue, ctx->mod);
        v = nmod_sub(Avalue, Bvalue, ctx->mod);
        u = nmod_mul(u, d0, ctx->mod);
        v = nmod_mul(v, d1, ctx->mod);

        FLINT_ASSERT(Fi < F->alloc);
        mpoly_monomial_zero(Fexp + N*Fi, N);
        (Fexp + N*Fi)[off] = e << shift;

        FLINT_ASSERT(u != 0 || v != 0);
        n_poly_fit_length(Fcoeff + Fi, 2);
        (Fcoeff + Fi)->coeffs[0] = u;
        (Fcoeff + Fi)->coeffs[1] = v;
        (Fcoeff + Fi)->length = 1 + (v != 0);
        lastdeg = FLINT_MAX(lastdeg, n_poly_degree(Fcoeff + Fi));
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
    const n_poly_t A,
    const n_poly_t B,
    const n_poly_t modulus,
    n_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0, Finc;
    ulong alpha = n_poly_get_coeff(alphapow, 1);
    slong lastdeg = -WORD(1);
    ulong u, v, FvalueA, FvalueB;
    slong Fi, Toff, Aexp, Bexp, e, fexp;
    ulong * Acoeff = A->coeffs;
    ulong * Bcoeff = B->coeffs;
    slong Flen = F->length;
    n_poly_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    n_poly_struct * Tcoeff;
    ulong * Texp;
    slong N, off, shift;

    FLINT_ASSERT(T->bits == F->bits);

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    Fi = 0;
    Aexp = n_poly_degree(A);
    Bexp = n_poly_degree(B);

    nmod_mpolyn_fit_length(T, Flen + FLINT_MAX(Aexp, Bexp) + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Toff = 0;

#if FLINT_WANT_ASSERT
    u = n_poly_mod_evaluate_nmod(modulus, alpha, ctx->mod);
    u = nmod_mul(u, alpha, ctx->mod);
    u = nmod_mul(u, 2, ctx->mod);
    FLINT_ASSERT(u == 1);
    u = n_poly_mod_evaluate_nmod(modulus, ctx->mod.n - alpha, ctx->mod);
    u = nmod_mul(u, alpha, ctx->mod);
    u = nmod_mul(u, 2, ctx->mod);
    FLINT_ASSERT(u == 1);
#endif

    while (Fi < Flen || Aexp >= 0 || Bexp >= 0)
    {
        FLINT_ASSERT(Toff < T->alloc);

        fexp = e = -WORD(1);
        if (Fi < Flen)
        {
            fexp = e = (Fexp + N*Fi)[off]>>shift;
            FLINT_ASSERT(!n_poly_is_zero(Fcoeff + Fi));
            FLINT_ASSERT(n_poly_degree(Fcoeff + Fi) < n_poly_degree(modulus));
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
            _nmod_poly_eval2_pow(&FvalueA, &FvalueB, Fcoeff + Fi, alphapow, ctx->mod);
        }

        if (e == Aexp)
        {
            FvalueA = nmod_sub(FvalueA, Acoeff[Aexp], ctx->mod);
        }
        if (e == Bexp)
        {
            FvalueB = nmod_sub(FvalueB, Bcoeff[Bexp], ctx->mod);
        }

        u = nmod_sub(FvalueB, FvalueA, ctx->mod);
        v = nmod_mul(ctx->mod.n - alpha, nmod_add(FvalueB, FvalueA, ctx->mod), ctx->mod);

        if (u != 0 || v != 0)
        {
            changed = 1;
            if (u != 0)
            {
                _n_poly_mod_scalar_mul_nmod(Tcoeff + Toff, modulus, u, ctx->mod);
                n_poly_shift_left(Tcoeff + Toff, Tcoeff + Toff, 1);
                _nmod_vec_scalar_addmul_nmod((Tcoeff + Toff)->coeffs,
                                modulus->coeffs, modulus->length, v, ctx->mod);
            }
            else
            {
                _n_poly_mod_scalar_mul_nmod(Tcoeff + Toff, modulus, v, ctx->mod);
            }

            if (Finc)
            {
                n_poly_mod_add(Tcoeff + Toff, Tcoeff + Toff, Fcoeff + Fi, ctx->mod);
            }
        }
        else
        {
            FLINT_ASSERT(Finc == 1);
            n_poly_set(Tcoeff + Toff, Fcoeff + Fi);
        }

        FLINT_ASSERT(!n_poly_is_zero(Tcoeff + Toff));
        lastdeg = FLINT_MAX(lastdeg, n_poly_degree(Tcoeff + Toff));
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


void nmod_mpolyn_interp_lift_sm_bpoly(
    nmod_mpolyn_t F,
    n_bpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    slong i, j, Fi;
    slong off0, shift0, off1, shift1;

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

    Fi = 0;
    for (i = A->length - 1; i >= 0; i--)
    {
        n_poly_struct * Ai = A->coeffs + i;
        for (j = Ai->length - 1; j >= 0; j--)
        {
            if (Ai->coeffs[j] == 0)
                continue;

            nmod_mpolyn_fit_length(F, Fi + 1, ctx);
            mpoly_monomial_zero(F->exps + N*Fi, N);
            (F->exps + N*Fi)[off0] += (i << shift0);
            (F->exps + N*Fi)[off1] += (j << shift1);
            n_poly_set_ui(F->coeffs + Fi, Ai->coeffs[j]);
            Fi++;
        }
    }

    F->length = Fi;
}

int nmod_mpolyn_interp_crt_sm_bpoly(
    slong * lastdeg,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    n_bpoly_t A,
    n_poly_t modulus,
    n_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    const dot_params_t params = _nmod_vec_dot_params(modulus->length, ctx->mod);
    slong N = mpoly_words_per_exp(F->bits, ctx->minfo);
    slong off0, shift0, off1, shift1;
    n_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    n_poly_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    n_poly_struct * Tcoeffs = T->coeffs;
    ulong v;
    ulong Fexpi, mask;

    mask = (-UWORD(1)) >> (FLINT_BITS - F->bits);
    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

    FLINT_ASSERT(T->bits == F->bits);

    *lastdeg = -1;

    Ti = Fi = 0;
    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(Acoeffs + Ai);

    while (Fi < Flen || Ai >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Ai);
            nmod_mpolyn_fit_length(T, Ti + extra + 1, ctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        if (Fi < Flen)
            Fexpi = pack_exp2(((Fexps + N*Fi)[off0]>>shift0)&mask,
                              ((Fexps + N*Fi)[off1]>>shift1)&mask);
        else
            Fexpi = 0;

        if (Fi < Flen && Ai >= 0 && Fexpi == pack_exp2(Ai, ai))
        {
            /* F term ok, A term ok */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            v = _n_poly_eval_pow(Fcoeffs + Fi, alphapow, params, ctx->mod);
            v = nmod_sub(Acoeffs[Ai].coeffs[ai], v, ctx->mod);
            if (v != 0)
            {
                changed = 1;
                n_poly_mod_scalar_addmul_nmod(Tcoeffs + Ti, Fcoeffs + Fi,
                                                         modulus, v, ctx->mod);
            }
            else
            {
                n_poly_set(Tcoeffs + Ti, Fcoeffs + Fi);
            }

            Fi++;
            do {
                ai--;
            } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Ai >= 0 && (Fi >= Flen || Fexpi < pack_exp2(Ai, ai)))
        {
            /* F term missing, A term ok */

            mpoly_monomial_zero(Texps + N*Ti, N);
            (Texps + N*Ti)[off0] += (Ai << shift0);
            (Texps + N*Ti)[off1] += (ai << shift1);

            changed = 1;
            _n_poly_mod_scalar_mul_nmod(Tcoeffs + Ti, modulus, Acoeffs[Ai].coeffs[ai], ctx->mod);

            do {
                ai--;
            } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
        else
        {
            FLINT_ASSERT(Fi < Flen && (Ai < 0 || Fexpi > pack_exp2(Ai, ai)));
            /* F term ok, Aterm missing */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            v = _n_poly_eval_pow(Fcoeffs + Fi, alphapow, params, ctx->mod);
            if (v != 0)
            {
                changed = 1;
                v = nmod_neg(v, ctx->mod);
                n_poly_mod_scalar_addmul_nmod(Tcoeffs + Ti, Fcoeffs + Fi,
                                                         modulus, v, ctx->mod);
            }
            else
            {
                n_poly_set(Tcoeffs + Ti, Fcoeffs + Fi);
            }

            Fi++;
        }

        FLINT_ASSERT(!n_poly_is_zero(Tcoeffs + Ti));
        *lastdeg = FLINT_MAX(*lastdeg, n_poly_degree(Tcoeffs + Ti));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        nmod_mpolyn_swap(T, F);

    return changed;
}

/*****************************************************************************/

/* The following functions are currently undocumented, unused. */

#if 0

/*
    E = A(x_var = alpha)
    A is in Fp[x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in Fp[x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyn_interp_reduce_sm_mpolyn(
    nmod_mpolyn_t E,
    nmod_mpolyn_t A,
    slong var,
    ulong alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    ulong v;
    n_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    n_poly_struct * Ecoeff;
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
        v = n_poly_mod_evaluate_nmod(Acoeff + Ai, alpha, ctx->mod);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (v == 0)
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            n_poly_set_coeff(Ecoeff + Ei - 1, k, v);
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
            n_poly_zero(Ecoeff + Ei);
            n_poly_set_coeff(Ecoeff + Ei, k, v);
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
    n_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;
    n_poly_struct * Acoeff;
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
                n_poly_zero(Acoeff + Ai);
                n_poly_set_coeff(Acoeff + Ai, 0, (Bcoeff + Bi)->coeffs[vi]);
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
    n_poly_t modulus,
    ulong alpha,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    slong vi;
    ulong v;
    n_poly_t tp;
    n_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;
    n_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;
    n_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    n_poly_init(tp);

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
        vi = n_poly_degree(Acoeff + Ai);
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
            v = n_poly_mod_evaluate_nmod(Fcoeff + Fi, alpha, ctx->mod);
            v = nmod_sub((Acoeff + Ai)->coeffs[vi], v, ctx->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                _n_poly_mod_scalar_mul_nmod(tp, modulus, v, ctx->mod);
                n_poly_mod_add(Tcoeff + Ti, Fcoeff + Fi, tp, ctx->mod);
            }
            else
            {
                n_poly_set(Tcoeff + Ti, Fcoeff + Fi);
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
                    vi = n_poly_degree(Acoeff + Ai);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            v = n_poly_mod_evaluate_nmod(Fcoeff + Fi, alpha, ctx->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                _n_poly_mod_scalar_mul_nmod(tp, modulus, v, ctx->mod);
                n_poly_mod_sub(Tcoeff + Ti, Fcoeff + Fi, tp, ctx->mod);
            }
            else
            {
                n_poly_set(Tcoeff + Ti, Fcoeff + Fi);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            changed = 1;
            _n_poly_mod_scalar_mul_nmod(Tcoeff + Ti, modulus, (Acoeff + Ai)->coeffs[vi], ctx->mod);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && (Acoeff + Ai)->coeffs[vi] == 0);
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = n_poly_degree(Acoeff + Ai);
                }
            }
        }

        FLINT_ASSERT(!n_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, n_poly_degree(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    n_poly_clear(tp);

    *lastdeg_ = FLINT_MAX(*lastdeg_, lastdeg);
    return changed;
}

#endif


/****************************************************************************/

void nmod_mpolyn_interp_reduce_2sm_mpolyn(
    nmod_mpolyn_t E,
    nmod_mpolyn_t F,
    nmod_mpolyn_t A,
    slong var,
    n_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    ulong e, f;
    n_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    n_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;
    n_poly_struct * Fcoeff;
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
        _nmod_poly_eval2_pow(&e, &f, Acoeff + Ai, alphapow, ctx->mod);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;

        if (e != 0)
        {
            if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
            {
                /* append to previous */
                n_poly_set_coeff(Ecoeff + Ei - 1, k, e);
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
                n_poly_zero(Ecoeff + Ei);
                n_poly_set_coeff(Ecoeff + Ei, k, e);
                Ei++;
            }
        }

        if (f != 0)
        {
            if (Fi > 0 && mpoly_monomial_equal_extra(Fexp + N*(Fi - 1), Aexp + N*Ai, N, offset, -(k << shift)))
            {
                /* append to previous */
                n_poly_set_coeff(Fcoeff + Fi - 1, k, f);
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
                n_poly_zero(Fcoeff + Fi);
                n_poly_set_coeff(Fcoeff + Fi, k, f);
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
    ulong alpha,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    n_poly_t tp;
    n_poly_t zero;
    n_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;
    n_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai, ai;
    n_poly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    slong Bi, bi;
    ulong u, v, Avalue, Bvalue, FvalueA, FvalueB;
    int cmp;
    ulong d0 = n_invmod(alpha + alpha, ctx->mod.n);

    n_poly_init(tp);
    n_poly_init(zero);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    nmod_mpolyn_fit_length(T, FLINT_MAX(Alen, Blen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;

    Ti = Ai = Bi = 0;
    ai = (Ai >= Alen) ? 0 : n_poly_degree(Acoeff + Ai);
    bi = (Bi >= Blen) ? 0 : n_poly_degree(Bcoeff + Bi);

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

        FvalueA = nmod_neg(Avalue, ctx->mod);
        FvalueB = nmod_neg(Bvalue, ctx->mod);
        u = nmod_sub(FvalueB, FvalueA, ctx->mod);
        v = nmod_mul(ctx->mod.n - alpha, nmod_add(FvalueB, FvalueA, ctx->mod), ctx->mod);

        FLINT_ASSERT(u != 0 || v != 0);
        n_poly_zero(Tcoeff + Ti);
        u = nmod_mul(u, d0, ctx->mod);
        v = nmod_mul(v, d0, ctx->mod);
        n_poly_set_coeff(Tcoeff + Ti, 0, v);
        n_poly_set_coeff(Tcoeff + Ti, 1, u);

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
                    ai = n_poly_degree(Acoeff + Ai);
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
                    bi = n_poly_degree(Bcoeff + Bi);
                }
            }
        }

        FLINT_ASSERT(!n_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, n_poly_degree(Tcoeff + Ti));
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
    n_poly_t modulus,
    n_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    n_poly_t tp;
    n_poly_t zero;
    n_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;
    n_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;
    n_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai, ai;
    n_poly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    slong Bi, bi;
    n_poly_struct * Fvalue;
    ulong u, v, Avalue, Bvalue, FvalueA, FvalueB;
    int texp_set, cmp;
    ulong alpha = n_poly_get_coeff(alphapow, 1);

#if FLINT_WANT_ASSERT
    u = n_poly_mod_evaluate_nmod(modulus, alpha, ctx->mod);
    u = nmod_mul(u, alpha, ctx->mod);
    u = nmod_mul(u, 2, ctx->mod);
    FLINT_ASSERT(u == 1);
    u = n_poly_mod_evaluate_nmod(modulus, ctx->mod.n - alpha, ctx->mod);
    u = nmod_mul(u, alpha, ctx->mod);
    u = nmod_mul(u, 2, ctx->mod);
    FLINT_ASSERT(u == 1);
#endif

    FLINT_ASSERT(nmod_mpolyn_is_canonical(A, ctx));
    FLINT_ASSERT(nmod_mpolyn_is_canonical(B, ctx));
    FLINT_ASSERT(nmod_mpolyn_is_canonical(F, ctx));

    n_poly_init(tp);
    n_poly_init(zero);

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
    ai = (Ai >= Alen) ? 0 : n_poly_degree(Acoeff + Ai);
    bi = (Bi >= Blen) ? 0 : n_poly_degree(Bcoeff + Bi);

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

        _nmod_poly_eval2_pow(&FvalueA, &FvalueB, Fvalue, alphapow, ctx->mod);
        FvalueA = nmod_sub(FvalueA, Avalue, ctx->mod);
        FvalueB = nmod_sub(FvalueB, Bvalue, ctx->mod);
        u = nmod_sub(FvalueB, FvalueA, ctx->mod);
        v = nmod_mul(ctx->mod.n - alpha, nmod_add(FvalueB, FvalueA, ctx->mod), ctx->mod);
        if (u != 0 || v != 0)
        {
            changed = 1;
            n_poly_set_coeff(tp, 0, v);
            n_poly_set_coeff(tp, 1, u);
            n_poly_mod_mul(Tcoeff + Ti, modulus, tp, ctx->mod);
            n_poly_mod_add(Tcoeff + Ti, Tcoeff + Ti, Fvalue, ctx->mod);

        }
        else
        {
            FLINT_ASSERT(!n_poly_is_zero(Fvalue));
            n_poly_set(Tcoeff + Ti, Fvalue);
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
                    ai = n_poly_degree(Acoeff + Ai);
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
                    bi = n_poly_degree(Bcoeff + Bi);
                }
            }
        }

        FLINT_ASSERT(!n_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, n_poly_degree(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    if (changed)
    {
        nmod_mpolyn_swap(T, F);
    }

    n_poly_clear(tp);
    n_poly_clear(zero);

    FLINT_ASSERT(nmod_mpolyn_is_canonical(F, ctx));

    *lastdeg_ = lastdeg;
    return changed;
}


/*****************************************************************************/

/* evaluate A at lastvar = alpha */
void nmod_mpolyn_interp_reduce_sm_mpoly(
    nmod_mpoly_t B,
    nmod_mpolyn_t A,
    ulong alpha,
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
        B->coeffs[k] = n_poly_mod_evaluate_nmod(A->coeffs + i, alpha, ctx->mod);
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
    ulong alpha,
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
    n_poly_struct * Acoeff;
    ulong * Bcoeff;
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
        n_poly_zero(Acoeff + i);
        n_poly_set_coeff(Acoeff + i, 0, Bcoeff[i]);
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
    n_poly_t modulus,
    ulong alpha,
    const nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    slong N;
    ulong v;
    flint_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    ulong * Acoeff = A->coeffs;
    n_poly_struct * Fcoeff = F->coeffs;
    n_poly_struct * Tcoeff;
    n_poly_t tp;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    n_poly_init(tp);

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
            FLINT_ASSERT(!n_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(n_poly_degree(Fcoeff + i) < n_poly_degree(modulus));

            /* F term ok, A term missing */
            v = n_poly_mod_evaluate_nmod(Fcoeff + i, alpha, ctx->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                _n_poly_mod_scalar_mul_nmod(tp, modulus, v, ctx->mod);
                n_poly_mod_sub(Tcoeff + k, Fcoeff + i, tp, ctx->mod);
            } else {
                n_poly_set(Tcoeff + k, Fcoeff + i);
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], n_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            FLINT_ASSERT(!n_poly_is_zero(Tcoeff + k));
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
                n_poly_zero(Tcoeff + k);
                _n_poly_mod_scalar_mul_nmod(Tcoeff + k, modulus, Acoeff[j], ctx->mod);
                lastdeg[0] = FLINT_MAX(lastdeg[0], n_poly_degree(Tcoeff + k));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!n_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(n_poly_degree(Fcoeff + i) < n_poly_degree(modulus));

            /* F term ok, A term ok */
            v = n_poly_mod_evaluate_nmod(Fcoeff + i, alpha, ctx->mod);
            v = nmod_sub(Acoeff[j], v, ctx->mod);
            if (v != UWORD(0))
            {
                changed = 1;
                _n_poly_mod_scalar_mul_nmod(tp, modulus, v, ctx->mod);
                n_poly_mod_add(Tcoeff + k, Fcoeff + i, tp, ctx->mod);
            } else {
                n_poly_set(Tcoeff + k, Fcoeff + i);
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], n_poly_degree(Tcoeff + k));
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            FLINT_ASSERT(!n_poly_is_zero(Tcoeff + k));
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

    n_poly_clear(tp);
    return changed;
}

int nmod_mpolyun_interp_crt_sm_mpolyu(
    slong * lastdeg,
    nmod_mpolyun_t F,
    nmod_mpolyun_t T,
    nmod_mpolyu_t A,
    n_poly_t modulus,
    ulong alpha,
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

    nmod_mpoly_init3(zero, 0, A->bits, ctx);

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


/*
    F = F + modulus*(A - F(alpha))
    monomials assumed to match
*/
int nmod_mpolyn_interp_mcrt_sm_mpoly(
    slong * lastdeg_,
    nmod_mpolyn_t F,
    const nmod_mpoly_t A,
    const n_poly_t modulus,
    n_poly_t alphapow,
    const nmod_mpoly_ctx_t ctx)
{
    slong lastdeg = -1;
    int changed = 0;
    slong i;
    ulong v;
    ulong * Acoeff = A->coeffs;
    slong Flen = F->length;

    FLINT_ASSERT(Flen == A->length);

    for (i = 0; i < Flen; i++)
    {
        /* F term ok, A term ok */
        v = n_poly_mod_eval_pow(F->coeffs + i, alphapow, ctx->mod);
        v = nmod_sub(Acoeff[i], v, ctx->mod);
        if (v != 0)
        {
            changed = 1;
            n_poly_mod_scalar_addmul_nmod(F->coeffs + i,
                                          F->coeffs + i, modulus, v, ctx->mod);
        }
        lastdeg = FLINT_MAX(lastdeg, n_poly_degree(F->coeffs + i));
    }

    *lastdeg_ = lastdeg;
    return changed;
}
