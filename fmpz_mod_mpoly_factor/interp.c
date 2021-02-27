/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"


/*
    E = A(v = alpha)
    A is in R[v][x]
    E is in R[x]
*/
void fmpz_mod_polyu1n_intp_reduce_sm_poly(
    fmpz_mod_poly_t E,
    const fmpz_mod_polyun_t A,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_t v;
    slong Ai;

    fmpz_init(v);
    fmpz_mod_poly_zero(E, ctx);
    for (Ai = 0; Ai < A->length; Ai++)
    {
        fmpz_mod_poly_evaluate_fmpz(v, A->coeffs + Ai, alpha, ctx);
        fmpz_mod_poly_set_coeff_fmpz(E, A->exps[Ai], v, ctx);
    }
    fmpz_clear(v);
}

/*
    A = B
    A, B are in R[X]
*/
void fmpz_mod_polyu1n_intp_lift_sm_poly(
    fmpz_mod_polyun_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong Bi;
    slong Blen = B->length;
    fmpz * Bcoeff = B->coeffs;
    slong Ai;

    fmpz_mod_polyun_fit_length(A, Blen, ctx);

    Ai = 0;
    for (Bi = Blen - 1; Bi >= 0; Bi--)
    {
        if (fmpz_is_zero(Bcoeff + Bi))
            continue;

        FLINT_ASSERT(Ai < A->alloc);

        fmpz_mod_poly_set_fmpz(A->coeffs + Ai, Bcoeff + Bi, ctx);
        A->exps[Ai] = Bi;
        Ai++;
    }
    A->length = Ai;
}


/*
    F = F + modulus*(A - F(v = alpha))
    no assumptions about matching monomials
    F is in Fp[X][v]
    A is in Fp[X]
    it is expected that modulus(alpha) == 1
*/
int fmpz_mod_polyu1n_intp_crt_sm_poly(
    slong * lastdeg,
    fmpz_mod_polyun_t F,
    fmpz_mod_polyun_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t modulus,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx)
{
    int changed = 0;
    slong lastlen = 0;
    fmpz_t v;
    slong Fi, Ti, Ai;
    fmpz * Acoeffs = A->coeffs;
    slong Flen = F->length;
    fmpz_mod_poly_struct * Fcoeffs = F->coeffs;
    ulong * Fexps = F->exps;
    fmpz_mod_poly_struct * Tcoeffs;
    ulong * Texps;

    Fi = 0;
    Ai = fmpz_mod_poly_degree(A, ctx);

    fmpz_init(v);

    fmpz_mod_polyun_fit_length(T, Flen + Ai + 1, ctx);
    Tcoeffs = T->coeffs;
    Texps = T->exps;
    Ti = 0;

    while (Fi < Flen || Ai >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Fcoeffs + Fi, ctx));
            FLINT_ASSERT(fmpz_mod_poly_degree(Fcoeffs + Fi, ctx) <
                                           fmpz_mod_poly_degree(modulus, ctx));
        }

        if (Ai >= 0)
        {
            FLINT_ASSERT(!fmpz_is_zero(Acoeffs + Ai));
        }

        if (Fi < Flen && Ai >= 0 && Fexps[Fi] == Ai)
        {
            /* F term ok, A term ok */
            fmpz_mod_poly_evaluate_fmpz(v, Fcoeffs + Fi, alpha, ctx);
            fmpz_mod_sub(v, Acoeffs + Ai, v, ctx);
            changed |= !fmpz_is_zero(v);
            fmpz_mod_poly_scalar_addmul_fmpz_mod(Tcoeffs + Ti,
                                                Fcoeffs + Fi, modulus, v, ctx);
            Texps[Ti] = Ai;
            Fi++;
            do {
                Ai--;
            } while (Ai >= 0 && fmpz_is_zero(Acoeffs + Ai));
        }
        else if (Fi < Flen && (Ai < 0 || Fexps[Fi] > Ai))
        {
            /* F term ok, A term missing */
            fmpz_mod_poly_evaluate_fmpz(v, Fcoeffs + Fi, alpha, ctx);
            fmpz_mod_neg(v, v, ctx);
            changed |= !fmpz_is_zero(v);
            fmpz_mod_poly_scalar_addmul_fmpz_mod(Tcoeffs + Ti,
                                                Fcoeffs + Fi, modulus, v, ctx);
            Texps[Ti] = Fexps[Fi];
            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai >= 0 && (Fi >= Flen || Fexps[Fi] < Ai));

            /* F term missing, A term ok */
            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz(Tcoeffs + Ti,
                                                   modulus, Acoeffs + Ai, ctx);
            Texps[Ti] = Ai;
            do {
                Ai--;
            } while (Ai >= 0 && fmpz_is_zero(Acoeffs + Ai));
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeffs + Ti, ctx));
        lastlen = FLINT_MAX(lastlen, Tcoeffs[Ti].length);

        Ti++;
    }
    T->length = Ti;

    fmpz_clear(v);

    if (changed)
        fmpz_mod_polyun_swap(T, F);

    *lastdeg = lastlen - 1;
    return changed;
}


void fmpz_mod_polyu1n_interp_reduce_2sm_poly(
    fmpz_mod_poly_t E,
    fmpz_mod_poly_t F,
    const fmpz_mod_polyun_t A,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_t u, v;

    fmpz_init(u);
    fmpz_init(v);

    fmpz_mod_poly_zero(E, ctx);
    fmpz_mod_poly_zero(F, ctx);
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_eval2_pow(u, v, A->coeffs + i, alphapow, ctx);
        fmpz_mod_poly_set_coeff_fmpz(E, A->exps[i], u, ctx);
        fmpz_mod_poly_set_coeff_fmpz(F, A->exps[i], v, ctx);
    }

    fmpz_clear(u);
    fmpz_clear(v);
}

void fmpz_mod_polyu1n_interp_lift_2sm_poly(
    slong * lastdeg,
    fmpz_mod_polyun_t F,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx)
{
    slong lastlen = 0;
    fmpz_t u, v, d0, d1, Avalue, Bvalue;
    slong Fi, Aexp, Bexp;
    const fmpz * Acoeff = A->coeffs;
    const fmpz * Bcoeff = B->coeffs;
    fmpz_mod_poly_struct * Fcoeffs;
    ulong * Fexps;
    slong e;

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(d0);
    fmpz_init(d1);
    fmpz_init(Avalue);
    fmpz_init(Bvalue);

    Aexp = _fmpz_mod_poly_degree(A);
    Bexp = _fmpz_mod_poly_degree(B);

    fmpz_mod_polyun_fit_length(F, FLINT_MAX(Aexp, Bexp) + 1, ctx);
    Fcoeffs = F->coeffs;
    Fexps = F->exps;

    fmpz_set_ui(d0, 2);
    fmpz_mod_inv(d0, d0, ctx);
    fmpz_mod_add(d1, alpha, alpha, ctx);
    fmpz_mod_inv(d1, d1, ctx);

    Fi = 0;
    while (Aexp >= 0 || Bexp >= 0)
    {
        e = Aexp;
        fmpz_zero(Avalue);
        fmpz_zero(Bvalue);
        if (Aexp == Bexp)
        {
            fmpz_set(Avalue, Acoeff + Aexp);
            fmpz_set(Bvalue, Bcoeff + Bexp);
        }
        else if (Aexp > Bexp)
        {
            fmpz_set(Avalue, Acoeff + Aexp);
        }
        else
        {
            FLINT_ASSERT(Bexp > Aexp);
            e = Bexp;
            fmpz_set(Bvalue, Bcoeff + Bexp);
        }
        FLINT_ASSERT(!fmpz_is_zero(Avalue) || !fmpz_is_zero(Bvalue));
        fmpz_mod_add(u, Avalue, Bvalue, ctx);
        fmpz_mod_sub(v, Avalue, Bvalue, ctx);
        fmpz_mod_mul(u, u, d0, ctx);
        fmpz_mod_mul(v, v, d1, ctx);

        FLINT_ASSERT(Fi < F->alloc);

        Fexps[Fi] = e;

        FLINT_ASSERT(!fmpz_is_zero(u) || !fmpz_is_zero(v));
        fmpz_mod_poly_fit_length(Fcoeffs + Fi, 2, ctx);
        fmpz_set(Fcoeffs[Fi].coeffs + 0, u);
        fmpz_set(Fcoeffs[Fi].coeffs + 1, v);
        Fcoeffs[Fi].length = 1 + !fmpz_is_zero(v);
        lastlen = FLINT_MAX(lastlen, Fcoeffs[Fi].length);
        Fi++;

        if (e == Aexp)
        {
            do {
                Aexp--;
            } while (Aexp >= 0 && fmpz_is_zero(Acoeff + Aexp));
        }
        if (e == Bexp)
        {
            do {
                Bexp--;
            } while (Bexp >= 0 && fmpz_is_zero(Bcoeff + Bexp));
        }
    }
    F->length = Fi;

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(d0);
    fmpz_clear(d1);
    fmpz_clear(Avalue);
    fmpz_clear(Bvalue);


    *lastdeg = lastlen - 1;
    return;
}

int fmpz_mod_polyu1n_interp_crt_2sm_poly(
    slong * lastdeg,
    fmpz_mod_polyun_t F,
    fmpz_mod_polyun_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx)
{
    int changed = 0, Finc;
    slong lastlen = 0;
    fmpz_mod_poly_struct * Fvalue;
    fmpz_t u, v, FvalueA, FvalueB;
    slong Fi, Ti, Aexp, Bexp, e, fexp;
    const fmpz * Acoeff = A->coeffs;
    const fmpz * Bcoeff = B->coeffs;
    slong Flen = F->length;
    fmpz_mod_poly_struct * Fcoeffs = F->coeffs;
    ulong * Fexps = F->exps;
    fmpz_mod_poly_struct * Tcoeffs;
    ulong * Texps;
    fmpz_mod_poly_t zero;

    zero->alloc = 0;
    zero->length = 0;
    zero->coeffs = NULL;

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(FvalueA);
    fmpz_init(FvalueB);

    Fi = 0;
    Aexp = _fmpz_mod_poly_degree(A);
    Bexp = _fmpz_mod_poly_degree(B);

    fmpz_mod_polyun_fit_length(T, Flen + FLINT_MAX(Aexp, Bexp) + 1, ctx);
    Tcoeffs = T->coeffs;
    Texps = T->exps;
    Ti = 0;

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_evaluate_fmpz(u, modulus, alphapow->coeffs + 1, ctx);
    fmpz_mod_mul(u, u, alphapow->coeffs + 1, ctx);
    fmpz_mod_add(u, u, u, ctx);
    FLINT_ASSERT(fmpz_is_one(u));

    fmpz_mod_neg(v, alphapow->coeffs + 1, ctx);
    fmpz_mod_poly_evaluate_fmpz(u, modulus, v, ctx);
    fmpz_mod_mul(u, u, alphapow->coeffs + 1, ctx);
    fmpz_mod_add(u, u, u, ctx);
    FLINT_ASSERT(fmpz_is_one(u));
#endif

    while (Fi < Flen || Aexp >= 0 || Bexp >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        fexp = e = -WORD(1);
        if (Fi < Flen)
        {
            fexp = e = Fexps[Fi];
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Fcoeffs + Fi, ctx));
            FLINT_ASSERT(_fmpz_mod_poly_degree(Fcoeffs + Fi) < _fmpz_mod_poly_degree(modulus));
        }
        if (Aexp >= 0)
        {
            e = FLINT_MAX(e, Aexp);
            FLINT_ASSERT(!fmpz_is_zero(Acoeff + Aexp));
        }
        if (Bexp >= 0)
        {
            e = FLINT_MAX(e, Bexp);
            FLINT_ASSERT(!fmpz_is_zero(Bcoeff + Bexp));
        }

        FLINT_ASSERT(e >= 0);

        Texps[Ti] = e;

        Fvalue = zero;
        fmpz_zero(FvalueA);
        fmpz_zero(FvalueB);
        Finc = 0;
        if (Fi < Flen && e == fexp)
        {
            Finc = 1;
            Fvalue = Fcoeffs + Fi;
            fmpz_mod_poly_eval2_pow(FvalueA, FvalueB, Fcoeffs + Fi, alphapow, ctx);
        }

        if (e == Aexp)
        {
            fmpz_mod_sub(FvalueA, FvalueA, Acoeff + Aexp, ctx);
        }
        if (e == Bexp)
        {
            fmpz_mod_sub(FvalueB, FvalueB, Bcoeff + Bexp, ctx);
        }

        fmpz_mod_sub(u, FvalueB, FvalueA, ctx);
        fmpz_mod_add(v, FvalueB, FvalueA, ctx);
        fmpz_mod_mul(v, v, alphapow->coeffs + 1, ctx);
        fmpz_mod_neg(v, v, ctx);
        changed |= !fmpz_is_zero(u) || !fmpz_is_zero(v);
        fmpz_mod_poly_addmul_linear(Tcoeffs + Ti, Fvalue, modulus, u, v, ctx);

        FLINT_ASSERT(Tcoeffs[Ti].length > 0);
        lastlen = FLINT_MAX(lastlen, Tcoeffs[Ti].length);
        Ti++;

        Fi += Finc;
        if (e == Aexp)
        {
            do {
                Aexp--;
            } while (Aexp >= 0 && fmpz_is_zero(Acoeff + Aexp));
        }
        if (e == Bexp)
        {
            do {
                Bexp--;
            } while (Bexp >= 0 && fmpz_is_zero(Bcoeff + Bexp));
        }
    }
    T->length = Ti;

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(FvalueA);
    fmpz_clear(FvalueB);

    if (changed)
        fmpz_mod_polyun_swap(T, F);

    *lastdeg = lastlen - 1;
    return changed;
}





/*
    E = A(v = alpha)
    A is in R[X][v]
    E is in R[X]
*/
void fmpz_mod_mpolyn_interp_reduce_sm_poly(
    fmpz_mod_poly_t E,
    const fmpz_mod_mpolyn_t A,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t v;
    slong Ai, Alen, k;
    fmpz_mod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N, off, shift;

    fmpz_init(v);

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    fmpz_mod_poly_zero(E, ctx->ffinfo);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fmpz_mod_poly_evaluate_fmpz(v, Acoeff + Ai, alpha, ctx->ffinfo);
        k = (Aexp + N*Ai)[off] >> shift;
        fmpz_mod_poly_set_coeff_fmpz(E, k, v, ctx->ffinfo);
    }

    fmpz_clear(v);
}

/*
    A = B
    A, B are in R[X]
*/
void fmpz_mod_mpolyn_interp_lift_sm_poly(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong Bi;
    slong Blen = fmpz_mod_poly_length(B, ctx->ffinfo);
    fmpz * Bcoeff = B->coeffs;
    fmpz_mod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fmpz_mod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bi = Blen - 1; Bi >= 0; Bi--)
    {
        if (!fmpz_is_zero(Bcoeff + Bi))
        {
            FLINT_ASSERT(Ai < A->alloc);

            fmpz_mod_poly_set_fmpz(Acoeff + Ai, Bcoeff + Bi, ctx->ffinfo);
            mpoly_monomial_zero(Aexp + N*Ai, N);
            (Aexp + N*Ai)[off] = Bi << shift;
            Ai++;
        }
    }
    A->length = Ai;
}


/*
    F = F + modulus*(A - F(v = alpha))
    no assumptions about matching monomials
    F is in Fp[X][v]
    A is in Fp[X]
    it is expected that modulus(alpha) == 1
*/
int fmpz_mod_mpolyn_interp_crt_sm_poly(
    slong * lastdeg_,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t modulus,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    fmpz_t u, v;
    slong Fi, Ti, Ai;
    fmpz * Acoeff = A->coeffs;
    slong Flen = F->length;
    fmpz_mod_poly_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    fmpz_mod_poly_struct * Tcoeff;
    ulong * Texp;
    fmpz_mod_poly_t tp;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    Fi = 0;
    Ai = fmpz_mod_poly_degree(A, ctx->ffinfo);

    fmpz_init(u);
    fmpz_init(v);
    fmpz_mod_poly_init(tp, ctx->ffinfo);

    fmpz_mod_mpolyn_fit_length(T, Flen + Ai + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    while (Fi < Flen || Ai >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Fcoeff + Fi, ctx->ffinfo));
            FLINT_ASSERT(fmpz_mod_poly_degree(Fcoeff + Fi, ctx->ffinfo) <
                                   fmpz_mod_poly_degree(modulus, ctx->ffinfo));
        }

        if (Ai >= 0)
        {
            FLINT_ASSERT(!fmpz_is_zero(Acoeff + Ai));
        }

        mpoly_monomial_zero(Texp + N*Ti, N);

        if (Fi < Flen && Ai >= 0 && ((Fexp + N*Fi)[off]>>shift) == Ai)
        {
            /* F term ok, A term ok */
            fmpz_mod_poly_evaluate_fmpz(u, Fcoeff + Fi, alpha, ctx->ffinfo);
            fmpz_mod_sub(v, Acoeff + Ai, u, ctx->ffinfo);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v, ctx->ffinfo);
                fmpz_mod_poly_add(Tcoeff + Ti, Fcoeff + Fi, tp, ctx->ffinfo);
            }
            else
            {
                fmpz_mod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->ffinfo);
            }

            (Texp + N*Ti)[off] = Ai << shift;
            Fi++;
            do {
                Ai--;
            } while (Ai >= 0 && fmpz_is_zero(Acoeff + Ai));
        }
        else if (Fi < Flen && (Ai < 0 || ((Fexp + N*Fi)[off]>>shift) > Ai))
        {
            /* F term ok, A term missing */
            fmpz_mod_poly_evaluate_fmpz(v, Fcoeff + Fi, alpha, ctx->ffinfo);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v, ctx->ffinfo);
                fmpz_mod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, tp, ctx->ffinfo);
            }
            else
            {
                fmpz_mod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->ffinfo);                
            }

            (Texp + N*Ti)[off] = (Fexp + N*Fi)[off];
            Fi++;
        }
        else if (Ai >= 0 && (Fi >= Flen || ((Fexp + N*Fi)[off]>>shift) < Ai))
        {
            /* F term missing, A term ok */
            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz(Tcoeff + Ti, modulus, Acoeff + Ai, ctx->ffinfo);

            (Texp + N*Ti)[off] = Ai << shift;
            do {
                Ai--;
            } while (Ai >= 0 && fmpz_is_zero(Acoeff + Ai));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + Ti, ctx->ffinfo));
        lastdeg = FLINT_MAX(lastdeg, fmpz_mod_poly_degree(Tcoeff + Ti, ctx->ffinfo));

        Ti++;
    }
    T->length = Ti;

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_mod_poly_clear(tp, ctx->ffinfo);

    if (changed)
    {
        fmpz_mod_mpolyn_swap(T, F, ctx);
    }

    *lastdeg_ = lastdeg;
    return changed;
}



void fmpz_mod_mpolyn_interp_reduce_2sm_mpolyn(
    fmpz_mod_mpolyn_t E,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t A,
    slong var,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    fmpz_t e, f;
    fmpz_mod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fmpz_mod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;
    fmpz_mod_poly_struct * Fcoeff;
    ulong * Fexp;
    slong Fi;

    fmpz_init(e);
    fmpz_init(f);

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
        fmpz_mod_poly_eval2_pow(e, f, Acoeff + Ai, alphapow, ctx->ffinfo);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;

        if (!fmpz_is_zero(e))
        {
            if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
            {
                /* append to previous */
                fmpz_mod_poly_set_coeff_fmpz(Ecoeff + Ei - 1, k, e, ctx->ffinfo);
            }
            else
            {
                FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

                /* create new */
                if (Ei >= E->alloc)
                {
                    fmpz_mod_mpolyn_fit_length(E, Ei + 1, ctx);
                    Ecoeff = E->coeffs;
                    Eexp = E->exps;
                }
                mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
                fmpz_mod_poly_zero(Ecoeff + Ei, ctx->ffinfo);
                fmpz_mod_poly_set_coeff_fmpz(Ecoeff + Ei, k, e, ctx->ffinfo);
                Ei++;
            }
        }

        if (!fmpz_is_zero(f))
        {
            if (Fi > 0 && mpoly_monomial_equal_extra(Fexp + N*(Fi - 1), Aexp + N*Ai, N, offset, -(k << shift)))
            {
                /* append to previous */
                fmpz_mod_poly_set_coeff_fmpz(Fcoeff + Fi - 1, k, f, ctx->ffinfo);
            }
            else
            {
                FLINT_ASSERT(Fi == 0 || mpoly_monomial_gt_nomask_extra(Fexp + N*(Fi - 1), Aexp + N*Ai, N, offset, -(k << shift)));

                /* create new */
                if (Fi >= F->alloc)
                {
                    fmpz_mod_mpolyn_fit_length(F, Fi + 1, ctx);
                    Fcoeff = F->coeffs;
                    Fexp = F->exps;
                }
                mpoly_monomial_set_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, -(k << shift));
                fmpz_mod_poly_zero(Fcoeff + Fi, ctx->ffinfo);
                fmpz_mod_poly_set_coeff_fmpz(Fcoeff + Fi, k, f, ctx->ffinfo);
                Fi++;
            }
        }
    }
    E->length = Ei;
    F->length = Fi;

    fmpz_clear(e);
    fmpz_clear(f);
}


void fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(
    slong * lastdeg,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong lastlen = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;

    fmpz_mod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    fmpz_mod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai, ai;

    fmpz_mod_poly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    slong Bi, bi;

    fmpz zero0 = 0;
    fmpz * Avalue, * Bvalue;
    fmpz_t u, v, FvalueA, FvalueB;
    int cmp;
    fmpz_t d0;

    fmpz_init(d0);
    fmpz_mod_add(d0, alpha, alpha, ctx->ffinfo);
    fmpz_mod_inv(d0, d0, ctx->ffinfo);

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(FvalueA);
    fmpz_init(FvalueB);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    fmpz_mod_mpolyn_fit_length(T, FLINT_MAX(Alen, Blen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;

    Ti = Ai = Bi = 0;
    ai = (Ai >= Alen) ? 0 : A->coeffs[Ai].length - 1;
    bi = (Bi >= Blen) ? 0 : B->coeffs[Bi].length - 1;

    while (Ai < Alen || Bi < Blen)
    {
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Alen - Ai, Blen - Bi);
            fmpz_mod_mpolyn_fit_length(T, Ti + extra + 1, ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        FLINT_ASSERT(Ai >= Alen || !fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
        FLINT_ASSERT(Bi >= Blen || !fmpz_is_zero((Bcoeff + Bi)->coeffs + bi));

        Avalue = &zero0;
        if (Ai < Alen)
        {
            Avalue = (Acoeff + Ai)->coeffs + ai;
            mpoly_monomial_set_extra(Texp + N*Ti,
                                     Aexp + N*Ai, N, offset, ai << shift);
        }

        Bvalue = &zero0;
        if (Bi < Blen)
        {
            cmp = (Avalue == &zero0) ? -1 : mpoly_monomial_cmp_nomask_extra(
                             Texp + N*Ti, Bexp + N*Bi, N, offset, bi << shift);
            if (cmp <= 0)
            {
                Bvalue = (Bcoeff + Bi)->coeffs + bi;
            }
            if (cmp < 0)
            {
                Avalue = &zero0;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Bexp + N*Bi, N, offset, bi << shift);
            }
        }

        fmpz_mod_neg(FvalueA, Avalue, ctx->ffinfo);
        fmpz_mod_neg(FvalueB, Bvalue, ctx->ffinfo);
        fmpz_mod_sub(u, FvalueB, FvalueA, ctx->ffinfo);
        fmpz_mod_add(v, FvalueB, FvalueA, ctx->ffinfo);
        fmpz_mod_mul(v, alpha, v, ctx->ffinfo);
        fmpz_mod_neg(v, v, ctx->ffinfo);

        FLINT_ASSERT(!fmpz_is_zero(u) || !fmpz_is_zero(v));
        fmpz_mod_poly_zero(Tcoeff + Ti, ctx->ffinfo);
        fmpz_mod_mul(u, u, d0, ctx->ffinfo);
        fmpz_mod_mul(v, v, d0, ctx->ffinfo);
        fmpz_mod_poly_set_coeff_fmpz(Tcoeff + Ti, 0, v, ctx->ffinfo);
        fmpz_mod_poly_set_coeff_fmpz(Tcoeff + Ti, 1, u, ctx->ffinfo);

        if (Avalue != &zero0)
        {
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                    ai = A->coeffs[Ai].length - 1;
            }
        }
        if (Bvalue != &zero0)
        {
            do {
                bi--;
            } while (bi >= 0 && fmpz_is_zero((Bcoeff + Bi)->coeffs + bi));
            if (bi < 0)
            {
                Bi++;
                if (Bi < Blen)
                    bi = B->coeffs[Bi].length - 1;
            }
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + Ti, ctx->ffinfo));
        lastlen = FLINT_MAX(lastlen, Tcoeff[Ti].length);
        Ti++;
    }
    T->length = Ti;

    *lastdeg = lastlen - 1;

    FLINT_ASSERT(fmpz_mod_mpolyn_is_canonical(T, ctx));

    fmpz_clear(d0);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(FvalueA);
    fmpz_clear(FvalueB);

    return;
}


int fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    slong var,
    fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong lastlen = 0;
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;
    fmpz_mod_poly_t zero;
    fmpz zero0 = 0;

    fmpz_mod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    fmpz_mod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    fmpz_mod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai, ai;

    fmpz_mod_poly_struct * Bcoeff = B->coeffs;
    slong Blen = B->length;
    ulong * Bexp = B->exps;
    slong Bi, bi;

    fmpz_mod_poly_struct * Fvalue;
    fmpz * Avalue, * Bvalue;
    fmpz_t u, v, FvalueA, FvalueB;
    int texp_set, cmp;

    FLINT_ASSERT(fmpz_mod_poly_degree(modulus, ctx->ffinfo) > 0);

    zero->coeffs = NULL;
    zero->alloc = 0;
    zero->length = 0;

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(FvalueA);
    fmpz_init(FvalueB);

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_evaluate_fmpz(u, modulus, alphapow->coeffs + 1, ctx->ffinfo);
    fmpz_mod_mul(u, u, alphapow->coeffs + 1, ctx->ffinfo);
    fmpz_mod_add(u, u, u, ctx->ffinfo);
    FLINT_ASSERT(fmpz_is_one(u));
    fmpz_mod_neg(v, alphapow->coeffs + 1, ctx->ffinfo);
    fmpz_mod_poly_evaluate_fmpz(u, modulus, v, ctx->ffinfo);
    fmpz_mod_mul(u, u, alphapow->coeffs + 1, ctx->ffinfo);
    fmpz_mod_add(u, u, u, ctx->ffinfo);
    FLINT_ASSERT(fmpz_is_one(u));
#endif

    FLINT_ASSERT(fmpz_mod_mpolyn_is_canonical(A, ctx));
    FLINT_ASSERT(fmpz_mod_mpolyn_is_canonical(B, ctx));
    FLINT_ASSERT(fmpz_mod_mpolyn_is_canonical(F, ctx));

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    fmpz_mod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;

    Ti = Fi = Ai = Bi = 0;
    ai = (Ai >= Alen) ? 0 : A->coeffs[Ai].length - 1;
    bi = (Bi >= Blen) ? 0 : B->coeffs[Bi].length - 1;

    while (Fi < Flen || Ai < Alen || Bi < Blen)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Alen - Ai);
            extra = FLINT_MAX(extra, Blen - Bi);
            fmpz_mod_mpolyn_fit_length(T, Ti + extra + 1, ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        FLINT_ASSERT(Fi >= Flen || (Fcoeff + Fi)->length != 0);
        FLINT_ASSERT(Ai >= Alen || !fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
        FLINT_ASSERT(Bi >= Blen || !fmpz_is_zero((Bcoeff + Bi)->coeffs + bi));

        Fvalue = zero;
        texp_set = 0;
        if (Fi < Flen)
        {
            Fvalue = Fcoeff + Fi;
            texp_set = 1;
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);
        }

        Avalue = &zero0;
        if (Ai < Alen)
        {
            cmp = (!texp_set) ? -1
                     : mpoly_monomial_cmp_nomask_extra(Texp + N*Ti,
                                      Aexp + N*Ai, N, offset, ai << shift);

            if (cmp <= 0)
            {
                Avalue = (Acoeff + Ai)->coeffs + ai;
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                texp_set = 1;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Aexp + N*Ai, N, offset, ai << shift);
            }
        }

        Bvalue = &zero0;
        if (Bi < Blen)
        {
            cmp = (!texp_set) ? -1 : mpoly_monomial_cmp_nomask_extra(
                             Texp + N*Ti, Bexp + N*Bi, N, offset, bi << shift);

            if (cmp <= 0)
            {
                Bvalue = (Bcoeff + Bi)->coeffs + bi;
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                Avalue = &zero0;
                texp_set = 1;
                mpoly_monomial_set_extra(Texp + N*Ti,
                                         Bexp + N*Bi, N, offset, bi << shift);
            }
        }

        FLINT_ASSERT(texp_set);

        fmpz_mod_poly_eval2_pow(FvalueA, FvalueB, Fvalue, alphapow, ctx->ffinfo);
        fmpz_mod_sub(FvalueA, FvalueA, Avalue, ctx->ffinfo);
        fmpz_mod_sub(FvalueB, FvalueB, Bvalue, ctx->ffinfo);
        fmpz_mod_sub(u, FvalueB, FvalueA, ctx->ffinfo);
        fmpz_mod_add(v, FvalueB, FvalueA, ctx->ffinfo);
        fmpz_mod_mul(v, alphapow->coeffs + 1, v, ctx->ffinfo);
        fmpz_mod_neg(v, v, ctx->ffinfo);
        changed |= !fmpz_is_zero(u) || !fmpz_is_zero(v);
        fmpz_mod_poly_addmul_linear(Tcoeff + Ti, Fvalue, modulus, u, v, ctx->ffinfo);

        Fi += (Fvalue != zero);
        if (Avalue != &zero0)
        {
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                    ai = A->coeffs[Ai].length - 1;
            }
        }
        if (Bvalue != &zero0)
        {
            do {
                bi--;
            } while (bi >= 0 && fmpz_is_zero((Bcoeff + Bi)->coeffs + bi));
            if (bi < 0)
            {
                Bi++;
                if (Bi < Blen)
                    bi = B->coeffs[Bi].length - 1;
            }
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + Ti, ctx->ffinfo));
        lastlen = FLINT_MAX(lastlen, Tcoeff[Ti].length);
        Ti++;
    }
    T->length = Ti;

    *lastdeg = lastlen - 1;

    if (changed)
        fmpz_mod_mpolyn_swap(T, F, ctx);

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(FvalueA);
    fmpz_clear(FvalueB);

    FLINT_ASSERT(fmpz_mod_mpolyn_is_canonical(F, ctx));

    return changed;
}



void fmpz_mod_mpolyn_interp_lift_sm_mpoly(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N;
    fmpz_mod_poly_struct * Acoeff;
    const fmpz * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Blen;

    FLINT_ASSERT(A->bits == B->bits);

    Blen = B->length;
    fmpz_mod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < Blen; i++)
    {
        fmpz_mod_poly_set_fmpz(Acoeff + i, Bcoeff + i, ctx->ffinfo);
        mpoly_monomial_set(Aexp + N*i, Bexp + N*i, N);
    }
    A->length = Blen;
}

int fmpz_mod_mpolyn_interp_crt_sm_mpoly(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_mpoly_t A,
    fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong lastlen = 0;
    int changed = 0;
    slong i, j, k;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    fmpz * Acoeff = A->coeffs;
    fmpz_mod_poly_struct * Fcoeff = F->coeffs;
    fmpz_mod_poly_struct * Tcoeff;
    fmpz_t v;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fmpz_mod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;


    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen
                       || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Fcoeff + i, ctx->ffinfo));
            FLINT_ASSERT(fmpz_mod_poly_degree(Fcoeff + i, ctx->ffinfo) <
                                   fmpz_mod_poly_degree(modulus, ctx->ffinfo));

            /* F term ok, A term missing */
            fmpz_mod_poly_eval_pow(v, Fcoeff + i, alphapow, ctx->ffinfo);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_neg(v, v, ctx->ffinfo);
                fmpz_mod_poly_scalar_addmul_fmpz_mod(Tcoeff + k,
                                          Fcoeff + i, modulus, v, ctx->ffinfo);
            }
            else
            {
                fmpz_mod_poly_set(Tcoeff + k, Fcoeff + i, ctx->ffinfo);                
            }

            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + k, ctx->ffinfo));
            lastlen = FLINT_MAX(lastlen, Tcoeff[k].length);

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            FLINT_ASSERT(!fmpz_is_zero(Acoeff + j));

            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz(Tcoeff + k, modulus, Acoeff + j, ctx->ffinfo);
            lastlen = FLINT_MAX(lastlen, Tcoeff[k].length);
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            k++;
            j++;
        }
        else if (i < Flen && j < Alen
                            && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Fcoeff + i, ctx->ffinfo));
            FLINT_ASSERT(fmpz_mod_poly_degree(Fcoeff + i, ctx->ffinfo) <
                                   fmpz_mod_poly_degree(modulus, ctx->ffinfo));

            /* F term ok, A term ok */
            fmpz_mod_poly_eval_pow(v, Fcoeff + i, alphapow, ctx->ffinfo);
            fmpz_mod_sub(v, Acoeff + j, v, ctx->ffinfo);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_addmul_fmpz_mod(Tcoeff + k,
                                          Fcoeff + i, modulus, v, ctx->ffinfo);
            }
            else
            {
                fmpz_mod_poly_set(Tcoeff + k, Fcoeff + i, ctx->ffinfo);            
            }
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + k, ctx->ffinfo));
            lastlen = FLINT_MAX(lastlen, Tcoeff[k].length);
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
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

    *lastdeg = lastlen - 1;

    if (changed)
        fmpz_mod_mpolyn_swap(T, F, ctx);

    return changed;
}

int fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong lastlen = 0;
    int changed = 0;
    slong i;
    fmpz_t v;
    fmpz * Acoeff = A->coeffs;
    slong Flen = F->length;

    FLINT_ASSERT(Flen == A->length);

    fmpz_init(v);

    for (i = 0; i < Flen; i++)
    {
        /* F term ok, A term ok */
        fmpz_mod_poly_eval_pow(v, F->coeffs + i, alphapow, ctx->ffinfo);
        fmpz_mod_sub(v, Acoeff + i, v, ctx->ffinfo);
        if (!fmpz_is_zero(v))
        {
            changed = 1;
            fmpz_mod_poly_scalar_addmul_fmpz_mod(F->coeffs + i,
                                       F->coeffs + i, modulus, v, ctx->ffinfo);
        }
        lastlen = FLINT_MAX(lastlen, F->coeffs[i].length);
    }

    fmpz_clear(v);

    *lastdeg = lastlen - 1;
    return changed;
}

