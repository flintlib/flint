/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/*
    interp_reduce: map from Fp[x] to Fp[x]/poly(x)
    interp_lift:   map from Fp[x]/poly(x) to Fp[x]
    interp_crt:    update element of Fp[x] with a new image in Fp[x]/poly(x)
    interp_mcrt:   same as interp_mcrt, but monomial match, thus easier

    ..._sm:    poly(x) is x - alpha
    ..._lg:    poly(x) is modulus of finite field
*/

/****************** bivar - image in F_q *************************************/

/*
    E = A mod alpha(v)
    A is in Fp[X][][v]
    E is in (Fp/alpha(v))[X]
*/

void nmod_mpolyn_interp_reduce_lg_poly(
    fq_nmod_poly_t E,
    const fq_nmod_ctx_t fqctx,
    const nmod_mpolyn_t A,
    const nmod_mpoly_ctx_t ctx)
{
    fq_nmod_t v;
    slong Ai, Alen, k;
    nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fq_nmod_init(v, fqctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;

    fq_nmod_poly_zero(E, fqctx);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        k = (Aexp + N*Ai)[off] >> shift;
        nmod_poly_rem(v, Acoeff + Ai, fqctx->modulus);
        fq_nmod_poly_set_coeff(E, k, v, fqctx);
    }

    fq_nmod_clear(v, fqctx);
}


/*
    Convert B to A using the lowest degree representative
    A is in           Fp [X][][v]
    B is in (Fp/alpha(v))[X]
*/

void nmod_mpolyn_interp_lift_lg_poly(
    slong * lastdeg_,
    nmod_mpolyn_t A,
    const nmod_mpoly_ctx_t ctx,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t fqctx)
{
    slong Bexp;
    slong Blen = fq_nmod_poly_length(B, fqctx);
    fq_nmod_struct * Bcoeff = B->coeffs;
    nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong lastdeg = -WORD(1);
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bexp = Blen - 1; Bexp >= 0; Bexp--)
    {
        if (!fq_nmod_is_zero(Bcoeff + Bexp, fqctx))
        {
            FLINT_ASSERT(Ai < A->alloc);
            nmod_poly_set(Acoeff + Ai, Bcoeff + Bexp);
            lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Bcoeff + Bexp));
            mpoly_monomial_zero(Aexp + N*Ai, N);
            (Aexp + N*Ai)[off] = Bexp << shift;
            Ai++;
        }
    }
    A->length = Ai;

    *lastdeg_ = lastdeg;
}


/*
    F = F + modulus*((A - F(alpha))/(modulus(alpha)))
    no assumptions about matching monomials
    F is in Fp[X][v]
    A is in (Fp/alpha(v))[X]    alpha(v) is fqctx->modulus(v)
*/

int nmod_mpolyn_interp_crt_lg_poly(
    slong * lastdeg_,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    nmod_poly_t modulus,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_poly_t A,
    const fq_nmod_ctx_t fqctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    fq_nmod_t u, v;
    nmod_poly_t w;
    slong Fi, Ti, Aexp;
    fq_nmod_struct * Acoeff = A->coeffs;
    slong Flen = F->length;
    nmod_poly_struct * Fcoeffs = F->coeffs;
    ulong * Fexps = F->exps;
    nmod_poly_struct * Tcoeffs;
    ulong * Texps;
    nmod_poly_t tp;
    fq_nmod_t inv_m_eval;
    slong N, off, shift;

    FLINT_ASSERT(T->bits == F->bits);

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    fq_nmod_init(inv_m_eval, fqctx);
    nmod_poly_rem(inv_m_eval, modulus, fqctx->modulus);
    fq_nmod_inv(inv_m_eval, inv_m_eval, fqctx);

    Fi = 0;

    fq_nmod_init(u, fqctx);
    fq_nmod_init(v, fqctx);
    nmod_poly_init(w, fqctx->modulus->mod.n);

    Aexp = fq_nmod_poly_degree(A, fqctx);

    nmod_poly_init(tp, ctx->ffinfo->mod.n);

    nmod_mpolyn_fit_length(T, Flen + Aexp + 1, ctx);
    Tcoeffs = T->coeffs;
    Texps = T->exps;
    Ti = 0;

    while (Fi < Flen || Aexp >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeffs + Fi));
            FLINT_ASSERT(nmod_poly_degree(Fcoeffs + Fi) < nmod_poly_degree(modulus));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(!fq_nmod_is_zero(Acoeff + Aexp, fqctx));
        }

        mpoly_monomial_zero(Texps + N*Ti, N);

        if (Fi < Flen && Aexp >= 0 && ((Fexps + N*Fi)[off]>>shift) == Aexp)
        {
            /* F term ok, A term ok */
            nmod_poly_rem(u, Fcoeffs + Fi, fqctx->modulus);
            fq_nmod_sub(v, Acoeff + Aexp, u, fqctx);
            if (!fq_nmod_is_zero(v, fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, fqctx);
                nmod_poly_mul(w, modulus, u);
                nmod_poly_add(Tcoeffs + Ti, Fcoeffs + Fi, w);
            }
            else
            {
                nmod_poly_set(Tcoeffs + Ti, Fcoeffs + Fi);                
            }

            (Texps + N*Ti)[off] = (Fexps + N*Fi)[off];
            Fi++;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, fqctx));
        }
        else if (Fi < Flen && (Aexp < 0 || ((Fexps + N*Fi)[off]>>shift) > Aexp))
        {
            /* F term ok, A term missing */
            nmod_poly_rem(v, Fcoeffs + Fi, fqctx->modulus);
            if (!fq_nmod_is_zero(v, fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, fqctx);
                nmod_poly_mul(w, u, modulus);
                nmod_poly_sub(Tcoeffs + Ti, Fcoeffs + Fi, w);
            }
            else
            {
                nmod_poly_set(Tcoeffs + Ti, Fcoeffs + Fi);                
            }

            (Texps + N*Ti)[off] = (Fexps + N*Fi)[off];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || ((Fexps + N*Fi)[off]>>shift) < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_mul(u, Acoeff + Aexp, inv_m_eval, fqctx);
            nmod_poly_mul(Tcoeffs + Ti, modulus, u);
            (Texps + N*Ti)[off] = Aexp << shift;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, fqctx));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        FLINT_ASSERT(!nmod_poly_is_zero(Tcoeffs + Ti));
        lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Tcoeffs + Ti));
        Ti++;
    }
    T->length = Ti;

    if (changed)
    {
        nmod_mpolyn_swap(T, F);
    }

    fq_nmod_clear(u, fqctx);
    fq_nmod_clear(v, fqctx);
    nmod_poly_clear(w);

    fq_nmod_clear(inv_m_eval, fqctx);

    *lastdeg_ = lastdeg;
    return changed;
}


/*********************** multivar - image in F_q *****************************/

/*
    E = A mod alpha(x_var)
    A is in                Fp[x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in (Fp/alpha(x_var))[x_0, ..., x_(var-2)][x_(var-1)]
    alpha is ectx->modulus
*/

void nmod_mpolyn_interp_reduce_lg_mpolyn(
    fq_nmod_mpolyn_t E,
    fq_nmod_mpoly_ctx_t ectx,
    nmod_mpolyn_t A,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    fq_nmod_t v;
    nmod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fq_nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    fq_nmod_init(v, ectx->fqctx);

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        nmod_poly_rem(v, Acoeff + Ai, ectx->fqctx->modulus);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (fq_nmod_is_zero(v, ectx->fqctx))
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            fq_nmod_poly_set_coeff(Ecoeff + Ei - 1, k, v, ectx->fqctx);
        }
        else
        {
            FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

            /* create new */
            if (Ei >= E->alloc)
            {
                fq_nmod_mpolyn_fit_length(E, Ei + 1, ectx);
                Ecoeff = E->coeffs;
                Eexp = E->exps;
            }
            mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
            fq_nmod_poly_zero(Ecoeff + Ei, ectx->fqctx);
            fq_nmod_poly_set_coeff(Ecoeff + Ei, k, v, ectx->fqctx);
            Ei++;
        }
    }
    E->length = Ei;

    fq_nmod_clear(v, ectx->fqctx);
}


/*
    A = B using lowest degree representative
    A is in                      Fp [x_0, ..., x_(var-1), x_(var-1)][x_var]
    B is in (Fp[x_var]/alpha(x_var))[x_0, ..., x_(var-2)][x_(var-1)]
*/
void nmod_mpolyn_interp_lift_lg_mpolyn(
    slong * lastdeg_,
    nmod_mpolyn_t A,
    slong var,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_mpolyn_t B,
    const fq_nmod_mpoly_ctx_t ectx)
{
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;
    fq_nmod_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;
    nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong lastdeg = -WORD(1);

    FLINT_ASSERT(A->bits == B->bits);

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
            if (!fq_nmod_is_zero((Bcoeff + Bi)->coeffs + vi, ectx->fqctx))
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                nmod_poly_set(Acoeff + Ai, (Bcoeff + Bi)->coeffs + vi);
                lastdeg = FLINT_MAX(lastdeg, nmod_poly_degree(Acoeff + Ai));
                Ai++;
            }
        }
    }
    A->length = Ai;

    *lastdeg_ = lastdeg;
}


int nmod_mpolyn_interp_crt_lg_mpolyn(
    slong * lastdeg_,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    nmod_poly_t modulus,
    slong var,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_mpolyn_t A,
    const fq_nmod_mpoly_ctx_t ectx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;
    slong lastdeg = -WORD(1);
    slong vi;
    fq_nmod_t u, v;
    nmod_poly_t w;

    nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    fq_nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    nmod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;
    fq_nmod_t inv_m_eval;

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    nmod_poly_rem(inv_m_eval, modulus, ectx->fqctx->modulus);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    nmod_poly_init(w, ctx->ffinfo->mod.n);

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
        vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
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
            nmod_poly_rem(u, Fcoeff + Fi, ectx->fqctx->modulus);
            fq_nmod_sub(v, (Acoeff + Ai)->coeffs + vi, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                nmod_poly_mul(w, modulus, u);
                nmod_poly_add(Tcoeff + Ti, Fcoeff + Fi, w);
            }
            else
            {
                nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ectx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            nmod_poly_rem(v, Fcoeff + Fi, ectx->fqctx->modulus);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                nmod_poly_mul(w, u, modulus);
                nmod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, w);
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
            fq_nmod_mul(u, (Acoeff + Ai)->coeffs + vi, inv_m_eval, ectx->fqctx);
            nmod_poly_mul(Tcoeff + Ti, modulus, u);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ectx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
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

    fq_nmod_clear(inv_m_eval, ectx->fqctx);
    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    nmod_poly_clear(w);

    FLINT_ASSERT(nmod_mpolyn_is_canonical(F, ctx));

    *lastdeg_ = lastdeg;

    return changed;
}


/******************** multivar - image in Fq *********************************/

/*
    A = B mod alpha(x_var)
    B is in               Fp [x_0, ..., x_(var-2), x_(var-1)][x_var]
    A is in (Fp/alpha(x_var))[x_0, ..., x_(var-2), x_(var-1)]
    alpha is ectx->modulus
*/
void nmod_mpolyn_interp_reduce_lg_mpoly(
    fq_nmod_mpoly_t A,
    nmod_mpolyn_t B,
    const fq_nmod_mpoly_ctx_t ffctx,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, k, N;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ffctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ffctx->minfo->nvars == ctx->minfo->nvars);

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    k = 0;
    fq_nmod_mpoly_fit_length(A, B->length, ffctx);
    for (i = 0; i < B->length; i++)
    {
        mpoly_monomial_set(A->exps + N*k, B->exps + N*i, N);
        nmod_poly_rem(A->coeffs + k, B->coeffs + i, ffctx->fqctx->modulus);
        k += !fq_nmod_is_zero(A->coeffs + k, ffctx->fqctx);
    }
    A->length = k;
}

/*
    A = B mod alpha(x_var)
    B is in               Fp [X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    A is in (Fp/alpha(x_var))[X][x_0, ..., x_(var-2), x_(var-1)]
    alpha is ectx->modulus
*/
void nmod_mpolyun_interp_reduce_lg_mpolyu(
    fq_nmod_mpolyu_t A,
    nmod_mpolyun_t B,
    const fq_nmod_mpoly_ctx_t ffctx,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, k, Blen;
    fq_nmod_mpoly_struct * Acoeff;
    nmod_mpolyn_struct * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    fq_nmod_mpolyu_fit_length(A, Blen, ffctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    k = 0;
    for (i = 0; i < Blen; i++)
    {
        nmod_mpolyn_interp_reduce_lg_mpoly(Acoeff + k, Bcoeff + i, ffctx, ctx);
        Aexp[k] = Bexp[i];
        k += !fq_nmod_mpoly_is_zero(Acoeff + k, ffctx);
    }
    A->length = k;  
}


/*
    A = Ap using lowest degree
    A  is in               Fp [x_0, ..., x_(var-2), x_(var-1)][x_var]
    Ap is in (Fp/alpha(x_var))[x_0, ..., x_(var-2), x_(var-1)]
    alpha is ectx->modulus
*/
void nmod_mpolyn_interp_lift_lg_mpoly(
    nmod_mpolyn_t A,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_mpoly_t Ap,
    const fq_nmod_mpoly_ctx_t ctxp)
{
    slong i, N;

    FLINT_ASSERT(Ap->bits == A->bits);
    nmod_mpolyn_fit_length(A, Ap->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    for (i = 0; i < Ap->length; i++)
    {
        mpoly_monomial_set(A->exps + N*i, Ap->exps + N*i, N);
        nmod_poly_set(A->coeffs + i, Ap->coeffs + i);
    }
    A->length = Ap->length;
}

/*
    A = Ap using lowest degree
    A  is in               Fp [X][x_0, ..., x_(var-2), x_(var-1)][x_var]
    Ap is in (Fp/alpha(x_var))[X][x_0, ..., x_(var-2), x_(var-1)]
    alpha is ectx->modulus
*/
void nmod_mpolyun_interp_lift_lg_mpolyu(
    nmod_mpolyun_t A,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_mpolyu_t Ap,
    const fq_nmod_mpoly_ctx_t ctxp)
{
    slong i;

    FLINT_ASSERT(Ap->bits == A->bits);
    nmod_mpolyun_fit_length(A, Ap->length, ctx);
    for (i = 0; i < Ap->length; i++)
    {
        A->exps[i] = Ap->exps[i];
        nmod_mpolyn_interp_lift_lg_mpoly(A->coeffs + i, ctx, Ap->coeffs + i, ctxp);
    }
    A->length = Ap->length;
}

/*
    F = F + modulus*((A - F(alpha))/(modulus(alpha)))
    no assumptions about matching monomials
*/
int nmod_mpolyn_interp_crt_lg_mpoly(
    slong * lastdeg,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    nmod_poly_t m,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_mpoly_t A,
    fq_nmod_t inv_m_eval,
    const fq_nmod_mpoly_ctx_t ffctx)
{
    int changed = 0;
    slong i, j, k;
    slong N;
    fq_nmod_t u, v;
    nmod_poly_t w;
    flint_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    fq_nmod_struct * Acoeff = A->coeffs;
    nmod_poly_struct * Fcoeff = F->coeffs;
    nmod_poly_struct * Tcoeff;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fq_nmod_init(u, ffctx->fqctx);
    fq_nmod_init(v, ffctx->fqctx);
    nmod_poly_init(w, ffctx->fqctx->modulus->mod.n);

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
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + i) < nmod_poly_degree(m));

            /* F term ok, A term missing */
            nmod_poly_rem(v, Fcoeff + i, ffctx->fqctx->modulus);
            if (!fq_nmod_is_zero(v, ffctx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ffctx->fqctx);
                nmod_poly_mul(w, u, m);
                nmod_poly_sub(Tcoeff + k, Fcoeff + i, w);
            } else {
                nmod_poly_set(Tcoeff + k, Fcoeff + i);
            }
            *lastdeg = FLINT_MAX(*lastdeg, nmod_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            FLINT_ASSERT(!nmod_poly_is_zero(Tcoeff + k));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (!fq_nmod_is_zero(Acoeff + j, ffctx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, Acoeff + j, inv_m_eval, ffctx->fqctx);
                nmod_poly_mul(Tcoeff + k, m, u);
                *lastdeg = FLINT_MAX(*lastdeg, nmod_poly_degree(Tcoeff + k));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!nmod_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(nmod_poly_degree(Fcoeff + i) < nmod_poly_degree(m));

            /* F term ok, A term ok */
            nmod_poly_rem(u, Fcoeff + i, ffctx->fqctx->modulus);
            fq_nmod_sub(v, Acoeff + j, u, ffctx->fqctx);
            if (!fq_nmod_is_zero(v, ffctx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ffctx->fqctx);
                nmod_poly_mul(w, m, u);
                nmod_poly_add(Tcoeff + k, Fcoeff + i, w);
            } else {
                nmod_poly_set(Tcoeff + k, Fcoeff + i);                
            }
            *lastdeg = FLINT_MAX(*lastdeg, nmod_poly_degree(Tcoeff + k));
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

    fq_nmod_clear(u, ffctx->fqctx);
    fq_nmod_clear(v, ffctx->fqctx);
    nmod_poly_clear(w);

    return changed;
}

/*
    F = F + modulus*((A - F(alpha))/(modulus(alpha)))
    no assumptions about matching monomials
*/
int nmod_mpolyun_interp_crt_lg_mpolyu(
    slong * lastdeg,
    nmod_mpolyun_t F,
    nmod_mpolyun_t T,
    nmod_poly_t m,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ffctx)
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
    fq_nmod_mpoly_struct  * Acoeff;
    fq_nmod_mpoly_t zero;
    fq_nmod_t inv_m_eval;

    *lastdeg = -WORD(1);

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

    fq_nmod_mpoly_init(zero, ffctx);
    fq_nmod_mpoly_fit_bits(zero, A->bits, ffctx);
    zero->bits = A->bits;

    fq_nmod_init(inv_m_eval, ffctx->fqctx);
    nmod_poly_rem(inv_m_eval, m, ffctx->fqctx->modulus);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ffctx->fqctx);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= nmod_mpolyn_interp_crt_lg_mpoly(lastdeg, Tcoeff + k,
                                           S, m, ctx, zero, inv_m_eval, ffctx);
            Texp[k] = Fexp[i];
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || Aexp[j] > Fexp[i]))
        {
            /* F term missing, A term ok */
            nmod_mpolyn_zero(Tcoeff + k, ctx);
            changed |= nmod_mpolyn_interp_crt_lg_mpoly(lastdeg, Tcoeff + k,
                                     S, m, ctx, Acoeff + j, inv_m_eval, ffctx);
            Texp[k] = Aexp[j];
            k++;
            j++;
        }
        else if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= nmod_mpolyn_interp_crt_lg_mpoly(lastdeg, Tcoeff + k,
                                     S, m, ctx, Acoeff + j, inv_m_eval, ffctx);
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

    fq_nmod_clear(inv_m_eval, ffctx->fqctx);

    nmod_mpolyn_clear(S, ctx);
    fq_nmod_mpoly_clear(zero, ffctx);
    return changed;    
}


/*
    Update H so that it does not change mod m, and is now A mod p
    It is asserted that the monomials in H and A match
*/
int nmod_mpolyn_CRT_fq_nmod_mpoly(
    slong * lastdeg,
    nmod_mpolyn_t H,
    const nmod_mpoly_ctx_t ctx,
    nmod_poly_t m,
    fq_nmod_t inv_m_eval,
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctxp)
{
    slong i;
#if WANT_ASSERT
    slong N;
#endif
    int changed = 0;
    fq_nmod_t u, v;
    nmod_poly_t w;

    fq_nmod_init(u, ctxp->fqctx);
    fq_nmod_init(v, ctxp->fqctx);
    nmod_poly_init(w, ctxp->fqctx->modulus->mod.n);

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);
#if WANT_ASSERT
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
#endif
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        nmod_poly_rem(u, H->coeffs + i, ctxp->fqctx->modulus);
        fq_nmod_sub(v, A->coeffs + i, u, ctxp->fqctx);
        if (!fq_nmod_is_zero(v, ctxp->fqctx))
        {
            changed = 1;
            fq_nmod_mul(u, v, inv_m_eval, ctxp->fqctx);
            nmod_poly_mul(w, u, m);
            nmod_poly_add(H->coeffs + i, H->coeffs + i, w);
        }

        *lastdeg = FLINT_MAX(*lastdeg, nmod_poly_degree(H->coeffs + i));

        FLINT_ASSERT(nmod_poly_degree(H->coeffs + i) < nmod_poly_degree(m)
                                                     + nmod_poly_degree(ctxp->fqctx->modulus));
    }

    fq_nmod_clear(u, ctxp->fqctx);
    fq_nmod_clear(v, ctxp->fqctx);
    nmod_poly_clear(w);

    return changed;
}

int nmod_mpolyun_interp_mcrt_lg_mpolyu(
    slong * lastdeg,
    nmod_mpolyun_t H,
    const nmod_mpoly_ctx_t ctx,
    nmod_poly_t m,
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ctxp)
{
    slong i;
    int changed = 0;
    fq_nmod_t inv_m_eval;

    *lastdeg = -WORD(1);

    fq_nmod_init(inv_m_eval, ctxp->fqctx);
    nmod_poly_rem(inv_m_eval, m, ctxp->fqctx->modulus);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ctxp->fqctx);

    FLINT_ASSERT(H->bits == A->bits);
    FLINT_ASSERT(H->length == A->length);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(H->exps[i] == A->exps[i]);
        changed |= nmod_mpolyn_CRT_fq_nmod_mpoly(lastdeg, H->coeffs + i, ctx, m,
                                              inv_m_eval, A->coeffs + i, ctxp);
    }
    H->length = A->length;
    fq_nmod_clear(inv_m_eval, ctxp->fqctx);
    return changed;
}

/*****************************************************************************/

/*
    E = A(v = alpha)
    A is in Fq[X][v]
    E is in Fq[X]
*/

void fq_nmod_mpolyn_interp_reduce_sm_poly(
    fq_nmod_poly_t E,
    const fq_nmod_mpolyn_t A,
    const fq_nmod_t alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_t v;
    slong Ai, Alen, k;
    fq_nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fq_nmod_init(v, ctx->fqctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    fq_nmod_poly_zero(E, ctx->fqctx);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fq_nmod_poly_evaluate_fq_nmod(v, Acoeff + Ai, alpha, ctx->fqctx);
        k = (Aexp + N*Ai)[off] >> shift;
        fq_nmod_poly_set_coeff(E, k, v, ctx->fqctx);
    }

    fq_nmod_clear(v, ctx->fqctx);
}

/*
    A = B
    A is in Fq[X][v]  (no v appears)
    B is in Fq[X]
*/

void fq_nmod_mpolyn_interp_lift_sm_poly(
    fq_nmod_mpolyn_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong Bi;
    slong Blen = fq_nmod_poly_length(B, ctx->fqctx);
    fq_nmod_struct * Bcoeff = B->coeffs;
    fq_nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fq_nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bi = Blen - 1; Bi >= 0; Bi--)
    {
        if (!fq_nmod_is_zero(Bcoeff + Bi, ctx->fqctx))
        {
            FLINT_ASSERT(Ai < A->alloc);

            fq_nmod_poly_set_fq_nmod(Acoeff + Ai, Bcoeff + Bi, ctx->fqctx);
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
    F is in Fq[X][v]
    A is in Fq[X]
    it is expected that modulus(alpha) == 1
*/

int fq_nmod_mpolyn_interp_crt_sm_poly(
    slong * lastdeg_,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    const fq_nmod_poly_t A,
    const fq_nmod_poly_t modulus,
    const fq_nmod_t alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    fq_nmod_t u, v;
    slong Fi, Ti, Ai;
    fq_nmod_struct * Acoeff = A->coeffs;
    slong Flen = F->length;
    fq_nmod_poly_struct * Fcoeffs = F->coeffs;
    ulong * Fexps = F->exps;
    fq_nmod_poly_struct * Tcoeffs;
    ulong * Texps;
    fq_nmod_poly_t tp;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    Fi = 0;
    Ai = fq_nmod_poly_degree(A, ctx->fqctx);

    fq_nmod_init(u, ctx->fqctx);
    fq_nmod_init(v, ctx->fqctx);
    fq_nmod_poly_init(tp, ctx->fqctx);

    fq_nmod_mpolyn_fit_length(T, Flen + Ai + 1, ctx);
    Tcoeffs = T->coeffs;
    Texps = T->exps;
    Ti = 0;

    while (Fi < Flen || Ai >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Fcoeffs + Fi, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree(Fcoeffs + Fi, ctx->fqctx) < fq_nmod_poly_degree(modulus, ctx->fqctx));
        }

        if (Ai >= 0)
        {
            FLINT_ASSERT(!fq_nmod_is_zero(Acoeff + Ai, ctx->fqctx));
        }

        mpoly_monomial_zero(Texps + N*Ti, N);

        if (Fi < Flen && Ai >= 0 && ((Fexps + N*Fi)[off]>>shift) == Ai)
        {
            /* F term ok, A term ok */
            fq_nmod_poly_evaluate_fq_nmod(u, Fcoeffs + Fi, alpha, ctx->fqctx);
            fq_nmod_sub(v, Acoeff + Ai, u, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_add(Tcoeffs + Ti, Fcoeffs + Fi, tp, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, ctx->fqctx);
            }

            (Texps + N*Ti)[off] = Ai << shift;
            Fi++;
            do {
                Ai--;
            } while (Ai >= 0 && fq_nmod_is_zero(Acoeff + Ai, ctx->fqctx));
        }
        else if (Fi < Flen && (Ai < 0 || ((Fexps + N*Fi)[off]>>shift) > Ai))
        {
            /* F term ok, A term missing */
            fq_nmod_poly_evaluate_fq_nmod(v, Fcoeffs + Fi, alpha, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_sub(Tcoeffs + Ti, Fcoeffs + Fi, tp, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, ctx->fqctx);
            }

            (Texps + N*Ti)[off] = (Fexps + N*Fi)[off];
            Fi++;
        }
        else if (Ai >= 0 && (Fi >= Flen || ((Fexps + N*Fi)[off]>>shift) < Ai))
        {
            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_poly_scalar_mul_fq_nmod(Tcoeffs + Ti, modulus, Acoeff + Ai, ctx->fqctx);

            (Texps + N*Ti)[off] = Ai << shift;
            do {
                Ai--;
            } while (Ai >= 0 && fq_nmod_is_zero(Acoeff + Ai, ctx->fqctx));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeffs + Ti, ctx->fqctx));
        lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree(Tcoeffs + Ti, ctx->fqctx));

        Ti++;
    }
    T->length = Ti;

    fq_nmod_clear(u, ctx->fqctx);
    fq_nmod_clear(v, ctx->fqctx);
    fq_nmod_poly_clear(tp, ctx->fqctx);

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    *lastdeg_ = lastdeg;
    return changed;
}


/****************************************************************************/

/*
    E = A(x_var = alpha)
    A is in Fq[x_0, ..., x_(var-2), x_(var-1)][x_var]
    E is in Fq[x_0, ..., x_(var-2)][x_(var-1)]
*/
void fq_nmod_mpolyn_interp_reduce_sm_mpolyn(
    fq_nmod_mpolyn_t E,
    fq_nmod_mpolyn_t A,
    slong var,
    fq_nmod_t alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    fq_nmod_t v;
    fq_nmod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fq_nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    fq_nmod_init(v, ctx->fqctx);

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fq_nmod_poly_evaluate_fq_nmod(v, Acoeff + Ai, alpha, ctx->fqctx);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (fq_nmod_is_zero(v, ctx->fqctx))
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            fq_nmod_poly_set_coeff(Ecoeff + Ei - 1, k, v, ctx->fqctx);
        }
        else
        {
            FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

            /* create new */
            if (Ei >= E->alloc)
            {
                fq_nmod_mpolyn_fit_length(E, Ei + 1, ctx);
                Ecoeff = E->coeffs;
                Eexp = E->exps;
            }
            mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
            fq_nmod_poly_zero(Ecoeff + Ei, ctx->fqctx);
            fq_nmod_poly_set_coeff(Ecoeff + Ei, k, v, ctx->fqctx);
            Ei++;
        }
    }
    E->length = Ei;

    fq_nmod_clear(v, ctx->fqctx);
}


/*
    A = B
    A is in Fq[x_0, ..., x_(var-1), x_(var-1)][x_var]
    B is in Fq[x_0, ..., x_(var-2)][x_(var-1)]
*/
void fq_nmod_mpolyn_interp_lift_sm_mpolyn(
    fq_nmod_mpolyn_t A,
    fq_nmod_mpolyn_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;

    fq_nmod_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;

    fq_nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;

    fq_nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Ai = 0;
    for (Bi = 0; Bi < Blen; Bi++)
    {
        if (Ai + (Bcoeff + Bi)->length >= A->alloc)
        {
            fq_nmod_mpolyn_fit_length(A, Ai + (Bcoeff + Bi)->length, ctx);
            Acoeff = A->coeffs;
            Aexp = A->exps;
        }
        for (vi = (Bcoeff + Bi)->length - 1; vi >= 0; vi--)
        {
            if (!fq_nmod_is_zero((Bcoeff + Bi)->coeffs + vi, ctx->fqctx))
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                fq_nmod_poly_zero(Acoeff + Ai, ctx->fqctx);
                fq_nmod_poly_set_coeff(Acoeff + Ai, 0, (Bcoeff + Bi)->coeffs + vi, ctx->fqctx);
                Ai++;
            }
        }
    }
    A->length = Ai;
}


/*
    T = F + modulus*(A - F(x_var = alpha))
    no assumptions about matching monomials
    F is in Fq[x_0, ..., x_(var-1), x_(var-1)][x_var]
    A is in Fq[x_0, ..., x_(var-2)][x_(var-1)]
    in order to fxn correctly, modulus(alpha) should be 1
*/
int fq_nmod_mpolyn_interp_crt_sm_mpolyn(
    slong * lastdeg_,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    fq_nmod_mpolyn_t A,
    slong var,
    fq_nmod_poly_t modulus,
    const fq_nmod_t alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    slong vi;
    fq_nmod_t v;
    fq_nmod_poly_t tp;

    fq_nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    fq_nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    fq_nmod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    fq_nmod_poly_init(tp, ctx->fqctx);
    fq_nmod_init(v, ctx->fqctx);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    fq_nmod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    Fi = Ai = vi = 0;
    if (Ai < Alen)
    {
        vi = fq_nmod_poly_degree(A->coeffs + Ai, ctx->fqctx);
    }
    while (Fi < Flen || Ai < Alen)
    {
        if (Ti >= T->alloc)
        {
            fq_nmod_mpolyn_fit_length(T, Ti + FLINT_MAX(Flen - Fi, Alen - Ai), ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        if (Ai < Alen)
        {
            FLINT_ASSERT(!fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ctx->fqctx));
        }

        if (Fi < Flen && Ai < Alen && mpoly_monomial_equal_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift))
        {
            /* F term ok, A term ok */
            fq_nmod_poly_evaluate_fq_nmod(v, Fcoeff + Fi, alpha, ctx->fqctx);
            fq_nmod_sub(v, (Acoeff + Ai)->coeffs + vi, v, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_add(Tcoeff + Ti, Fcoeff + Fi, tp, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ctx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ctx->fqctx);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            fq_nmod_poly_evaluate_fq_nmod(v, Fcoeff + Fi, alpha, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, tp, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_poly_scalar_mul_fq_nmod(Tcoeff + Ti, modulus, (Acoeff + Ai)->coeffs + vi, ctx->fqctx);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ctx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ctx->fqctx);
                }
            }
        }

        FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + Ti, ctx->fqctx));
        lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree(Tcoeff + Ti, ctx->fqctx));
        Ti++;
    }
    T->length = Ti;

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    fq_nmod_poly_clear(tp, ctx->fqctx);
    fq_nmod_clear(v, ctx->fqctx);

    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(F, ctx));

    *lastdeg_ = lastdeg;
    return changed;
}

/*****************************************************************************/

/*
    A is in              Fq [X][v]
    E is in (Fq[v]/alpha(v))[X]
*/
void fq_nmod_mpolyn_interp_reduce_lg_poly(
    fq_nmod_poly_t E,
    const fq_nmod_mpoly_ctx_t ectx,
    fq_nmod_mpolyn_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    const bad_fq_nmod_embed_t emb)
{
    fq_nmod_t v;
    slong Ai, Alen, k;
    fq_nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fq_nmod_init(v, ectx->fqctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;

    fq_nmod_poly_zero(E, ectx->fqctx);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        k = (Aexp + N*Ai)[off] >> shift;
        bad_fq_nmod_embed_sm_to_lg(v, Acoeff + Ai, emb);
        fq_nmod_poly_set_coeff(E, k, v, ectx->fqctx);
    }

    fq_nmod_clear(v, ectx->fqctx);
}


/*
    A = B
    A is in              Fq [X][][v]
    B is in (Fq[v]/alpha(v))[X]
*/
void fq_nmod_mpolyn_interp_lift_lg_poly(
    slong * lastdeg_,
    fq_nmod_mpolyn_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_poly_t B,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    slong Bexp;
    slong Blen = fq_nmod_poly_length(B, ectx->fqctx);
    fq_nmod_struct * Bcoeff = B->coeffs;
    fq_nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong lastdeg = -WORD(1);
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fq_nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bexp = Blen - 1; Bexp >= 0; Bexp--)
    {
        if (!fq_nmod_is_zero(Bcoeff + Bexp, ectx->fqctx))
        {
            FLINT_ASSERT(Ai < A->alloc);
            bad_fq_nmod_embed_lg_to_sm(Acoeff + Ai, Bcoeff + Bexp, emb);
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Acoeff + Ai, ctx->fqctx));
            lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree(Acoeff + Ai, ctx->fqctx));
            mpoly_monomial_zero(Aexp + N*Ai, N);
            (Aexp + N*Ai)[off] = Bexp << shift;
            Ai++;
        }
    }
    A->length = Ai;

    *lastdeg_ = lastdeg;
}


/*
    F = F + modulus*((A - F(alpha))/(modulus(alpha)))
    F is in              Fq [X][][v]
    A is in (Fq[v]/alpha(v))[X]
    alpha(v) is ectx->fqctx->modulus(v)
*/
int fq_nmod_mpolyn_interp_crt_lg_poly(
    slong * lastdeg_,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    fq_nmod_poly_t modulus,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_poly_t A,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    fq_nmod_t u, v;
    fq_nmod_poly_t u_sm, w;
    slong Fi, Ti, Aexp;
    fq_nmod_struct * Acoeff = A->coeffs;
    slong Flen = F->length;
    fq_nmod_poly_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    fq_nmod_poly_struct * Tcoeff;
    ulong * Texp;
    fq_nmod_poly_t tp;
    fq_nmod_t inv_m_eval;
    slong N, off, shift;

    FLINT_ASSERT(T->bits == F->bits);

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    bad_fq_nmod_embed_sm_to_lg(inv_m_eval, modulus, emb);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    fq_nmod_poly_init(u_sm, ctx->fqctx);
    fq_nmod_poly_init(w, ctx->fqctx);

    Fi = 0;

    Aexp = fq_nmod_poly_degree(A, ectx->fqctx);

    fq_nmod_poly_init(tp, ctx->fqctx);

    fq_nmod_mpolyn_fit_length(T, Flen + Aexp + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    while (Fi < Flen || Aexp >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Fcoeff + Fi, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree(Fcoeff + Fi, ctx->fqctx) < fq_nmod_poly_degree(modulus, ctx->fqctx));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(!fq_nmod_is_zero(Acoeff + Aexp, ectx->fqctx));
        }

        mpoly_monomial_zero(Texp + N*Ti, N);

        if (Fi < Flen && Aexp >= 0 && ((Fexp + N*Fi)[off]>>shift) == Aexp)
        {
            /* F term ok, A term ok */
            bad_fq_nmod_embed_sm_to_lg(u, Fcoeff + Fi, emb);
            fq_nmod_sub(v, Acoeff + Aexp, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, modulus, u_sm, ctx->fqctx);
                fq_nmod_poly_add(Tcoeff + Ti, Fcoeff + Fi, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }

            (Texp + N*Ti)[off] = (Fexp + N*Fi)[off];
            Fi++;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, ectx->fqctx));
        }
        else if (Fi < Flen && (Aexp < 0 || ((Fexp + N*Fi)[off]>>shift) > Aexp))
        {
            /* F term ok, A term missing */
            bad_fq_nmod_embed_sm_to_lg(v, Fcoeff + Fi, emb);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, u_sm, modulus, ctx->fqctx);
                fq_nmod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }

            (Texp + N*Ti)[off] = (Fexp + N*Fi)[off];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || ((Fexp + N*Fi)[off]>>shift) < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_mul(u, Acoeff + Aexp, inv_m_eval, ectx->fqctx);
            bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
            fq_nmod_poly_mul(Tcoeff + Ti, modulus, u_sm, ctx->fqctx);

            (Texp + N*Ti)[off] = Aexp << shift;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, ectx->fqctx));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + Ti, ctx->fqctx));
        lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree(Tcoeff + Ti, ctx->fqctx));
        Ti++;
    }
    T->length = Ti;

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_poly_clear(u_sm, ctx->fqctx);
    fq_nmod_poly_clear(w, ctx->fqctx);

    fq_nmod_clear(inv_m_eval, ectx->fqctx);

    *lastdeg_ = lastdeg;
    return changed;
}



/*****************************************************************************/

/*
    A is in                      Fq [x_0,...,x_(var-2), x_(var-1)][x_var]
    E is in (Fq[x_var]/alpha(x_var))[x_0,...,x_(var-2)][x_(var-1)]
*/

void fq_nmod_mpolyn_interp_reduce_lg_mpolyn(
    fq_nmod_mpolyn_t E,
    const fq_nmod_mpoly_ctx_t ectx,
    fq_nmod_mpolyn_t A,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx,
    const bad_fq_nmod_embed_t emb)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    fq_nmod_t v;
    fq_nmod_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    fq_nmod_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;

    fq_nmod_init(v, ectx->fqctx);

    FLINT_ASSERT(var > 0);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        bad_fq_nmod_embed_sm_to_lg(v, Acoeff + Ai, emb);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (fq_nmod_is_zero(v, ectx->fqctx))
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            fq_nmod_poly_set_coeff(Ecoeff + Ei - 1, k, v, ectx->fqctx);
        }
        else
        {
            FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

            /* create new */
            if (Ei >= E->alloc)
            {
                fq_nmod_mpolyn_fit_length(E, Ei + 1, ectx);
                Ecoeff = E->coeffs;
                Eexp = E->exps;
            }
            mpoly_monomial_set_extra(Eexp + N*Ei, Aexp + N*Ai, N, offset, -(k << shift));
            fq_nmod_poly_zero(Ecoeff + Ei, ectx->fqctx);
            fq_nmod_poly_set_coeff(Ecoeff + Ei, k, v, ectx->fqctx);
            Ei++;
        }
    }
    E->length = Ei;

    fq_nmod_clear(v, ectx->fqctx);
}


/*
    A = B using lowest degree representative
    A is in                      Fq [x_0,...,x_(var-2), x_(var-1)][x_var]
    B is in (Fq[x_var]/alpha(x_var))[x_0,...,x_(var-2)][x_(var-1)]
    alpha is emb->h
*/
void fq_nmod_mpolyn_interp_lift_lg_mpolyn(
    slong * lastdeg_,
    fq_nmod_mpolyn_t A,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_mpolyn_t B,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;
    fq_nmod_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;
    fq_nmod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong lastdeg = -WORD(1);

    FLINT_ASSERT(A->bits == B->bits);

    fq_nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Ai = 0;
    for (Bi = 0; Bi < Blen; Bi++)
    {
        if (Ai + (Bcoeff + Bi)->length >= A->alloc)
        {
            fq_nmod_mpolyn_fit_length(A, Ai + (Bcoeff + Bi)->length, ctx);
            Acoeff = A->coeffs;
            Aexp = A->exps;
        }
        for (vi = (Bcoeff + Bi)->length - 1; vi >= 0; vi--)
        {
            if (!fq_nmod_is_zero((Bcoeff + Bi)->coeffs + vi, ectx->fqctx))
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                bad_fq_nmod_embed_lg_to_sm(Acoeff + Ai, (Bcoeff + Bi)->coeffs + vi, emb);
                FLINT_ASSERT(!fq_nmod_poly_is_zero(Acoeff + Ai, ctx->fqctx));
                lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree(Acoeff + Ai, ctx->fqctx));
                Ai++;
            }
        }
    }
    A->length = Ai;

    *lastdeg_ = lastdeg;
}


/*
    F = F + modulus*((A - F(alpha))/modulus(alpha))
    F is in                      Fq [x_0,...,x_(var-2), x_(var-1)][x_var]
    A is in (Fq[x_var]/alpha(x_var))[x_0,...,x_(var-2)][x_(var-1)]
*/
int fq_nmod_mpolyn_interp_crt_lg_mpolyn(
    slong * lastdeg_,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    fq_nmod_poly_t modulus,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_mpolyn_t A,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;
    slong lastdeg = -WORD(1);
    slong vi;
    fq_nmod_t u, v;
    fq_nmod_poly_t u_sm, w;
    fq_nmod_t inv_m_eval;

    fq_nmod_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    fq_nmod_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    fq_nmod_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    bad_fq_nmod_embed_sm_to_lg(inv_m_eval, modulus, emb);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    fq_nmod_poly_init(u_sm, ctx->fqctx);
    fq_nmod_poly_init(w, ctx->fqctx);

    FLINT_ASSERT(var > 0);
    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Flen = F->length;

    fq_nmod_mpolyn_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    Fi = Ai = vi = 0;
    if (Ai < Alen)
    {
        vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
    }
    while (Fi < Flen || Ai < Alen)
    {
        if (Ti >= T->alloc)
        {
            fq_nmod_mpolyn_fit_length(T, Ti + FLINT_MAX(Flen - Fi, Alen - Ai), ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        if (Fi < Flen && Ai < Alen && mpoly_monomial_equal_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift))
        {
            /* F term ok, A term ok */
            bad_fq_nmod_embed_sm_to_lg(u, Fcoeff + Fi, emb);
            fq_nmod_sub(v, (Acoeff + Ai)->coeffs + vi, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, modulus, u_sm, ctx->fqctx);
                fq_nmod_poly_add(Tcoeff + Ti, Fcoeff + Fi, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ectx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            bad_fq_nmod_embed_sm_to_lg(v, Fcoeff + Fi, emb);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, u_sm, modulus, ctx->fqctx);
                fq_nmod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_mul(u, (Acoeff + Ai)->coeffs + vi, inv_m_eval, ectx->fqctx);
            bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
            fq_nmod_poly_mul(Tcoeff + Ti, modulus, u_sm, ctx->fqctx);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && fq_nmod_is_zero((Acoeff + Ai)->coeffs + vi, ectx->fqctx));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = fq_nmod_poly_degree(A->coeffs + Ai, ectx->fqctx);
                }
            }
        }

        FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + Ti, ctx->fqctx));
        lastdeg = FLINT_MAX(lastdeg, fq_nmod_poly_degree(Tcoeff + Ti, ctx->fqctx));
        Ti++;
    }
    T->length = Ti;

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    fq_nmod_clear(inv_m_eval, ectx->fqctx);
    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_poly_clear(u_sm, ctx->fqctx);
    fq_nmod_poly_clear(w, ctx->fqctx);

    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(F, ctx));

    *lastdeg_ = lastdeg;
    return changed;
}



/*****************************************************************************/

/* evaluate A at lastvar = alpha */
void fq_nmod_mpolyn_interp_reduce_sm_mpoly(
    fq_nmod_mpoly_t B,
    fq_nmod_mpolyn_t A,
    fq_nmod_t alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, k, N;
    FLINT_ASSERT(B->bits == A->bits);

    fq_nmod_mpoly_fit_length(B, A->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    k = 0;
    for (i = 0; i < A->length; i++)
    {
        mpoly_monomial_set(B->exps + N*k, A->exps + N*i, N);
        fq_nmod_poly_evaluate_fq_nmod(B->coeffs + k, A->coeffs + i, alpha,
                                                                   ctx->fqctx);
        k += !fq_nmod_is_zero(B->coeffs + k, ctx->fqctx);
    }
    B->length = k;
}

void fq_nmod_mpolyun_interp_reduce_sm_mpolyu(
    fq_nmod_mpolyu_t B,
    fq_nmod_mpolyun_t A,
    fq_nmod_t alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, k;

    FLINT_ASSERT(B->bits == A->bits);

    fq_nmod_mpolyu_fit_length(B, A->length, ctx);
    k = 0;
    for (i = 0; i < A->length; i++)
    {
        B->exps[k] = A->exps[i];
        fq_nmod_mpolyn_interp_reduce_sm_mpoly(B->coeffs + k, A->coeffs + i, alpha, ctx);
        k += !fq_nmod_mpoly_is_zero(B->coeffs + k, ctx);
    }
    B->length = k;
}


void fq_nmod_mpolyn_interp_lift_sm_mpoly(
    fq_nmod_mpolyn_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong N;
    fq_nmod_poly_struct * Acoeff;
    fq_nmod_struct * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Blen;

    fq_nmod_mpolyn_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    Blen = B->length;
    fq_nmod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < Blen; i++)
    {
        fq_nmod_poly_zero(Acoeff + i, ctx->fqctx);
        fq_nmod_poly_set_coeff(Acoeff + i, 0, Bcoeff + i, ctx->fqctx);
        mpoly_monomial_set(Aexp + N*i, Bexp + N*i, N);
    }
    A->length = Blen;
}

void fq_nmod_mpolyun_interp_lift_sm_mpolyu(
    fq_nmod_mpolyun_t A,
    const fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(A->bits == B->bits);
    fq_nmod_mpolyun_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fq_nmod_mpolyn_interp_lift_sm_mpoly(A->coeffs + i, B->coeffs + i, ctx);
        FLINT_ASSERT((A->coeffs + i)->bits == B->bits);
    }
    A->length = B->length;
}


/*
    F = F + modulus*(A - F(alpha))
*/
int fq_nmod_mpolyn_interp_crt_sm_mpoly(
    slong * lastdeg,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    fq_nmod_mpoly_t A,
    fq_nmod_poly_t modulus,
    fq_nmod_t alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    slong N;
    fq_nmod_t v;
    flint_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    fq_nmod_struct * Acoeff = A->coeffs;
    fq_nmod_poly_struct * Fcoeff = F->coeffs;
    fq_nmod_poly_struct * Tcoeff;
    fq_nmod_poly_t tp;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fq_nmod_init(v, ctx->fqctx);
    fq_nmod_poly_init(tp, ctx->fqctx);

    fq_nmod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen
                        || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Fcoeff + i, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree(Fcoeff + i, ctx->fqctx)
                                   < fq_nmod_poly_degree(modulus, ctx->fqctx));

            /* F term ok, A term missing */
            fq_nmod_poly_evaluate_fq_nmod(v, Fcoeff + i, alpha, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_sub(Tcoeff + k, Fcoeff + i, tp, ctx->fqctx);
            } else {
                fq_nmod_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);                
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0],
                                  fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + k, ctx->fqctx));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (!fq_nmod_is_zero(Acoeff + j, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_zero(Tcoeff + k, ctx->fqctx);
                fq_nmod_poly_scalar_mul_fq_nmod(Tcoeff + k, modulus, Acoeff + j,
                                                                   ctx->fqctx);
                lastdeg[0] = FLINT_MAX(lastdeg[0],
                                  fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Fcoeff + i, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree(Fcoeff + i, ctx->fqctx) < fq_nmod_poly_degree(modulus, ctx->fqctx));

            /* F term ok, A term ok */
            fq_nmod_poly_evaluate_fq_nmod(v, Fcoeff + i, alpha, ctx->fqctx);
            fq_nmod_sub(v, Acoeff + j, v, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                fq_nmod_poly_add(Tcoeff + k, Fcoeff + i, tp, ctx->fqctx);
            } else {
                fq_nmod_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0],
                                  fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + k, ctx->fqctx));
            k++;
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }
    T->length = k;

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    fq_nmod_poly_clear(tp, ctx->fqctx);
    fq_nmod_clear(v, ctx->fqctx);

    return changed;
}


int fq_nmod_mpolyun_interp_crt_sm_mpolyu(
    slong * lastdeg,
    fq_nmod_mpolyun_t F,
    fq_nmod_mpolyun_t T,
    fq_nmod_mpolyu_t A,
    fq_nmod_poly_t modulus,
    fq_nmod_t alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    fq_nmod_mpolyn_t S;
    fq_nmod_mpolyn_struct * Tcoeff;
    fq_nmod_mpolyn_struct * Fcoeff;
    fq_nmod_mpoly_struct  * Acoeff;
    fq_nmod_mpoly_t zero;

    lastdeg[0] = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    fq_nmod_mpolyn_init(S, F->bits, ctx);

    Flen = F->length;
    Alen = A->length;
    fq_nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    fq_nmod_mpoly_init(zero, ctx);
    fq_nmod_mpoly_fit_bits(zero, A->bits, ctx);
    zero->bits = A->bits;

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            fq_nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= fq_nmod_mpolyn_interp_crt_sm_mpoly(lastdeg, Tcoeff + k, S, zero, modulus, alpha, ctx);
            Texp[k] = Fexp[i];
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || Aexp[j] > Fexp[i]))
        {
            /* F term missing, A term ok */
            fq_nmod_mpolyn_zero(Tcoeff + k, ctx);
            changed |= fq_nmod_mpolyn_interp_crt_sm_mpoly(lastdeg, Tcoeff + k, S, Acoeff + j, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            k++;
            j++;
        }
        else if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            fq_nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= fq_nmod_mpolyn_interp_crt_sm_mpoly(lastdeg, Tcoeff + k, S, Acoeff + j, modulus, alpha, ctx);
            Texp[k] = Aexp[j];
            FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            i++;
            j++;
        } else 
        {
            FLINT_ASSERT(0);
        }
    }
    T->length = k;

    if (changed)
    {
        fq_nmod_mpolyun_swap(T, F);
    }

    fq_nmod_mpolyn_clear(S, ctx);
    fq_nmod_mpoly_clear(zero, ctx);
    return changed;    
}


/*****************************************************************************/

/* reduce B via the map F_q[x] -> F_q^n */
void fq_nmod_mpolyn_interp_reduce_lg_mpoly(
    fq_nmod_mpoly_t A,
    fq_nmod_mpolyn_t B,
    const fq_nmod_mpoly_ctx_t ectx,
    const fq_nmod_mpoly_ctx_t ctx,
    const bad_fq_nmod_embed_t emb)
{
    slong i, k, N;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ectx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ectx->minfo->nvars == ctx->minfo->nvars);

    N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);

    k = 0;
    fq_nmod_mpoly_fit_length(A, k + 1, ectx);
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_fit_length(A, k + 1, ectx);
        mpoly_monomial_set(A->exps + N*k, B->exps + N*i, N);
        bad_fq_nmod_embed_sm_to_lg(A->coeffs + k, B->coeffs + i, emb);
        k += !fq_nmod_is_zero(A->coeffs + k, ectx->fqctx);
    }

    _fq_nmod_mpoly_set_length(A, k, ectx);
}

/* reduce B via the map from F_q[x] -> F_q^n */
void fq_nmod_mpolyun_interp_reduce_lg_mpolyu(
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyun_t B,
    const fq_nmod_mpoly_ctx_t ectx,
    const fq_nmod_mpoly_ctx_t ctx,
    const bad_fq_nmod_embed_t emb)
{
    slong i, k, Blen;
    fq_nmod_mpoly_struct * Acoeff;
    fq_nmod_mpolyn_struct * Bcoeff;
    ulong * Aexp, * Bexp;

    Blen = B->length;
    fq_nmod_mpolyu_fit_length(A, Blen, ectx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    k = 0;
    for (i = 0; i < Blen; i++)
    {
        fq_nmod_mpolyn_interp_reduce_lg_mpoly(Acoeff + k, Bcoeff + i, ectx, ctx, emb);
        Aexp[k] = Bexp[i];
        k += !fq_nmod_mpoly_is_zero(Acoeff + k, ectx);
    }
    A->length = k;  
}

/* Convert B to A using the map  F_q^n -> F_q[x] */
void fq_nmod_mpolyn_interp_lift_lg_mpoly(
    fq_nmod_mpolyn_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    slong i, N;

    FLINT_ASSERT(B->bits == A->bits);
    fq_nmod_mpolyn_fit_length(A, B->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
        bad_fq_nmod_embed_lg_to_sm(A->coeffs + i, B->coeffs + i, emb);
        FLINT_ASSERT(!fq_nmod_poly_is_zero(A->coeffs + i, ctx->fqctx));
    }
    A->length = B->length;
}

void fq_nmod_mpolyun_interp_lift_lg_mpolyu(
    fq_nmod_mpolyun_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    slong i;

    FLINT_ASSERT(B->bits == A->bits);
    fq_nmod_mpolyun_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fq_nmod_mpolyn_interp_lift_lg_mpoly(A->coeffs + i, ctx, B->coeffs + i, ectx, emb);
        FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(A->coeffs + i, ctx));
    }
    A->length = B->length;
}


/*
    update F so that it doesn't change mod m and is A mod emb->h

    F = F + m*((A - F(alpha))/(m(alpha)))

    no assumptions about matching monomials
*/
int fq_nmod_mpolyn_interp_crt_lg_mpoly(
    slong * lastdeg,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    fq_nmod_poly_t m,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_mpoly_t A,
    fq_nmod_t inv_m_eval,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    int changed = 0;
    slong i, j, k;
    slong N;
    fq_nmod_t u, v;
    fq_nmod_poly_t u_sm, w;
    flint_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    fq_nmod_struct * Acoeff = A->coeffs;
    fq_nmod_poly_struct * Fcoeff = F->coeffs;
    fq_nmod_poly_struct * Tcoeff;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    fq_nmod_poly_init(u_sm, ctx->fqctx);
    fq_nmod_poly_init(w, ctx->fqctx);

    fq_nmod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen
                        || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Fcoeff + i, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree(Fcoeff + i, ctx->fqctx)
                                  < fq_nmod_poly_degree(m, ctx->fqctx));

            /* F term ok, A term missing */
            bad_fq_nmod_embed_sm_to_lg(v, Fcoeff + i, emb);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, u_sm, m, ctx->fqctx);
                fq_nmod_poly_sub(Tcoeff + k, Fcoeff + i, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + k, ctx->fqctx));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (!fq_nmod_is_zero(Acoeff + j, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, Acoeff + j, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(Tcoeff + k, m, u_sm, ctx->fqctx);
                lastdeg[0] = FLINT_MAX(lastdeg[0], fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Fcoeff + i, ctx->fqctx));
            FLINT_ASSERT(fq_nmod_poly_degree(Fcoeff + i, ctx->fqctx)
                                        < fq_nmod_poly_degree(m, ctx->fqctx));

            /* F term ok, A term ok */
            bad_fq_nmod_embed_sm_to_lg(u, Fcoeff + i, emb);
            fq_nmod_sub(v, Acoeff + j, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, m, u_sm, ctx->fqctx);
                fq_nmod_poly_add(Tcoeff + k, Fcoeff + i, w, ctx->fqctx);
            }
            else
            {
                fq_nmod_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);
            }
            lastdeg[0] = FLINT_MAX(lastdeg[0], fq_nmod_poly_degree(Tcoeff + k, ctx->fqctx));
            mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
            FLINT_ASSERT(!fq_nmod_poly_is_zero(Tcoeff + k, ctx->fqctx));
            k++;
            i++;
            j++;
        }
        else 
        {
            FLINT_ASSERT(0);
        }
    }

    fq_nmod_mpolyn_set_length(T, k, ctx);

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_poly_clear(u_sm, ctx->fqctx);
    fq_nmod_poly_clear(w, ctx->fqctx);

    return changed;
}


int fq_nmod_mpolyun_interp_crt_lg_mpolyu(
    slong * lastdeg,
    fq_nmod_mpolyun_t F,
    fq_nmod_mpolyun_t T,
    fq_nmod_poly_t m,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    int changed = 0;
    slong i, j, k;
    ulong * Texp;
    ulong * Fexp;
    ulong * Aexp;
    slong Flen;
    slong Alen;
    fq_nmod_mpolyn_t S;
    fq_nmod_mpolyn_struct * Tcoeff;
    fq_nmod_mpolyn_struct * Fcoeff;
    fq_nmod_mpoly_struct  * Acoeff;
    fq_nmod_mpoly_t zero;
    fq_nmod_t inv_m_eval;

    lastdeg[0] = -WORD(1);

    FLINT_ASSERT(F->bits == T->bits);
    FLINT_ASSERT(T->bits == A->bits);

    fq_nmod_mpolyn_init(S, F->bits, ctx);

    Flen = F->length;
    Alen = A->length;
    fq_nmod_mpolyun_fit_length(T, Flen + Alen, ctx);

    Tcoeff = T->coeffs;
    Fcoeff = F->coeffs;
    Acoeff = A->coeffs;
    Texp = T->exps;
    Fexp = F->exps;
    Aexp = A->exps;   

    fq_nmod_mpoly_init(zero, ectx);
    fq_nmod_mpoly_fit_bits(zero, A->bits, ectx);
    zero->bits = A->bits;

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    bad_fq_nmod_embed_sm_to_lg(inv_m_eval, m, emb);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen || Fexp[i] > Aexp[j]))
        {
            /* F term ok, A term missing */
            fq_nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= fq_nmod_mpolyn_interp_crt_lg_mpoly(lastdeg, Tcoeff + k,
                                       S, m, ctx, zero, inv_m_eval, ectx, emb);
            Texp[k] = Fexp[i];
            FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen || Aexp[j] > Fexp[i]))
        {
            /* F term missing, A term ok */
            fq_nmod_mpolyn_zero(Tcoeff + k, ctx);
            changed |= fq_nmod_mpolyn_interp_crt_lg_mpoly(lastdeg, Tcoeff + k,
                                 S, m, ctx, Acoeff + j, inv_m_eval, ectx, emb);
            Texp[k] = Aexp[j];
            FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(Tcoeff + k, ctx));
            k++;
            j++;
        }
        else if (i < Flen && j < Alen && (Fexp[i] == Aexp[j]))
        {
            /* F term ok, A term ok */
            fq_nmod_mpolyn_set(Tcoeff + k, Fcoeff + i, ctx);
            changed |= fq_nmod_mpolyn_interp_crt_lg_mpoly(lastdeg, Tcoeff + k,
                                 S, m, ctx, Acoeff + j, inv_m_eval, ectx, emb);
            Texp[k] = Aexp[j];
            FLINT_ASSERT(!fq_nmod_mpolyn_is_zero(Tcoeff + k, ctx));
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
        fq_nmod_mpolyun_swap(T, F);
    }

    fq_nmod_clear(inv_m_eval, ectx->fqctx);

    fq_nmod_mpolyn_clear(S, ctx);
    fq_nmod_mpoly_clear(zero, ectx);
    return changed;    
}

/*
    Update H so that it does not change mod m, and is now A mod p
    It is asserted that the monomials in H and A match
*/
int fq_nmod_mpolyn_interp_mcrt_lg_mpoly(
    slong * lastdeg,
    fq_nmod_mpolyn_t H,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_poly_t m,
    fq_nmod_t inv_m_eval,
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    slong i;
#if WANT_ASSERT
    slong N;
#endif
    int changed = 0;
    fq_nmod_t u, v;
    fq_nmod_poly_t w, u_sm;

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    fq_nmod_poly_init(w, ctx->fqctx);
    fq_nmod_poly_init(u_sm, ctx->fqctx);

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);
#if WANT_ASSERT
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
#endif
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        bad_fq_nmod_embed_sm_to_lg(u, H->coeffs + i, emb);
        fq_nmod_sub(v, A->coeffs + i, u, ectx->fqctx);
        if (!fq_nmod_is_zero(v, ectx->fqctx))
        {
            changed = 1;
            fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
            bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
            fq_nmod_poly_mul(w, u_sm, m, ctx->fqctx);
            fq_nmod_poly_add(H->coeffs + i, H->coeffs + i, w, ctx->fqctx);
        }

        *lastdeg = FLINT_MAX(*lastdeg,
                               fq_nmod_poly_degree(H->coeffs + i, ctx->fqctx));

        FLINT_ASSERT(fq_nmod_poly_degree(H->coeffs + i, ctx->fqctx)
                         <  fq_nmod_poly_degree(m, ctx->fqctx)
                          + fq_nmod_poly_degree(emb->h, ctx->fqctx));
    }

    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_poly_clear(w, ctx->fqctx);
    fq_nmod_poly_clear(u_sm, ctx->fqctx);

    return changed;
}

int fq_nmod_mpolyun_interp_mcrt_lg_mpolyu(
    slong * lastdeg,
    fq_nmod_mpolyun_t H,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_poly_t m,
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ectx,
    bad_fq_nmod_embed_t emb)
{
    slong i;
    int changed = 0;
    fq_nmod_t inv_m_eval;

    *lastdeg = -WORD(1);

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    bad_fq_nmod_embed_sm_to_lg(inv_m_eval, m, emb);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    FLINT_ASSERT(H->bits == A->bits);
    FLINT_ASSERT(H->length == A->length);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(H->exps[i] == A->exps[i]);
        changed |= fq_nmod_mpolyn_interp_mcrt_lg_mpoly(lastdeg, H->coeffs + i,
                                 ctx, m, inv_m_eval, A->coeffs + i, ectx, emb);
    }
    H->length = A->length;
    fq_nmod_clear(inv_m_eval, ectx->fqctx);
    return changed;
}


