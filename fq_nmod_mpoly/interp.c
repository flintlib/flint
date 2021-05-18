/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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

static void _n_poly_mul_n_fq(
    n_poly_t a,
    const n_poly_t b,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    n_poly_t C;
    C->coeffs = (mp_limb_t *) c;
    C->length = d;
    C->alloc = d;
    _n_poly_normalise(C);
    n_poly_mod_mul(a, b, C, ctx->mod);
}

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
    n_poly_struct * Acoeff;
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
        n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(v), Acoeff + Ai,
                evil_const_cast_nmod_poly_to_n_poly(fqctx->modulus), ctx->mod);
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
    n_poly_struct * Acoeff;
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
            n_poly_set_nmod_poly(Acoeff + Ai, Bcoeff + Bexp);
            lastdeg = FLINT_MAX(lastdeg, n_poly_degree(Acoeff + Ai));
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
    n_poly_t modulus,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_poly_t A,
    const fq_nmod_ctx_t fqctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    fq_nmod_t u, v;
    n_poly_t w;
    slong Fi, Ti, Aexp;
    fq_nmod_struct * Acoeff = A->coeffs;
    slong Flen = F->length;
    n_poly_struct * Fcoeffs = F->coeffs;
    ulong * Fexps = F->exps;
    n_poly_struct * Tcoeffs;
    ulong * Texps;
    fq_nmod_t inv_m_eval;
    slong N, off, shift;

    FLINT_ASSERT(T->bits == F->bits);

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    fq_nmod_init(inv_m_eval, fqctx);
    n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(inv_m_eval), modulus,
                evil_const_cast_nmod_poly_to_n_poly(fqctx->modulus), ctx->mod);
    fq_nmod_inv(inv_m_eval, inv_m_eval, fqctx);

    Fi = 0;

    fq_nmod_init(u, fqctx);
    fq_nmod_init(v, fqctx);
    n_poly_init(w);

    Aexp = fq_nmod_poly_degree(A, fqctx);

    nmod_mpolyn_fit_length(T, Flen + Aexp + 1, ctx);
    Tcoeffs = T->coeffs;
    Texps = T->exps;
    Ti = 0;

    while (Fi < Flen || Aexp >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!n_poly_is_zero(Fcoeffs + Fi));
            FLINT_ASSERT(n_poly_degree(Fcoeffs + Fi) < n_poly_degree(modulus));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(!fq_nmod_is_zero(Acoeff + Aexp, fqctx));
        }

        mpoly_monomial_zero(Texps + N*Ti, N);

        if (Fi < Flen && Aexp >= 0 && ((Fexps + N*Fi)[off]>>shift) == Aexp)
        {
            /* F term ok, A term ok */
            n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(u), Fcoeffs + Fi,
                evil_const_cast_nmod_poly_to_n_poly(fqctx->modulus), ctx->mod);
            fq_nmod_sub(v, Acoeff + Aexp, u, fqctx);
            if (!fq_nmod_is_zero(v, fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, fqctx);
                n_poly_mod_mul(w, modulus, evil_cast_nmod_poly_to_n_poly(u), ctx->mod);
                n_poly_mod_add(Tcoeffs + Ti, Fcoeffs + Fi, w, ctx->mod);
            }
            else
            {
                n_poly_set(Tcoeffs + Ti, Fcoeffs + Fi);
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
            n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(v), Fcoeffs + Fi,
                evil_const_cast_nmod_poly_to_n_poly(fqctx->modulus), ctx->mod);
            if (!fq_nmod_is_zero(v, fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, fqctx);
                n_poly_mod_mul(w, evil_const_cast_nmod_poly_to_n_poly(u), modulus, ctx->mod);
                n_poly_mod_sub(Tcoeffs + Ti, Fcoeffs + Fi, w, ctx->mod);
            }
            else
            {
                n_poly_set(Tcoeffs + Ti, Fcoeffs + Fi);
            }

            (Texps + N*Ti)[off] = (Fexps + N*Fi)[off];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || ((Fexps + N*Fi)[off]>>shift) < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_mul(u, Acoeff + Aexp, inv_m_eval, fqctx);
            n_poly_mod_mul(Tcoeffs + Ti, modulus,
                             evil_const_cast_nmod_poly_to_n_poly(u), ctx->mod);
            (Texps + N*Ti)[off] = Aexp << shift;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, fqctx));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        FLINT_ASSERT(!n_poly_is_zero(Tcoeffs + Ti));
        lastdeg = FLINT_MAX(lastdeg, n_poly_degree(Tcoeffs + Ti));
        Ti++;
    }
    T->length = Ti;

    if (changed)
    {
        nmod_mpolyn_swap(T, F);
    }

    fq_nmod_clear(u, fqctx);
    fq_nmod_clear(v, fqctx);
    n_poly_clear(w);

    fq_nmod_clear(inv_m_eval, fqctx);

    *lastdeg_ = lastdeg;
    return changed;
}

void nmod_mpolyn_interp_lift_lg_bpoly(
    slong * lastdeg_,
    nmod_mpolyn_t F,
    const nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    slong N = mpoly_words_per_exp_sp(F->bits, smctx->minfo);
    slong i, j, Fi;
    slong lastdeg = -WORD(1);
    slong off0, shift0, off1, shift1;

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, smctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, smctx->minfo);

    Fi = 0;
    for (i = A->length - 1; i >= 0; i--)
    {
        n_fq_poly_struct * Ai = A->coeffs + i;
        for (j = Ai->length - 1; j >= 0; j--)
        {
            if (_n_fq_is_zero(Ai->coeffs + lgd*j, lgd))
                continue;

            nmod_mpolyn_fit_length(F, Fi + 1, smctx);
            mpoly_monomial_zero(F->exps + N*Fi, N);
            (F->exps + N*Fi)[off0] += (i << shift0);
            (F->exps + N*Fi)[off1] += (j << shift1);
            n_fq_get_n_poly(F->coeffs + Fi, Ai->coeffs + lgd*j, lgctx->fqctx);
            lastdeg = FLINT_MAX(lastdeg, n_poly_degree(F->coeffs + Fi));

            Fi++;
        }
    }

    F->length = Fi;

    *lastdeg_ = lastdeg;
}

int nmod_mpolyn_interp_crt_lg_bpoly(
    slong * lastdeg,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    n_fq_poly_t modulus,
    const nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    int changed = 0;
    slong N = mpoly_words_per_exp(T->bits, smctx->minfo);
    slong off0, shift0, off1, shift1;
    n_fq_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    n_poly_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    n_poly_struct * Tcoeffs = T->coeffs;
    mp_limb_t * u = FLINT_ARRAY_ALLOC(3*lgd, mp_limb_t);
    mp_limb_t * v = u + lgd;
    mp_limb_t * inv_m_eval = v + lgd;
    ulong Fexpi, mask;

    mask = (-UWORD(1)) >> (FLINT_BITS - F->bits);
    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, smctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, smctx->minfo);

    _n_fq_set_n_poly(u, modulus->coeffs, modulus->length, lgctx->fqctx);

    n_fq_inv(inv_m_eval, u, lgctx->fqctx);

    FLINT_ASSERT(nmod_mpolyn_is_canonical(F, smctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, lgctx->fqctx));

    FLINT_ASSERT(T->bits == F->bits);

    *lastdeg = -1;

    Ti = Fi = 0;
    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_fq_poly_degree(Acoeffs + Ai);

    while (Fi < Flen || Ai >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Flen - Fi, Ai);
            nmod_mpolyn_fit_length(T, Ti + extra + 1, smctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        if (Fi < Flen)
        {
            Fexpi = pack_exp2(((Fexps + N*Fi)[off0]>>shift0)&mask,
                              ((Fexps + N*Fi)[off1]>>shift1)&mask);
        }
        else
        {
            Fexpi = 0;
        }

        if (Fi < Flen && Ai >= 0 && Fexpi == pack_exp2(Ai, ai))
        {
            /* F term ok, A term ok */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            _n_fq_set_n_poly(u, Fcoeffs[Fi].coeffs, Fcoeffs[Fi].length, lgctx->fqctx);
            n_fq_sub(v, Acoeffs[Ai].coeffs + lgd*ai, u, lgctx->fqctx);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                _n_poly_mul_n_fq(Tcoeffs + Ti, modulus, u, lgctx->fqctx);
                n_poly_mod_add(Tcoeffs + Ti, Tcoeffs + Ti, Fcoeffs + Fi, smctx->mod);
            }
            else
            {
                n_poly_set(Tcoeffs + Ti, Fcoeffs + Fi);
            }

            Fi++;
            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + lgd*ai, lgd));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Fi < Flen && (Ai < 0 || Fexpi > pack_exp2(Ai, ai)))
        {
            /* F term ok, Aterm missing */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            _n_fq_set_n_poly(v, Fcoeffs[Fi].coeffs, Fcoeffs[Fi].length, lgctx->fqctx);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                _n_poly_mul_n_fq(Tcoeffs + Ti, modulus, u, lgctx->fqctx);
                n_poly_mod_sub(Tcoeffs + Ti, Fcoeffs + Fi, Tcoeffs + Ti, smctx->mod);
            }
            else
            {
                n_poly_set(Tcoeffs + Ti, Fcoeffs + Fi);
            }

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai >= 0 && (Fi >= Flen || Fexpi < pack_exp2(Ai, ai)));

            /* F term missing, A term ok */
            mpoly_monomial_zero(Texps + N*Ti, N);
            (Texps + N*Ti)[off0] += (Ai << shift0);
            (Texps + N*Ti)[off1] += (ai << shift1);

            changed = 1;
            n_fq_mul(u, Acoeffs[Ai].coeffs + lgd*ai, inv_m_eval, lgctx->fqctx);
            _n_poly_mul_n_fq(Tcoeffs + Ti, modulus, u, lgctx->fqctx);

            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + lgd*ai, lgd));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }

        FLINT_ASSERT(!n_poly_is_zero(Tcoeffs + Ti));
        *lastdeg = FLINT_MAX(*lastdeg, n_poly_degree(Tcoeffs + Ti));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        nmod_mpolyn_swap(T, F);

    FLINT_ASSERT(nmod_mpolyn_is_canonical(F, smctx));

    flint_free(u);

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
    n_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    n_poly_struct * Ecoeff;
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
        n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(v), Acoeff + Ai,
          evil_const_cast_nmod_poly_to_n_poly(ectx->fqctx->modulus), ctx->mod);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (fq_nmod_is_zero(v, ectx->fqctx))
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            n_fq_poly_set_coeff_fq_nmod(Ecoeff + Ei - 1, k, v, ectx->fqctx);
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
            n_poly_zero(Ecoeff + Ei);
            n_fq_poly_set_coeff_fq_nmod(Ecoeff + Ei, k, v, ectx->fqctx);
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
    slong d = fq_nmod_ctx_degree(ectx->fqctx);
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;
    n_fq_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;
    n_poly_struct * Acoeff;
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
            if (!_n_fq_is_zero((Bcoeff + Bi)->coeffs + d*vi, d))
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                n_fq_get_n_poly(Acoeff + Ai, (Bcoeff + Bi)->coeffs + d*vi, ectx->fqctx);
                lastdeg = FLINT_MAX(lastdeg, n_poly_degree(Acoeff + Ai));
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
    n_poly_t modulus,
    slong var,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_mpolyn_t A,
    const fq_nmod_mpoly_ctx_t ectx)
{
    slong d = fq_nmod_ctx_degree(ectx->fqctx);
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;
    slong lastdeg = -WORD(1);
    slong vi;
    fq_nmod_t u, v;
    n_poly_t w;
    n_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;
    n_fq_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;
    n_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;
    fq_nmod_t inv_m_eval;
    fq_nmod_t at;

    fq_nmod_init(inv_m_eval, ectx->fqctx);
    n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(inv_m_eval), modulus,
          evil_const_cast_nmod_poly_to_n_poly(ectx->fqctx->modulus), ctx->mod);
    fq_nmod_inv(inv_m_eval, inv_m_eval, ectx->fqctx);

    fq_nmod_init(at, ectx->fqctx);
    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    n_poly_init(w);

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
        vi = n_poly_degree(A->coeffs + Ai);
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
            n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(u), Fcoeff + Fi,
               evil_const_cast_nmod_poly_to_n_poly(ectx->fqctx->modulus), ctx->mod);
            n_fq_get_fq_nmod(at, Acoeff[Ai].coeffs + d*vi, ectx->fqctx);
            fq_nmod_sub(v, at, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                n_poly_mod_mul(w, modulus,
                             evil_const_cast_nmod_poly_to_n_poly(u), ctx->mod);
                n_poly_mod_add(Tcoeff + Ti, Fcoeff + Fi, w, ctx->mod);
            }
            else
            {
                n_poly_set(Tcoeff + Ti, Fcoeff + Fi);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && _n_fq_is_zero(Acoeff[Ai].coeffs + d*vi, d));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = n_poly_degree(A->coeffs + Ai);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(v), Fcoeff + Fi,
                evil_const_cast_nmod_poly_to_n_poly(ectx->fqctx->modulus), ctx->mod);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                n_poly_mod_mul(w, evil_const_cast_nmod_poly_to_n_poly(u),
                                                            modulus, ctx->mod);
                n_poly_mod_sub(Tcoeff + Ti, Fcoeff + Fi, w, ctx->mod);
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
            n_fq_get_fq_nmod(at, Acoeff[Ai].coeffs + d*vi, ectx->fqctx);
            fq_nmod_mul(u, at, inv_m_eval, ectx->fqctx);
            n_poly_mod_mul(Tcoeff + Ti, evil_const_cast_nmod_poly_to_n_poly(u),
                                                            modulus, ctx->mod);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && _n_fq_is_zero(Acoeff[Ai].coeffs + d*vi, d));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = n_poly_degree(A->coeffs + Ai);
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

    fq_nmod_clear(inv_m_eval, ectx->fqctx);
    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_clear(at, ectx->fqctx);
    n_poly_clear(w);

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
    slong d = fq_nmod_ctx_degree(ffctx->fqctx);
    slong N = N = mpoly_words_per_exp(B->bits, ctx->minfo);
    slong i, Alen;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ffctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ffctx->minfo->nvars == ctx->minfo->nvars);

    Alen = 0;
    fq_nmod_mpoly_fit_length(A, B->length, ffctx);
    for (i = 0; i < B->length; i++)
    {
        mpoly_monomial_set(A->exps + N*Alen, B->exps + N*i, N);

        if (B->coeffs[i].length > d)
        {
            _nmod_poly_rem(A->coeffs + d*Alen, B->coeffs[i].coeffs,
                   B->coeffs[i].length, ffctx->fqctx->modulus->coeffs, d + 1,
                                                fq_nmod_ctx_mod(ffctx->fqctx));
        }
        else
        {
            _n_fq_set_n_poly(A->coeffs + d*Alen, B->coeffs[i].coeffs,
                                           B->coeffs[i].length, ffctx->fqctx);
        }

        Alen += !_n_fq_is_zero(A->coeffs + d*Alen, d);
    }
    A->length = Alen;
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
    slong d = fq_nmod_ctx_degree(ctxp->fqctx);
    slong i, N;

    FLINT_ASSERT(Ap->bits == A->bits);
    nmod_mpolyn_fit_length(A, Ap->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    for (i = 0; i < Ap->length; i++)
    {
        mpoly_monomial_set(A->exps + N*i, Ap->exps + N*i, N);
        n_fq_get_n_poly(A->coeffs + i, Ap->coeffs + d*i, ctxp->fqctx);
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
static int nmod_mpolyn_interp_crt_lg_mpoly(
    slong * lastdeg,
    nmod_mpolyn_t F,
    nmod_mpolyn_t T,
    n_poly_t m,
    const nmod_mpoly_ctx_t ctx,
    fq_nmod_mpoly_t A,
    fq_nmod_t inv_m_eval,
    const fq_nmod_mpoly_ctx_t ffctx)
{
    slong d = fq_nmod_ctx_degree(ffctx->fqctx);
    int changed = 0;
    slong i, j, k;
    slong N;
    fq_nmod_t u, v;
    n_poly_t w;
    flint_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    mp_limb_t * Acoeff = A->coeffs;
    n_poly_struct * Fcoeff = F->coeffs;
    n_poly_struct * Tcoeff;
    fq_nmod_t at;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fq_nmod_init(u, ffctx->fqctx);
    fq_nmod_init(v, ffctx->fqctx);
    n_poly_init(w);
    fq_nmod_init(at, ffctx->fqctx);

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
            FLINT_ASSERT(n_poly_degree(Fcoeff + i) < n_poly_degree(m));

            /* F term ok, A term missing */
            n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(v), Fcoeff + i,
               evil_const_cast_nmod_poly_to_n_poly(ffctx->fqctx->modulus), ctx->mod);
            if (!fq_nmod_is_zero(v, ffctx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ffctx->fqctx);
                n_poly_mod_mul(w, evil_const_cast_nmod_poly_to_n_poly(u),
                                                                  m, ctx->mod);
                n_poly_mod_sub(Tcoeff + k, Fcoeff + i, w, ctx->mod);
            } else {
                n_poly_set(Tcoeff + k, Fcoeff + i);
            }
            FLINT_ASSERT(!n_poly_is_zero(Tcoeff + k));
            *lastdeg = FLINT_MAX(*lastdeg, n_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (!_n_fq_is_zero(Acoeff + d*j, d))
            {
                changed = 1;
                n_fq_get_fq_nmod(at, Acoeff + d*j, ffctx->fqctx);
                fq_nmod_mul(u, at, inv_m_eval, ffctx->fqctx);
                n_poly_mod_mul(Tcoeff + k, m,
                             evil_const_cast_nmod_poly_to_n_poly(u), ctx->mod);
                *lastdeg = FLINT_MAX(*lastdeg, n_poly_degree(Tcoeff + k));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!n_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(n_poly_degree(Fcoeff + i) < n_poly_degree(m));

            /* F term ok, A term ok */
            n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(u), Fcoeff + i,
                evil_const_cast_nmod_poly_to_n_poly(ffctx->fqctx->modulus), ctx->mod);
            n_fq_get_fq_nmod(at, Acoeff + d*j, ffctx->fqctx);
            fq_nmod_sub(v, at, u, ffctx->fqctx);
            if (!fq_nmod_is_zero(v, ffctx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ffctx->fqctx);
                n_poly_mod_mul(w, m, evil_const_cast_nmod_poly_to_n_poly(u), ctx->mod);
                n_poly_mod_add(Tcoeff + k, Fcoeff + i, w, ctx->mod);
            } else {
                n_poly_set(Tcoeff + k, Fcoeff + i);                
            }
            FLINT_ASSERT(!n_poly_is_zero(Tcoeff + k));
            *lastdeg = FLINT_MAX(*lastdeg, n_poly_degree(Tcoeff + k));
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

    nmod_mpolyn_set_length(T, k, ctx);

    if (changed)
    {
        nmod_mpolyn_swap(T, F);
    }

    fq_nmod_clear(u, ffctx->fqctx);
    fq_nmod_clear(v, ffctx->fqctx);
    n_poly_clear(w);
    fq_nmod_clear(at, ffctx->fqctx);

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
    n_poly_t m,
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
    fq_nmod_mpoly_fit_length_reset_bits(zero, 0, A->bits, ffctx);

    fq_nmod_init(inv_m_eval, ffctx->fqctx);
    n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(inv_m_eval), m,
          evil_const_cast_nmod_poly_to_n_poly(ffctx->fqctx->modulus), ctx->mod);
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

int nmod_mpolyn_interp_mcrt_lg_mpoly(
    slong * lastdeg_,
    nmod_mpolyn_t H,
    const nmod_mpoly_ctx_t smctx,
    const n_poly_t m,
    const mp_limb_t * inv_m_eval,
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx)
{
    slong lastdeg = *lastdeg_;
    slong i, lgd = fq_nmod_ctx_degree(lgctx->fqctx);
#if FLINT_WANT_ASSERT
    slong N = mpoly_words_per_exp(A->bits, smctx->minfo);
#endif
    int changed = 0;
    mp_limb_t * u = FLINT_ARRAY_ALLOC(lgd, mp_limb_t);
    n_poly_t w;

    n_poly_init(w);

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);

    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        _n_fq_set_n_poly(u, H->coeffs[i].coeffs, H->coeffs[i].length, lgctx->fqctx);
        n_fq_sub(u, A->coeffs + lgd*i, u, lgctx->fqctx);
        if (!_n_fq_is_zero(u, lgd))
        {
            changed = 1;
            n_fq_mul(u, u, inv_m_eval, lgctx->fqctx);
            _n_poly_mul_n_fq(w, m, u, lgctx->fqctx);
            n_poly_mod_add(H->coeffs + i, H->coeffs + i, w, smctx->mod);
        }

        lastdeg = FLINT_MAX(lastdeg, n_poly_degree(H->coeffs + i));

        FLINT_ASSERT(n_poly_degree(H->coeffs + i) < n_poly_degree(m) +
                                      nmod_poly_degree(lgctx->fqctx->modulus));
    }

    *lastdeg_ = lastdeg;

    flint_free(u);
    n_poly_clear(w);
    
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
    n_poly_t m,
    fq_nmod_t inv_m_eval,
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctxp)
{
    slong d = fq_nmod_ctx_degree(ctxp->fqctx);
    slong i;
#if FLINT_WANT_ASSERT
    slong N;
#endif
    int changed = 0;
    fq_nmod_t u, v, at;
    n_poly_t w;

    fq_nmod_init(u, ctxp->fqctx);
    fq_nmod_init(v, ctxp->fqctx);
    fq_nmod_init(at, ctxp->fqctx);
    n_poly_init(w);

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);
#if FLINT_WANT_ASSERT
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
#endif
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(u), H->coeffs + i,
           evil_const_cast_nmod_poly_to_n_poly(ctxp->fqctx->modulus), ctx->mod);
        n_fq_get_fq_nmod(at, A->coeffs + d*i, ctxp->fqctx);
        fq_nmod_sub(v, at, u, ctxp->fqctx);
        if (!fq_nmod_is_zero(v, ctxp->fqctx))
        {
            changed = 1;
            fq_nmod_mul(u, v, inv_m_eval, ctxp->fqctx);
            n_poly_mod_mul(w, evil_const_cast_nmod_poly_to_n_poly(u), m, ctx->mod);
            n_poly_mod_add(H->coeffs + i, H->coeffs + i, w, ctx->mod);
        }

        *lastdeg = FLINT_MAX(*lastdeg, n_poly_degree(H->coeffs + i));

        FLINT_ASSERT(n_poly_degree(H->coeffs + i) < n_poly_degree(m) +
                                       nmod_poly_degree(ctxp->fqctx->modulus));
    }

    fq_nmod_clear(u, ctxp->fqctx);
    fq_nmod_clear(v, ctxp->fqctx);
    fq_nmod_clear(at, ctxp->fqctx);
    n_poly_clear(w);

    return changed;
}

int nmod_mpolyun_interp_mcrt_lg_mpolyu(
    slong * lastdeg,
    nmod_mpolyun_t H,
    const nmod_mpoly_ctx_t ctx,
    n_poly_t m,
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ctxp)
{
    slong i;
    int changed = 0;
    fq_nmod_t inv_m_eval;

    *lastdeg = -WORD(1);

    fq_nmod_init(inv_m_eval, ctxp->fqctx);
    n_poly_mod_rem(evil_cast_nmod_poly_to_n_poly(inv_m_eval), m,
          evil_const_cast_nmod_poly_to_n_poly(ctxp->fqctx->modulus), ctx->mod);
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
    n_poly_struct * Acoeff;
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
        n_fq_poly_evaluate_fq_nmod(v, Acoeff + Ai, alpha, ctx->fqctx);
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
    n_poly_struct * Acoeff;
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

            n_fq_poly_set_fq_nmod(Acoeff + Ai, Bcoeff + Bi, ctx->fqctx);
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
    n_poly_struct * Fcoeffs = F->coeffs;
    ulong * Fexps = F->exps;
    n_poly_struct * Tcoeffs;
    ulong * Texps;
    fq_nmod_poly_t tp;
    slong N, off, shift;
    n_poly_t tpt;

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    Fi = 0;
    Ai = fq_nmod_poly_degree(A, ctx->fqctx);

    fq_nmod_init(u, ctx->fqctx);
    fq_nmod_init(v, ctx->fqctx);
    fq_nmod_poly_init(tp, ctx->fqctx);
    n_poly_init(tpt);

    fq_nmod_mpolyn_fit_length(T, Flen + Ai + 1, ctx);
    Tcoeffs = T->coeffs;
    Texps = T->exps;
    Ti = 0;

    while (Fi < Flen || Ai >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!n_poly_is_zero(Fcoeffs + Fi));
            FLINT_ASSERT(n_poly_degree(Fcoeffs + Fi) < fq_nmod_poly_degree(modulus, ctx->fqctx));
        }

        if (Ai >= 0)
        {
            FLINT_ASSERT(!fq_nmod_is_zero(Acoeff + Ai, ctx->fqctx));
        }

        mpoly_monomial_zero(Texps + N*Ti, N);

        if (Fi < Flen && Ai >= 0 && ((Fexps + N*Fi)[off]>>shift) == Ai)
        {
            /* F term ok, A term ok */
            n_fq_poly_evaluate_fq_nmod(u, Fcoeffs + Fi, alpha, ctx->fqctx);
            fq_nmod_sub(v, Acoeff + Ai, u, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(tpt, tp, ctx->fqctx);
                n_fq_poly_add(Tcoeffs + Ti, Fcoeffs + Fi, tpt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, ctx->fqctx);
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
            n_fq_poly_evaluate_fq_nmod(v, Fcoeffs + Fi, alpha, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(tpt, tp, ctx->fqctx);
                n_fq_poly_sub(Tcoeffs + Ti, Fcoeffs + Fi, tpt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, ctx->fqctx);
            }

            (Texps + N*Ti)[off] = (Fexps + N*Fi)[off];
            Fi++;
        }
        else if (Ai >= 0 && (Fi >= Flen || ((Fexps + N*Fi)[off]>>shift) < Ai))
        {
            /* F term missing, A term ok */
            changed = 1;
            fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, Acoeff + Ai, ctx->fqctx);
            n_fq_poly_set_fq_nmod_poly(Tcoeffs + Ti, tp, ctx->fqctx);

            (Texps + N*Ti)[off] = Ai << shift;
            do {
                Ai--;
            } while (Ai >= 0 && fq_nmod_is_zero(Acoeff + Ai, ctx->fqctx));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        FLINT_ASSERT(!n_poly_is_zero(Tcoeffs + Ti));
        lastdeg = FLINT_MAX(lastdeg, n_poly_degree(Tcoeffs + Ti));

        Ti++;
    }
    T->length = Ti;

    fq_nmod_clear(u, ctx->fqctx);
    fq_nmod_clear(v, ctx->fqctx);
    fq_nmod_poly_clear(tp, ctx->fqctx);
    n_poly_clear(tpt);

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    *lastdeg_ = lastdeg;
    return changed;
}


void fq_nmod_mpolyn_interp_lift_sm_bpoly(
    fq_nmod_mpolyn_t F,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    slong i, j, Fi;
    slong off0, shift0, off1, shift1;

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

    Fi = 0;
    for (i = A->length - 1; i >= 0; i--)
    {
        n_fq_poly_struct * Ai = A->coeffs + i;
        for (j = Ai->length - 1; j >= 0; j--)
        {
            if (_n_fq_is_zero(Ai->coeffs + d*j, d))
                continue;

            fq_nmod_mpolyn_fit_length(F, Fi + 1, ctx);
            mpoly_monomial_zero(F->exps + N*Fi, N);
            (F->exps + N*Fi)[off0] += (i << shift0);
            (F->exps + N*Fi)[off1] += (j << shift1);
            n_fq_poly_set_n_fq(F->coeffs + Fi, Ai->coeffs + d*j, ctx->fqctx);
            Fi++;
        }
    }

    F->length = Fi;
}


int fq_nmod_mpolyn_interp_crt_sm_bpoly(
    slong * lastdeg,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    const n_fq_bpoly_t A,
    const n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(F->bits, ctx->minfo);
    slong off0, shift0, off1, shift1;
    n_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    n_fq_poly_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    n_fq_poly_struct * Tcoeffs = T->coeffs;
    mp_limb_t * v = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    ulong Fexpi, mask;

    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(F, ctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, ctx->fqctx));

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
            slong extra = FLINT_MAX(Flen - Fi, Ai);
            fq_nmod_mpolyn_fit_length(T, Ti + extra + 1, ctx);
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

            n_fq_poly_eval_pow(v, Fcoeffs + Fi, alphapow, ctx->fqctx);
            n_fq_sub(v, Acoeffs[Ai].coeffs + d*ai, v, ctx->fqctx);
            if (!_n_fq_is_zero(v, d))
            {
                changed = 1;
                n_fq_poly_scalar_addmul_n_fq(Tcoeffs + Ti,
                                         Fcoeffs + Fi, modulus, v, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, ctx->fqctx);
            }

            Fi++;
            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Ai >= 0 && (Fi >= Flen || Fexpi < pack_exp2(Ai, ai)))
        {
            /* F term missing, A term ok */

            mpoly_monomial_zero(Texps + N*Ti, N);
            (Texps + N*Ti)[off0] += (Ai << shift0);
            (Texps + N*Ti)[off1] += (ai << shift1);

            changed = 1;
            n_fq_poly_scalar_mul_n_fq(Tcoeffs + Ti, modulus,
                                        Acoeffs[Ai].coeffs + d*ai, ctx->fqctx);

            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }
        else
        {
            FLINT_ASSERT(Fi < Flen && (Ai < 0 || Fexpi > pack_exp2(Ai, ai)));
            /* F term ok, Aterm missing */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            n_fq_poly_eval_pow(v, Fcoeffs + Fi, alphapow, ctx->fqctx);
            if (!_n_fq_is_zero(v, d))
            {
                changed = 1;
                _n_fq_neg(v, v, d, ctx->fqctx->mod);
                n_fq_poly_scalar_addmul_n_fq(Tcoeffs + Ti, 
                                         Fcoeffs + Fi, modulus, v, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, ctx->fqctx);                
            }

            Fi++;
        }

        FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeffs + Ti));
        *lastdeg = FLINT_MAX(*lastdeg, n_fq_poly_degree(Tcoeffs + Ti));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        fq_nmod_mpolyn_swap(T, F);

    flint_free(v);

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
    n_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    n_poly_struct * Ecoeff;
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
        n_fq_poly_evaluate_fq_nmod(v, Acoeff + Ai, alpha, ctx->fqctx);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (fq_nmod_is_zero(v, ctx->fqctx))
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            n_fq_poly_set_coeff_fq_nmod(Ecoeff + Ei - 1, k, v, ctx->fqctx);
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
            n_poly_zero(Ecoeff + Ei);
            n_fq_poly_set_coeff_fq_nmod(Ecoeff + Ei, k, v, ctx->fqctx);
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
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;

    n_fq_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;

    n_fq_poly_struct * Acoeff;
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
            if (!_n_fq_is_zero(Bcoeff[Bi].coeffs + d*vi, d))
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                n_fq_poly_zero(Acoeff + Ai);
                n_fq_poly_set_coeff_n_fq(Acoeff + Ai, 0, Bcoeff[Bi].coeffs + d*vi, ctx->fqctx);
                Ai++;
            }
        }
    }
    A->length = Ai;
}


/* F = F + modulus*(A - F(alpha)) with matching monomials */
int fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
    slong * lastdeg_,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpoly_t A,  /* could have zero coeffs */
    const n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
#if FLINT_WANT_ASSERT
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
#endif
    slong lastdeg = *lastdeg_;
    int changed = 0;
    mp_limb_t * v = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    slong i, Alen = A->length;
    mp_limb_t * Acoeffs = A->coeffs;
    n_fq_poly_struct * Fcoeffs = F->coeffs;

    FLINT_ASSERT(F->bits == A->bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    FLINT_ASSERT(F->length == Alen);
    for (i = 0; i < Alen; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(F->exps + N*i, A->exps + N*i, N));
        FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeffs + i));
        FLINT_ASSERT(n_fq_poly_degree(Fcoeffs + i) < n_fq_poly_degree(modulus));

        n_fq_poly_eval_pow(v, Fcoeffs + i, alphapow, ctx->fqctx);
        n_fq_sub(v, Acoeffs + d*i, v, ctx->fqctx);
        if (!_n_fq_is_zero(v, d))
        {
            changed = 1;
            n_fq_poly_scalar_addmul_n_fq(Fcoeffs + i, Fcoeffs + i,
                                                       modulus, v, ctx->fqctx);
        }

        FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeffs + i));
        lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(Fcoeffs + i));
    }

    flint_free(v);

    *lastdeg_ = lastdeg;
    return changed;
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
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong lastdeg = -WORD(1);
    slong offset, shift;
    slong vi;
    fq_nmod_t v, at;
    fq_nmod_poly_t tp;
    n_fq_poly_t tpt;

    n_fq_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    n_fq_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    n_fq_poly_struct * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;

    fq_nmod_poly_init(tp, ctx->fqctx);
    fq_nmod_init(v, ctx->fqctx);
    fq_nmod_init(at, ctx->fqctx);
    n_fq_poly_init(tpt);

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
        vi = n_fq_poly_degree(A->coeffs + Ai);
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
            FLINT_ASSERT(!_n_fq_is_zero(Acoeff[Ai].coeffs + d*vi, d));
        }

        if (Fi < Flen && Ai < Alen && mpoly_monomial_equal_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift))
        {
            /* F term ok, A term ok */
            n_fq_poly_evaluate_fq_nmod(v, Fcoeff + Fi, alpha, ctx->fqctx);
            n_fq_get_fq_nmod(at, Acoeff[Ai].coeffs + d*vi, ctx->fqctx);
            fq_nmod_sub(v, at, v, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(tpt, tp, ctx->fqctx);
                n_fq_poly_add(Tcoeff + Ti, Fcoeff + Fi, tpt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && _n_fq_is_zero(Acoeff[Ai].coeffs + d*vi, d));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = n_fq_poly_degree(A->coeffs + Ai);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            n_fq_poly_evaluate_fq_nmod(v, Fcoeff + Fi, alpha, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(tpt, tp, ctx->fqctx);
                n_fq_poly_sub(Tcoeff + Ti, Fcoeff + Fi, tpt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            changed = 1;
            n_fq_get_fq_nmod(at, Acoeff[Ai].coeffs + d*vi, ctx->fqctx);
            fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, at, ctx->fqctx);
            n_fq_poly_set_fq_nmod_poly(Tcoeff + Ti, tp, ctx->fqctx);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && _n_fq_is_zero(Acoeff[Ai].coeffs + d*vi, d));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = n_fq_poly_degree(A->coeffs + Ai);
                }
            }
        }

        FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    fq_nmod_poly_clear(tp, ctx->fqctx);
    fq_nmod_clear(v, ctx->fqctx);
    fq_nmod_clear(at, ctx->fqctx);
    n_fq_poly_clear(tpt);

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
    n_fq_poly_struct * Acoeff;
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
        bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(v, Acoeff + Ai, emb);
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
    n_fq_poly_struct * Acoeff;
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
            bad_fq_nmod_embed_fq_nmod_lg_to_n_fq_sm(Acoeff + Ai, Bcoeff + Bexp, emb);
            FLINT_ASSERT(!n_fq_poly_is_zero(Acoeff + Ai));
            lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(Acoeff + Ai));
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
    n_fq_poly_t wt;
    slong Fi, Ti, Aexp;
    fq_nmod_struct * Acoeff = A->coeffs;
    slong Flen = F->length;
    n_fq_poly_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    n_fq_poly_struct * Tcoeff;
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
    n_fq_poly_init(wt);

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
            FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeff + Fi));
            FLINT_ASSERT(n_fq_poly_degree(Fcoeff + Fi) < fq_nmod_poly_degree(modulus, ctx->fqctx));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(!fq_nmod_is_zero(Acoeff + Aexp, ectx->fqctx));
        }

        mpoly_monomial_zero(Texp + N*Ti, N);

        if (Fi < Flen && Aexp >= 0 && ((Fexp + N*Fi)[off]>>shift) == Aexp)
        {
            /* F term ok, A term ok */
            bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(u, Fcoeff + Fi, emb);
            fq_nmod_sub(v, Acoeff + Aexp, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, modulus, u_sm, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(wt, w, ctx->fqctx);
                n_fq_poly_add(Tcoeff + Ti, Fcoeff + Fi, wt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
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
            bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(v, Fcoeff + Fi, emb);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, u_sm, modulus, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(wt, w, ctx->fqctx);
                n_fq_poly_add(Tcoeff + Ti, Fcoeff + Fi, wt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
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
            fq_nmod_poly_mul(w, modulus, u_sm, ctx->fqctx);
            n_fq_poly_set_fq_nmod_poly(Tcoeff + Ti, w, ctx->fqctx);

            (Texp + N*Ti)[off] = Aexp << shift;
            do {
                Aexp--;
            } while (Aexp >= 0 && fq_nmod_is_zero(Acoeff + Aexp, ectx->fqctx));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(Tcoeff + Ti));
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
    n_fq_poly_clear(wt);

    fq_nmod_clear(inv_m_eval, ectx->fqctx);

    *lastdeg_ = lastdeg;
    return changed;
}


void fq_nmod_mpolyn_interp_lift_lg_bpoly(
    slong * lastdeg_,
    fq_nmod_mpolyn_t F,
    const fq_nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx,
    const bad_fq_nmod_embed_t emb)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    slong N = mpoly_words_per_exp_sp(F->bits, smctx->minfo);
    slong i, j, Fi;
    slong lastdeg = -WORD(1);
    slong off0, shift0, off1, shift1;

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, smctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, smctx->minfo);

    Fi = 0;
    for (i = A->length - 1; i >= 0; i--)
    {
        n_fq_poly_struct * Ai = A->coeffs + i;
        for (j = Ai->length - 1; j >= 0; j--)
        {
            if (_n_fq_is_zero(Ai->coeffs + lgd*j, lgd))
                continue;

            fq_nmod_mpolyn_fit_length(F, Fi + 1, smctx);
            mpoly_monomial_zero(F->exps + N*Fi, N);
            (F->exps + N*Fi)[off0] += (i << shift0);
            (F->exps + N*Fi)[off1] += (j << shift1);
            bad_n_fq_embed_lg_to_sm(F->coeffs + Fi, Ai->coeffs + lgd*j, emb);
            lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(F->coeffs + Fi));

            Fi++;
        }
    }

    F->length = Fi;

    *lastdeg_ = lastdeg;
}

int fq_nmod_mpolyn_interp_crt_lg_bpoly(
    slong * lastdeg,
    fq_nmod_mpolyn_t F,
    fq_nmod_mpolyn_t T,
    n_fq_poly_t modulus,
    const fq_nmod_mpoly_ctx_t smctx,
    n_fq_bpoly_t A,
    const fq_nmod_mpoly_ctx_t lgctx,
    const bad_fq_nmod_embed_t emb)
{
    slong lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    int changed = 0;
    slong N = mpoly_words_per_exp(T->bits, smctx->minfo);
    slong off0, shift0, off1, shift1;
    n_fq_poly_struct * Acoeffs = A->coeffs;
    slong Fi, Ti, Ai, ai;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    n_fq_poly_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    n_fq_poly_struct * Tcoeffs = T->coeffs;
    mp_limb_t * u = FLINT_ARRAY_ALLOC(3*lgd, mp_limb_t);
    mp_limb_t * v = u + lgd;
    mp_limb_t * inv_m_eval = v + lgd;
    n_fq_poly_t u_sm;
    ulong Fexpi, mask;

    mask = (-UWORD(1)) >> (FLINT_BITS - F->bits);
    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, smctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, smctx->minfo);

    n_fq_poly_init(u_sm);

    bad_n_fq_embed_sm_to_lg(u, modulus, emb);
    n_fq_inv(inv_m_eval, u, lgctx->fqctx);

    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(F, smctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, lgctx->fqctx));

    FLINT_ASSERT(T->bits == F->bits);

    *lastdeg = -1;

    Ti = Fi = 0;
    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_fq_poly_degree(Acoeffs + Ai);

    while (Fi < Flen || Ai >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Flen - Fi, Ai);
            fq_nmod_mpolyn_fit_length(T, Ti + extra + 1, smctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        if (Fi < Flen)
        {
            Fexpi = pack_exp2(((Fexps + N*Fi)[off0]>>shift0)&mask,
                              ((Fexps + N*Fi)[off1]>>shift1)&mask);
        }
        else
        {
            Fexpi = 0;
        }

        if (Fi < Flen && Ai >= 0 && Fexpi == pack_exp2(Ai, ai))
        {
            /* F term ok, A term ok */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            bad_n_fq_embed_sm_to_lg(u, Fcoeffs + Fi, emb);
            n_fq_sub(v, Acoeffs[Ai].coeffs + lgd*ai, u, lgctx->fqctx);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
                n_fq_poly_mul(Tcoeffs + Ti, modulus, u_sm, smctx->fqctx);
                n_fq_poly_add(Tcoeffs + Ti, Tcoeffs + Ti, Fcoeffs + Fi, smctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, smctx->fqctx);
            }

            Fi++;
            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + lgd*ai, lgd));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Fi < Flen && (Ai < 0 || Fexpi > pack_exp2(Ai, ai)))
        {
            /* F term ok, Aterm missing */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            bad_n_fq_embed_sm_to_lg(v, Fcoeffs + Fi, emb);
            if (!_n_fq_is_zero(v, lgd))
            {
                changed = 1;
                n_fq_mul(u, v, inv_m_eval, lgctx->fqctx);
                bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
                n_fq_poly_mul(Tcoeffs + Ti, modulus, u_sm, smctx->fqctx);
                n_fq_poly_sub(Tcoeffs + Ti, Fcoeffs + Fi, Tcoeffs + Ti, smctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeffs + Ti, Fcoeffs + Fi, smctx->fqctx);
            }

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai >= 0 && (Fi >= Flen || Fexpi < pack_exp2(Ai, ai)));

            /* F term missing, A term ok */
            mpoly_monomial_zero(Texps + N*Ti, N);
            (Texps + N*Ti)[off0] += (Ai << shift0);
            (Texps + N*Ti)[off1] += (ai << shift1);

            changed = 1;
            n_fq_mul(u, Acoeffs[Ai].coeffs + lgd*ai, inv_m_eval, lgctx->fqctx);
            bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
            n_fq_poly_mul(Tcoeffs + Ti, modulus, u_sm, smctx->fqctx);

            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + lgd*ai, lgd));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }

        FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeffs + Ti));
        *lastdeg = FLINT_MAX(*lastdeg, n_fq_poly_degree(Tcoeffs + Ti));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        fq_nmod_mpolyn_swap(T, F);

    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(F, smctx));

    n_fq_poly_clear(u_sm);
    flint_free(u);

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
    n_fq_poly_struct * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    n_fq_poly_struct * Ecoeff;
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
        bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(v, Acoeff + Ai, emb);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (fq_nmod_is_zero(v, ectx->fqctx))
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)))
        {
            /* append to previous */
            n_fq_poly_set_coeff_fq_nmod(Ecoeff + Ei - 1, k, v, ectx->fqctx);
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
            n_fq_poly_zero(Ecoeff + Ei);
            n_fq_poly_set_coeff_fq_nmod(Ecoeff + Ei, k, v, ectx->fqctx);
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
    slong lgd = fq_nmod_ctx_degree(ectx->fqctx);
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;
    n_fq_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;
    n_fq_poly_struct * Acoeff;
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
        if (Ai + Bcoeff[Bi].length >= A->alloc)
        {
            fq_nmod_mpolyn_fit_length(A, Ai + Bcoeff[Bi].length, ctx);
            Acoeff = A->coeffs;
            Aexp = A->exps;
        }
        for (vi = Bcoeff[Bi].length - 1; vi >= 0; vi--)
        {
            if (_n_fq_is_zero(Bcoeff[Bi].coeffs + lgd*vi, lgd))
                continue;

            mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
            bad_n_fq_embed_lg_to_sm(Acoeff + Ai, Bcoeff[Bi].coeffs + lgd*vi, emb);
            FLINT_ASSERT(!n_fq_poly_is_zero(Acoeff + Ai));
            lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(Acoeff + Ai));
            Ai++;
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
    slong lgd = fq_nmod_ctx_degree(ectx->fqctx);
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift;
    slong lastdeg = -WORD(1);
    slong vi;
    fq_nmod_t u, v;
    fq_nmod_poly_t u_sm, w;
    fq_nmod_t inv_m_eval;
    fq_nmod_t at;
    n_fq_poly_t wt;

    n_fq_poly_struct * Tcoeff;
    ulong * Texp;
    slong Ti;

    n_fq_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;

    n_fq_poly_struct * Fcoeff = F->coeffs;
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
    fq_nmod_init(at, ectx->fqctx);
    n_fq_poly_init(wt);

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
        vi = n_fq_poly_degree(A->coeffs + Ai);
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
            bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(u, Fcoeff + Fi, emb);
            n_fq_get_fq_nmod(at, Acoeff[Ai].coeffs + lgd*vi, ectx->fqctx);
            fq_nmod_sub(v, at, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, modulus, u_sm, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(wt, w, ctx->fqctx);
                n_fq_poly_add(Tcoeff + Ti, Fcoeff + Fi, wt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && _n_fq_is_zero(Acoeff[Ai].coeffs + lgd*vi, lgd));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = n_fq_poly_degree(A->coeffs + Ai);
                }
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(v, Fcoeff + Fi, emb);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, u_sm, modulus, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(wt, w, ctx->fqctx);
                n_fq_poly_sub(Tcoeff + Ti, Fcoeff + Fi, wt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->fqctx);
            }
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen || mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)));

            /* F term missing, A term ok */
            changed = 1;
            n_fq_get_fq_nmod(at, Acoeff[Ai].coeffs + lgd*vi, ectx->fqctx);
            fq_nmod_mul(u, at, inv_m_eval, ectx->fqctx);
            bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
            fq_nmod_poly_mul(w, modulus, u_sm, ctx->fqctx);
            n_fq_poly_set_fq_nmod_poly(Tcoeff + Ti, w, ctx->fqctx);
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && _n_fq_is_zero(Acoeff[Ai].coeffs + lgd*vi, lgd));
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    vi = n_fq_poly_degree(A->coeffs + Ai);
                }
            }
        }

        FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + Ti));
        lastdeg = FLINT_MAX(lastdeg, n_fq_poly_degree(Tcoeff + Ti));
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
    fq_nmod_clear(at, ectx->fqctx);
    n_fq_poly_clear(wt);

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
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, k, N;
    fq_nmod_t bt;
    FLINT_ASSERT(B->bits == A->bits);

    fq_nmod_init(bt, ctx->fqctx);

    fq_nmod_mpoly_fit_length(B, A->length, ctx);
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    k = 0;
    for (i = 0; i < A->length; i++)
    {
        mpoly_monomial_set(B->exps + N*k, A->exps + N*i, N);
        n_fq_poly_evaluate_fq_nmod(bt, A->coeffs + i, alpha, ctx->fqctx);
        n_fq_set_fq_nmod(B->coeffs + d*k, bt, ctx->fqctx);
        k += !_n_fq_is_zero(B->coeffs + d*k, d);
    }
    B->length = k;

    fq_nmod_clear(bt, ctx->fqctx);
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
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i;
    slong N;
    n_fq_poly_struct * Acoeff;
    mp_limb_t * Bcoeff;
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
        n_fq_poly_set_n_fq(Acoeff + i, Bcoeff + d*i, ctx->fqctx);
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
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    int changed = 0;
    slong i, j, k;
    fq_nmod_t v;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    mp_limb_t * Acoeff = A->coeffs;
    n_fq_poly_struct * Fcoeff = F->coeffs;
    n_fq_poly_struct * Tcoeff;
    fq_nmod_poly_t tp;
    n_fq_poly_t tpt;
    fq_nmod_t at;

    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fq_nmod_init(v, ctx->fqctx);
    fq_nmod_poly_init(tp, ctx->fqctx);
    n_fq_poly_init(tpt);
    fq_nmod_init(at, ctx->fqctx);

    fq_nmod_mpolyn_fit_length(T, Flen + Alen, ctx);
    Texp = T->exps;
    Tcoeff = T->coeffs;


    i = j = k = 0;
    while (i < Flen || j < Alen)
    {
        if (i < Flen && (j >= Alen
                        || mpoly_monomial_gt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(n_fq_poly_degree(Fcoeff + i) < fq_nmod_poly_degree(modulus, ctx->fqctx));

            /* F term ok, A term missing */
            n_fq_poly_evaluate_fq_nmod(v, Fcoeff + i, alpha, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(tpt, tp, ctx->fqctx);
                n_fq_poly_sub(Tcoeff + k, Fcoeff + i, tpt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);                
            }

            FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + k));
            lastdeg[0] = FLINT_MAX(lastdeg[0], n_fq_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (!_n_fq_is_zero(Acoeff + d*j, d))
            {
                changed = 1;
                n_fq_get_fq_nmod(at, Acoeff + d*j, ctx->fqctx);
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, at, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(Tcoeff + k, tp, ctx->fqctx);

                FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + k));
                lastdeg[0] = FLINT_MAX(lastdeg[0], n_fq_poly_degree(Tcoeff + k));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(n_fq_poly_degree(Fcoeff + i) < fq_nmod_poly_degree(modulus, ctx->fqctx));

            /* F term ok, A term ok */
            n_fq_poly_evaluate_fq_nmod(v, Fcoeff + i, alpha, ctx->fqctx);
            n_fq_get_fq_nmod(at, Acoeff + d*j, ctx->fqctx);
            fq_nmod_sub(v, at, v, ctx->fqctx);
            if (!fq_nmod_is_zero(v, ctx->fqctx))
            {
                changed = 1;
                fq_nmod_poly_scalar_mul_fq_nmod(tp, modulus, v, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(tpt, tp, ctx->fqctx);
                n_fq_poly_add(Tcoeff + k, Fcoeff + i, tpt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);
            }

            FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + k));
            lastdeg[0] = FLINT_MAX(lastdeg[0], n_fq_poly_degree(Tcoeff + k));

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

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    fq_nmod_poly_clear(tp, ctx->fqctx);
    fq_nmod_clear(v, ctx->fqctx);
    n_fq_poly_clear(tpt);
    fq_nmod_clear(at, ctx->fqctx);

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
    fq_nmod_mpoly_fit_length_reset_bits(zero, 0, A->bits, ctx);

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
    slong lgd = fq_nmod_ctx_degree(ectx->fqctx);
    slong N = mpoly_words_per_exp_sp(B->bits, ctx->minfo);
    slong i, k;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ectx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ectx->minfo->nvars == ctx->minfo->nvars);

    k = 0;
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_fit_length(A, k + 1, ectx);
        mpoly_monomial_set(A->exps + N*k, B->exps + N*i, N);
        bad_n_fq_embed_sm_to_lg(A->coeffs + lgd*k, B->coeffs + i, emb);
        k += !_n_fq_is_zero(A->coeffs + lgd*k, lgd);
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
    slong lgd = fq_nmod_ctx_degree(ectx->fqctx);
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong i;

    FLINT_ASSERT(B->bits == A->bits);

    fq_nmod_mpolyn_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
        bad_n_fq_embed_lg_to_sm(A->coeffs + i, B->coeffs + lgd*i, emb);
        FLINT_ASSERT(!n_fq_poly_is_zero(A->coeffs + i));
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
    slong lgd = fq_nmod_ctx_degree(ectx->fqctx);
    int changed = 0;
    slong i, j, k;
    slong N;
    fq_nmod_t u, v;
    fq_nmod_poly_t u_sm, w;
    flint_bitcnt_t bits = A->bits;
    slong Flen = F->length, Alen = A->length;
    ulong * Fexp = F->exps, * Aexp = A->exps;
    ulong * Texp;
    mp_limb_t * Acoeff = A->coeffs;
    n_fq_poly_struct * Fcoeff = F->coeffs;
    n_fq_poly_struct * Tcoeff;
    fq_nmod_t at;
    n_fq_poly_t wt;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(F->bits == bits);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    fq_nmod_poly_init(u_sm, ctx->fqctx);
    fq_nmod_poly_init(w, ctx->fqctx);
    n_fq_poly_init(wt);
    fq_nmod_init(at, ectx->fqctx);

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
            FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(n_fq_poly_degree(Fcoeff + i) < fq_nmod_poly_degree(m, ctx->fqctx));

            /* F term ok, A term missing */
            bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(v, Fcoeff + i, emb);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, u_sm, m, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(wt, w, ctx->fqctx);
                n_fq_poly_sub(Tcoeff + k, Fcoeff + i, wt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);
            }

            FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + k));
            lastdeg[0] = FLINT_MAX(lastdeg[0], n_fq_poly_degree(Tcoeff + k));

            mpoly_monomial_set(Texp + N*k, Fexp + N*i, N);
            k++;
            i++;
        }
        else if (j < Alen && (i >= Flen
                        || mpoly_monomial_lt_nomask(Fexp + N*i, Aexp + N*j, N)))
        {
            /* F term missing, A term ok */
            if (!_n_fq_is_zero(Acoeff + lgd*j, lgd))
            {
                changed = 1;
                n_fq_get_fq_nmod(at, Acoeff + lgd*j, ectx->fqctx);
                fq_nmod_mul(u, at, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, m, u_sm, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(Tcoeff + k, w, ctx->fqctx);

                FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + k));
                lastdeg[0] = FLINT_MAX(lastdeg[0], n_fq_poly_degree(Tcoeff + k));
                mpoly_monomial_set(Texp + N*k, Aexp + N*j, N);
                k++;
            }
            j++;
        }
        else if (i < Flen && j < Alen
                             && mpoly_monomial_equal(Fexp + N*i, Aexp + N*j, N))
        {
            FLINT_ASSERT(!n_fq_poly_is_zero(Fcoeff + i));
            FLINT_ASSERT(n_fq_poly_degree(Fcoeff + i) < fq_nmod_poly_degree(m, ctx->fqctx));

            /* F term ok, A term ok */
            bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(u, Fcoeff + i, emb);
            n_fq_get_fq_nmod(at, Acoeff + lgd*j, ectx->fqctx);
            fq_nmod_sub(v, at, u, ectx->fqctx);
            if (!fq_nmod_is_zero(v, ectx->fqctx))
            {
                changed = 1;
                fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
                bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
                fq_nmod_poly_mul(w, m, u_sm, ctx->fqctx);
                n_fq_poly_set_fq_nmod_poly(wt, w, ctx->fqctx);
                n_fq_poly_add(Tcoeff + k, Fcoeff + i, wt, ctx->fqctx);
            }
            else
            {
                n_fq_poly_set(Tcoeff + k, Fcoeff + i, ctx->fqctx);
            }

            FLINT_ASSERT(!n_fq_poly_is_zero(Tcoeff + k));
            lastdeg[0] = FLINT_MAX(lastdeg[0], n_fq_poly_degree(Tcoeff + k));
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

    if (changed)
    {
        fq_nmod_mpolyn_swap(T, F);
    }

    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_poly_clear(u_sm, ctx->fqctx);
    fq_nmod_poly_clear(w, ctx->fqctx);
    n_fq_poly_clear(wt);
    fq_nmod_clear(at, ectx->fqctx);

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
    fq_nmod_mpoly_fit_length_reset_bits(zero, 0, A->bits, ectx);

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
    slong lgd = fq_nmod_ctx_degree(ectx->fqctx);
    slong i;
#if FLINT_WANT_ASSERT
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
#endif
    int changed = 0;
    fq_nmod_t u, v, at;
    fq_nmod_poly_t w, u_sm;
    n_fq_poly_t wt;

    fq_nmod_init(u, ectx->fqctx);
    fq_nmod_init(v, ectx->fqctx);
    fq_nmod_poly_init(w, ctx->fqctx);
    n_fq_poly_init(wt);
    fq_nmod_poly_init(u_sm, ctx->fqctx);
    fq_nmod_init(at, ectx->fqctx);

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);

    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(u, H->coeffs + i, emb);
        n_fq_get_fq_nmod(at, A->coeffs + lgd*i, ectx->fqctx);
        fq_nmod_sub(v, at, u, ectx->fqctx);
        if (!fq_nmod_is_zero(v, ectx->fqctx))
        {
            changed = 1;
            fq_nmod_mul(u, v, inv_m_eval, ectx->fqctx);
            bad_fq_nmod_embed_lg_to_sm(u_sm, u, emb);
            fq_nmod_poly_mul(w, u_sm, m, ctx->fqctx);
            n_fq_poly_set_fq_nmod_poly(wt, w, ctx->fqctx);
            n_fq_poly_add(H->coeffs + i, H->coeffs + i, wt, ctx->fqctx);
        }

        *lastdeg = FLINT_MAX(*lastdeg, n_fq_poly_degree(H->coeffs + i));

        FLINT_ASSERT(n_fq_poly_degree(H->coeffs + i)
                         <  fq_nmod_poly_degree(m, ctx->fqctx)
                          + fq_nmod_poly_degree(emb->h, ctx->fqctx));
    }

    fq_nmod_clear(u, ectx->fqctx);
    fq_nmod_clear(v, ectx->fqctx);
    fq_nmod_poly_clear(w, ctx->fqctx);
    n_fq_poly_clear(wt);
    fq_nmod_poly_clear(u_sm, ctx->fqctx);
    fq_nmod_clear(at, ectx->fqctx);

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


