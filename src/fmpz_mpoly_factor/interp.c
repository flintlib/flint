/*
    Copyright (C) 2019, 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "nmod_mpoly.h"

/*
    interp_reduce: map from Z to Z/pZ
    interp_lift:   map from Z/pZ to Z
    interp_crt:    update element of Z with a new image in Z/pZ
    interp_mcrt:   same as interp_crt, but monomial match, thus easier
*/

/*****************************************************************************/

/*
    Ap = A mod p
    Ap is in Fp[x_0, ..., x_(n-1)]
    A  is in ZZ[x_0, ..., x_(n-1)], n = ctx->minfo->nvars
*/
void fmpz_mpoly_interp_reduce_p(
    nmod_mpoly_t Ap,
    const nmod_mpoly_ctx_t ctxp,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, k, N;

    FLINT_ASSERT(Ap->bits == A->bits);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    nmod_mpoly_fit_length(Ap, A->length, ctxp);
    k = 0;
    for (i = 0; i < A->length; i++)
    {
        mpoly_monomial_set(Ap->exps + N*k, A->exps + N*i, N);
        Ap->coeffs[k] = fmpz_get_nmod(A->coeffs + i, ctxp->mod);
        k += (Ap->coeffs[k] != 0);
    }
    Ap->length = k;
}

/*
    Convert Ap to A using the symmetric range (-p/2, p/2]
    A  is in ZZ[x_0, ..., x_(n-1)], n = ctx->minfo->nvars
    Ap is in Fp[x_0, ..., x_(n-1)]
*/
void fmpz_mpoly_interp_lift_p(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    nmod_mpoly_t Ap,
    const nmod_mpoly_ctx_t ctxp)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    FLINT_ASSERT(Ap->bits == A->bits);

    fmpz_mpoly_fit_length(A, Ap->length, ctx);
    mpoly_copy_monomials(A->exps, Ap->exps, Ap->length, N);
    _fmpz_vec_set_nmod_vec(A->coeffs, Ap->coeffs, Ap->length, ctxp->mod);
    _fmpz_mpoly_set_length(A, Ap->length, ctx);
}

/*
    Update H so that it does not change mod m, and is now A mod p
    It is asserted that the monomials in H and A match
    H is in ZZ[x_0, ..., x_(n-1)], n = ctx->minfo->nvars
    A is in Fp[x_0, ..., x_(n-1)]
*/
int fmpz_mpoly_interp_mcrt_p(
    flint_bitcnt_t * coeffbits_,
    fmpz_mpoly_t H,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_t m,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctxp)
{
    slong i;
#if FLINT_WANT_ASSERT
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
#endif
    int changed = 0;
    flint_bitcnt_t coeffbits = 0;
    fmpz_t t;

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);

    fmpz_init(t);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        fmpz_CRT_ui(t, H->coeffs + i, m, A->coeffs[i], ctxp->mod.n, 1);
        coeffbits = FLINT_MAX(coeffbits, fmpz_bits(t));
        changed |= !fmpz_equal(t, H->coeffs + i);
        fmpz_swap(t, H->coeffs + i);
    }
    fmpz_clear(t);
    *coeffbits_ = coeffbits;
    return changed;
}

/*****************************************************************************/

/*
    E = A mod p
    E is in Fp[x_0,...,x_(m-2)][x_(m-1)]
    A is in ZZ[x_0,...,x_(m-2), x_(m-1)], m = ctx->minfo->nvars
*/
void fmpz_mpoly_interp_reduce_p_mpolyn(
    nmod_mpolyn_t E,
    const nmod_mpoly_ctx_t pctx,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong offset, shift, k;
    ulong mask;
    mp_limb_t v;
    fmpz * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    slong Alen = A->length;
    slong Ai;
    n_poly_struct * Ecoeff;
    ulong * Eexp;
    slong Ei;
    slong m = ctx->minfo->nvars;

    mpoly_gen_offset_shift_sp(&offset, &shift, m - 1, A->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);

    Ecoeff = E->coeffs;
    Eexp = E->exps;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        v = fmpz_fdiv_ui(Acoeff + Ai, pctx->mod.n);
        k = ((Aexp + N*Ai)[offset] >> shift) & mask;
        if (v == 0)
        {
            continue;
        }

        if (Ei > 0 && mpoly_monomial_equal_extra(Eexp + N*(Ei - 1),
                                        Aexp + N*Ai, N, offset, -(k << shift)))
        {
            n_poly_set_coeff(Ecoeff + Ei - 1, k, v);
        }
        else
        {
            FLINT_ASSERT(Ei == 0 || mpoly_monomial_gt_nomask_extra(
                    Eexp + N*(Ei - 1), Aexp + N*Ai, N, offset, -(k << shift)));

            /* create new */
            if (Ei >= E->alloc)
            {
                nmod_mpolyn_fit_length(E, Ei + 1, pctx);
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
    A = B using symmetric range
    A is in ZZ[x_0, ..., x_(var-2), x_(m-1)]
    B is in Fp[x_0, ..., x_(var-2)][x_(m-1)]
*/
void fmpz_mpoly_interp_lift_p_mpolyn(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t pctx)
{
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    slong offset, shift;
    slong vi;
    n_poly_struct * Bcoeff = B->coeffs;
    ulong * Bexp = B->exps;
    slong Blen = B->length;
    slong Bi;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong var = ctx->minfo->nvars;

    FLINT_ASSERT(var == pctx->minfo->nvars);

    fmpz_mpoly_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, A->bits, ctx->minfo);

    Ai = 0;
    for (Bi = 0; Bi < Blen; Bi++)
    {
        if (Ai + (Bcoeff + Bi)->length >= A->alloc)
        {
            fmpz_mpoly_fit_length(A, Ai + (Bcoeff + Bi)->length, ctx);
            Acoeff = A->coeffs;
            Aexp = A->exps;
        }
        for (vi = (Bcoeff + Bi)->length - 1; vi >= 0; vi--)
        {
            if ((Bcoeff + Bi)->coeffs[vi] != 0)
            {
                mpoly_monomial_set_extra(Aexp + N*Ai, Bexp + N*Bi, N, offset, vi << shift);
                fmpz_set_ui_smod(Acoeff + Ai, (Bcoeff + Bi)->coeffs[vi], pctx->mod.n);
                Ai++;
            }
        }
    }
    A->length = Ai;
}


/*
    F = F + modulus*((A - F(mod p))/(modulus (mod p)))
    no assumptions about matching monomials
    F is in ZZ[x_0, ..., x_(m-1), x_(m-1)]
    A is in Fp[x_0, ..., x_(m-2)][x_(m-1)]
*/

int fmpz_mpoly_interp_crt_p_mpolyn(
    fmpz_mpoly_t F,
    fmpz_mpoly_t T,
    const fmpz_mpoly_ctx_t ctx,
    fmpz_t modulus,
    const nmod_mpolyn_t A,
    const nmod_mpoly_ctx_t pctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp_sp(T->bits, ctx->minfo);
    slong offset, shift;
    slong vi;
    fmpz * Tcoeff;
    ulong * Texp;
    slong Ti;
    n_poly_struct * Acoeff = A->coeffs;
    slong Alen = A->length;
    ulong * Aexp = A->exps;
    slong Ai;
    fmpz * Fcoeff = F->coeffs;
    slong Flen = F->length;
    ulong * Fexp = F->exps;
    slong Fi;
    fmpz_t zero;

    fmpz_init(zero);

    FLINT_ASSERT(T->bits == A->bits);
    FLINT_ASSERT(F->bits == A->bits);
    FLINT_ASSERT(pctx->minfo->nvars == ctx->minfo->nvars);

    FLINT_ASSERT(fmpz_mpoly_is_canonical(F, ctx));
    FLINT_ASSERT(nmod_mpolyn_is_canonical(A, pctx));

    mpoly_gen_offset_shift_sp(&offset, &shift, pctx->minfo->nvars - 1, A->bits, ctx->minfo);

    Flen = F->length;

    fmpz_mpoly_fit_length(T, FLINT_MAX(Flen, Alen), ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    Fi = Ai = vi = 0;
    if (Ai < Alen)
        vi = n_poly_degree(A->coeffs + Ai);

    while (Fi < Flen || Ai < Alen)
    {
        if (Ti >= T->alloc)
        {
            fmpz_mpoly_fit_length(T, Ti + FLINT_MAX(Flen - Fi, Alen - Ai), ctx);
            Tcoeff = T->coeffs;
            Texp = T->exps;
        }

        if (Fi < Flen && Ai < Alen && mpoly_monomial_equal_extra(Fexp + N*Fi,
                                          Aexp + N*Ai, N, offset, vi << shift))
        {
            /* F term ok, A term ok */
            fmpz_CRT_ui(Tcoeff + Ti, Fcoeff + Fi, modulus,
                                    (Acoeff + Ai)->coeffs[vi], pctx->mod.n, 1);
            changed |= !fmpz_equal(Tcoeff + Ti, Fcoeff + Fi);
            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
            do {
                vi--;
            } while (vi >= 0 && (Acoeff + Ai)->coeffs[vi] == 0);
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                    vi = n_poly_degree(A->coeffs + Ai);
            }
        }
        else if (Fi < Flen && (Ai >= Alen || mpoly_monomial_gt_nomask_extra(
                            Fexp + N*Fi, Aexp + N*Ai, N, offset, vi << shift)))
        {
            /* F term ok, A term missing */
            fmpz_CRT_ui(Tcoeff + Ti, Fcoeff + Fi, modulus, 0, pctx->mod.n, 1);
            changed |= !fmpz_equal(Tcoeff + Ti, Fcoeff + Fi);

            mpoly_monomial_set(Texp + N*Ti, Fexp + N*Fi, N);

            Fi++;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen && (Fi >= Flen ||
                     mpoly_monomial_lt_nomask_extra(Fexp + N*Fi, Aexp + N*Ai,
                                                     N, offset, vi << shift)));

            /* F term missing, A term ok */
            fmpz_CRT_ui(Tcoeff + Ti, zero, modulus, (Acoeff + Ai)->coeffs[vi],
                                                               pctx->mod.n, 1);            
            FLINT_ASSERT(!fmpz_is_zero(Tcoeff + Ti));
            changed = 1;
            mpoly_monomial_set_extra(Texp + N*Ti, Aexp + N*Ai, N, offset, vi << shift);

            do {
                vi--;
            } while (vi >= 0 && (Acoeff + Ai)->coeffs[vi] == 0);
            if (vi < 0)
            {
                Ai++;
                if (Ai < Alen)
                    vi = n_poly_degree(A->coeffs + Ai);
            }
        }

        FLINT_ASSERT(!fmpz_is_zero(Tcoeff + Ti));
        Ti++;
    }
    T->length = Ti;

    if (changed)
    {
        fmpz_mpoly_swap(F, T, ctx);
    }

    fmpz_clear(zero);
    return changed;
}

