/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "fmpz_mod_mpoly_factor.h"
#include "fmpz_mod_vec.h"
#include "n_poly.h"
#include "nmod_mpoly_factor.h"
#include "ulong_extras.h"

typedef struct {
    nmod_berlekamp_massey_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} nmod_bma_mpoly_struct;

typedef nmod_bma_mpoly_struct nmod_bma_mpoly_t[1];

typedef struct {
    fmpz_mod_berlekamp_massey_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
    slong pointcount;
} fmpz_mod_bma_mpoly_struct;

typedef fmpz_mod_bma_mpoly_struct fmpz_mod_bma_mpoly_t[1];

typedef struct
{
    slong * degbounds;
    ulong * subdegs;
    fmpz_mod_discrete_log_pohlig_hellman_t dlogenv;
    nmod_discrete_log_pohlig_hellman_t dlogenv_sp;
} mpoly_bma_interpolate_ctx_struct;
typedef mpoly_bma_interpolate_ctx_struct mpoly_bma_interpolate_ctx_t[1];


static void mpoly_bma_interpolate_ctx_init(mpoly_bma_interpolate_ctx_t I, slong nvars)
{
    I->degbounds = (slong *) flint_malloc(nvars*sizeof(slong));
    I->subdegs   = (ulong *) flint_malloc(nvars*sizeof(ulong));
    fmpz_mod_discrete_log_pohlig_hellman_init(I->dlogenv);
    nmod_discrete_log_pohlig_hellman_init(I->dlogenv_sp);
}

static void mpoly_bma_interpolate_ctx_clear(mpoly_bma_interpolate_ctx_t I)
{
    flint_free(I->degbounds);
    flint_free(I->subdegs);
    fmpz_mod_discrete_log_pohlig_hellman_clear(I->dlogenv);
    nmod_discrete_log_pohlig_hellman_clear(I->dlogenv_sp);
}

static void mpoly2_nmod_monomial_evals(
    n_polyun_t EH,
    const ulong * Aexps, flint_bitcnt_t Abits, ulong * Amarks, slong Amarkslen,
    n_poly_struct * caches,
    const mpoly_ctx_t mctx,
    nmod_t fpctx)
{
    slong start, stop, i, j, k, n;
    slong e0, e1;
    slong nvars = mctx->nvars;
    mp_limb_t * p;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    TMP_INIT;

    FLINT_ASSERT(nvars > 2);
    FLINT_ASSERT(Amarkslen > 0);

    TMP_START;

    off = (slong *) TMP_ALLOC(2*nvars*sizeof(slong));
    shift = off + nvars;
    for (k = 0; k < nvars; k++)
        mpoly_gen_offset_shift_sp(&off[k], &shift[k], k, Abits, mctx);

    n_polyun_fit_length(EH, Amarkslen);

    for (i = 0; i < Amarkslen; i++)
    {
        start = Amarks[i];
        stop = Amarks[i + 1];
        FLINT_ASSERT(start < stop);
        n = stop - start;

        e0 = (Aexps[N*start + off[0]] >> shift[0]) & mask;
        e1 = (Aexps[N*start + off[1]] >> shift[1]) & mask;

        EH->exps[i] = pack_exp2(e0, e1);
        n_poly_fit_length(EH->coeffs + i, n);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            p[j] = 1;
            for (k = 2; k < nvars; k++)
            {
                ulong ei = (Aexps[N*(start + j) + off[k]] >> shift[k]) & mask;
                p[j] = nmod_pow_cache_mulpow_ui(p[j], ei, caches + 3*k + 0,
                                    caches + 3*k + 1, caches + 3*k + 2, fpctx);
            }
        }
    }

    EH->length = Amarkslen;

    TMP_END;
}


/*
    evaluation at
        gen(start) -> caches[0]
        gen(start+1) -> caches[1]
        ...
        gen(stop-1) -> caches[stop-start-1]

    the other gen are assumed to not appear in A
*/
static void mpoly_nmod_monomial_evals(
    n_poly_t EH,
    const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
    n_poly_struct * caches,
    slong start,
    slong stop,
    const mpoly_ctx_t mctx,
    nmod_t fpctx)
{
    slong i, k;
    mp_limb_t * p;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    slong num = stop - start;
    TMP_INIT;

    TMP_START;

    off = TMP_ARRAY_ALLOC(2*num, slong);
    shift = off + num;
    for (k = 0; k < num; k++)
        mpoly_gen_offset_shift_sp(&off[k], &shift[k], k + start, Abits, mctx);

    n_poly_fit_length(EH, Alen);
    EH->length = Alen;
    p = EH->coeffs;

    for (i = 0; i < Alen; i++)
    {
        p[i] = 1;
        for (k = 0; k < num; k++)
        {
            ulong ei = (Aexps[N*i + off[k]] >> shift[k]) & mask;
            p[i] = nmod_pow_cache_mulpow_ui(p[i], ei, caches + 3*k + 0,
                                    caches + 3*k + 1, caches + 3*k + 2, fpctx);
        }
    }

    TMP_END;
}


static void fmpz_mpoly2_nmod_coeffs(
    n_polyun_t EH,
    const fmpz * Acoeffs, ulong * Amarks, slong Amarkslen,
    nmod_t fpctx)
{
    slong start, stop, i, n;

    n_polyun_fit_length(EH, Amarkslen);

    for (i = 0; i < Amarkslen; i++)
    {
        start = Amarks[i];
        stop = Amarks[i + 1];
        FLINT_ASSERT(start < stop);
        n = stop - start;

        EH->exps[i] = 0;
        n_poly_fit_length(EH->coeffs + i, n);
        EH->coeffs[i].length = n;
        _fmpz_vec_get_nmod_vec(EH->coeffs[i].coeffs, Acoeffs + start, n, fpctx);
    }

    EH->length = Amarkslen;
}

static void fmpz_mpoly2_fmpz_mod_coeffs(
    fmpz_mod_polyun_t EH,
    const fmpz * Acoeffs, ulong * Amarks, slong Amarkslen,
    const fmpz_mod_ctx_t fpctx)
{
    slong start, stop, i, n;

    fmpz_mod_polyun_fit_length(EH, Amarkslen, fpctx);

    for (i = 0; i < Amarkslen; i++)
    {
        start = Amarks[i];
        stop = Amarks[i + 1];
        FLINT_ASSERT(start < stop);
        n = stop - start;

        EH->exps[i] = 0;
        fmpz_mod_poly_fit_length(EH->coeffs + i, n, fpctx);
        EH->coeffs[i].length = n;
        _fmpz_mod_vec_set_fmpz_vec(EH->coeffs[i].coeffs, Acoeffs + start, n, fpctx);
    }

    EH->length = Amarkslen;
}



static void fmpz_mpoly_nmod_coeffs(
    n_poly_t EH,
    const fmpz * Acoeffs, slong Alen,
    nmod_t fpctx)
{
    n_poly_fit_length(EH, Alen);
    EH->length = Alen;
    _fmpz_vec_get_nmod_vec(EH->coeffs, Acoeffs, Alen, fpctx);
}

static void fmpz_mpoly_fmpz_mod_coeffs(
    fmpz_mod_poly_t EH,
    const fmpz * Acoeffs, slong Alen,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mod_poly_fit_length(EH, Alen, fpctx);
    EH->length = Alen;
    _fmpz_mod_vec_set_fmpz_vec(EH->coeffs, Acoeffs, Alen, fpctx);
}

mp_limb_t n_poly_mod_zip_eval_cur_inc_coeff(
    n_poly_t Acur,
    n_poly_t Ainc,
    n_poly_t Acoeff,
    nmod_t ctx)
{
    return _nmod_zip_eval_step(Acur->coeffs, Ainc->coeffs, Acoeff->coeffs,
                                                            Acur->length, ctx);
}

void fmpz_mod_poly_zip_eval_cur_inc_coeff(
    fmpz_t e,
    fmpz_mod_poly_t Acur,
    fmpz_mod_poly_t Ainc,
    fmpz_mod_poly_t Acoeff,
    const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_zip_eval_step(e, Acur->coeffs, Ainc->coeffs, Acoeff->coeffs,
                                                            Acur->length, ctx);
}


static void n_polyun_mod_zip_eval_cur_inc_coeff(
    n_polyun_t E,
    n_polyun_t Acur,
    const n_polyun_t Ainc,
    const n_polyun_t Acoeff,
    nmod_t ctx)
{
    slong i, Ei;
    slong e0, e1;
    mp_limb_t c;
    n_poly_struct * Ec;

    FLINT_ASSERT(Acur->length > 0);
    FLINT_ASSERT(Acur->length == Ainc->length);
    FLINT_ASSERT(Acur->length == Acoeff->length);

    e0 = extract_exp(Acur->exps[0], 1, 2);
    e1 = extract_exp(Acur->exps[0], 0, 2);

    n_polyun_fit_length(E, 4);
    Ei = 0;
    E->exps[Ei] = e1;
    Ec = E->coeffs + Ei;
    n_poly_zero(Ec);

    for (i = 0; i < Acur->length; i++)
    {
        slong this_len = Acur->coeffs[i].length;
        FLINT_ASSERT(this_len == Ainc->coeffs[i].length);
        FLINT_ASSERT(this_len == Acoeff->coeffs[i].length);

        c = _nmod_zip_eval_step(Acur->coeffs[i].coeffs, Ainc->coeffs[i].coeffs,
                                      Acoeff->coeffs[i].coeffs, this_len, ctx);

        e0 = extract_exp(Acur->exps[i], 1, 2);
        e1 = extract_exp(Acur->exps[i], 0, 2);

        if (E->exps[Ei] != e0)
        {
            n_polyun_fit_length(E, Ei + 2);
            Ei += !n_poly_is_zero(E->coeffs + Ei);
            E->exps[Ei] = e0;
            Ec = E->coeffs + Ei;
            n_poly_zero(Ec);
        }

        n_poly_set_coeff(Ec, e1, c);
    }

    Ei += !n_poly_is_zero(E->coeffs + Ei);
    E->length = Ei;

    FLINT_ASSERT(n_polyun_mod_is_canonical(E, ctx));
}


static void fmpz_mpoly2_eval_nmod(
    n_polyun_t E,
    n_polyun_t EHcur,
    n_polyun_t EHinc,
    n_polyun_t EHcoeffs,
    const fmpz_mpoly_t A, ulong * Amarks, slong Amarkslen,
    n_poly_struct * alpha_caches,
    const fmpz_mpoly_ctx_t ctx,
    nmod_t fpctx)
{
    mpoly2_nmod_monomial_evals(EHcur, A->exps, A->bits, Amarks, Amarkslen,
                                              alpha_caches, ctx->minfo, fpctx);
    n_polyun_set(EHinc, EHcur);
    fmpz_mpoly2_nmod_coeffs(EHcoeffs, A->coeffs, Amarks, Amarkslen, fpctx);
    n_polyun_mod_zip_eval_cur_inc_coeff(E, EHcur, EHinc, EHcoeffs, fpctx);
}

static void fmpz_mpoly2_eval_fmpz_mod(
    fmpz_mod_polyun_t E,
    fmpz_mod_polyun_t EHcur,
    fmpz_mod_polyun_t EHinc,
    fmpz_mod_polyun_t EHcoeffs,
    const fmpz_mpoly_t A, ulong * Amarks, slong Amarkslen,
    fmpz_mod_poly_struct * alpha_caches,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    mpoly2_monomial_evals_fmpz_mod(EHcur, A->exps, A->bits, Amarks, Amarkslen,
                       alpha_caches + 2, ctx->minfo->nvars, ctx->minfo, fpctx);
    fmpz_mod_polyun_set(EHinc, EHcur, fpctx);
    fmpz_mpoly2_fmpz_mod_coeffs(EHcoeffs, A->coeffs, Amarks, Amarkslen, fpctx);
    fmpz_mod_polyu2n_zip_eval_cur_inc_coeff(E, EHcur, EHinc, EHcoeffs, fpctx);
}


/*
    set out to the evaluation of variables after ksub at alpha^w

    out[m+0]   = alpha ^ (w * subdegs[n-1] * subdegs[n-2] * ... * * subdegs[m+1])
      ...
    out[n-3] = alpha ^ (w * subdegs[n-1] * subdegs[n-2])
    out[n-2] = alpha ^ (w * subdegs[n-1])
    out[n-1] = alpha ^ (w)

    secret: subdegs[0] is not used
*/

void nmod_mpoly_bma_interpolate_alpha_powers(
    mp_limb_t * out,
    ulong w,
    slong m,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mpoly_ctx_t ctx,
    nmod_t fpctx)
{
    slong j = ctx->minfo->nvars - 1;
    out[j] = nmod_pow_ui(Ictx->dlogenv_sp->alpha, w, fpctx);
    for (; j > m; j--)
        out[j - 1] = nmod_pow_ui(out[j], Ictx->subdegs[j], fpctx);
}

void fmpz_mod_mpoly_bma_interpolate_alpha_powers(
    fmpz * out,
    const fmpz_t w,
    slong m,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong j = ctx->minfo->nvars - 1;
    fmpz_mod_pow_fmpz(out + j, Ictx->dlogenv->alpha, w, fpctx);
    for (; j > m; j--)
        fmpz_mod_pow_ui(out + j - 1, out + j, Ictx->subdegs[j], fpctx);
}


void nmod_bma_mpoly_init(nmod_bma_mpoly_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->exps = NULL;
    A->coeffs = NULL;
    A->pointcount = 0;
}

void nmod_bma_mpoly_reset_prime(
    nmod_bma_mpoly_t A,
    nmod_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_berlekamp_massey_set_prime(A->coeffs + i, fpctx.n);
}


void nmod_bma_mpoly_clear(nmod_bma_mpoly_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        nmod_berlekamp_massey_clear(A->coeffs + i);
    }

    if (A->exps)
        flint_free(A->exps);
    if (A->coeffs)
        flint_free(A->coeffs);
}


void nmod_bma_mpoly_fit_length(
    nmod_bma_mpoly_t A,
    slong length,
    nmod_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        A->exps = FLINT_ARRAY_REALLOC(A->exps, new_alloc, ulong);
        A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, new_alloc,
                                                 nmod_berlekamp_massey_struct);
        for (i = old_alloc; i < new_alloc; i++)
            nmod_berlekamp_massey_init(A->coeffs + i, fpctx.n);

        A->alloc = new_alloc;
    }
}

void nmod_bma_mpoly_zero(nmod_bma_mpoly_t L)
{
    L->length = 0;
    L->pointcount = 0;
}


int nmod_bma_mpoly_reduce(nmod_bma_mpoly_t L)
{
    slong i;
    int changed;

    changed = 0;

    for (i = 0; i < L->length; i++)
    {
        FLINT_ASSERT(L->pointcount == nmod_berlekamp_massey_point_count(L->coeffs + i));
        changed |= nmod_berlekamp_massey_reduce(L->coeffs + i);
    }

    return changed;
}


void nmod_bma_mpoly_add_point(
    nmod_bma_mpoly_t L,
    const n_polyun_t A,
    nmod_t fpctx)
{
    slong j;
    slong Alen = A->length;
    slong Li, Ai, ai;
    nmod_berlekamp_massey_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;
    ulong Aexp;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        for (ai = A->coeffs[Ai].length - 1; ai >= 0; ai--)
            tot += (0 != A->coeffs[Ai].coeffs[ai]);
        nmod_bma_mpoly_fit_length(L, tot, fpctx);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = n_poly_degree(A->coeffs + Ai);
        Aexp = pack_exp2(A->exps[Ai], ai);
    }

    while (Li < Llen || Ai < Alen)
    {
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L term present, A term present */
add_same_exp:
            nmod_berlekamp_massey_add_point(Lcoeff + Li, A->coeffs[Ai].coeffs[ai]);
            Li++;

            do {
                ai--;
            } while (ai >= 0 && A->coeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = n_poly_degree(A->coeffs + Ai);
                    Aexp = pack_exp2(A->exps[Ai], ai);
                }
            }
            else
            {
                FLINT_ASSERT(Ai < A->length);
                Aexp = pack_exp2(A->exps[Ai], ai);
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            /* L term present, A term missing */
            nmod_berlekamp_massey_add_zeros(Lcoeff + Li, 1);
            Li++;
        }
        else
        {
            /* L term missing, A term present */
            FLINT_ASSERT(Ai < Alen && (Li >= Llen || Lexp[Li] < Aexp));
            {
                ulong texp;
                nmod_berlekamp_massey_struct tcoeff;

                nmod_bma_mpoly_fit_length(L, Llen + 1, fpctx);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Llen - 1; j >= Li; j--)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            nmod_berlekamp_massey_start_over(Lcoeff + Li);
            nmod_berlekamp_massey_add_zeros(Lcoeff + Li, L->pointcount);
            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;

            goto add_same_exp;
        }
    }

    L->pointcount++;
}


static int n_polyu2n_add_zipun_must_match(
    n_polyun_t Z,
    const n_polyun_t A,
    slong cur_length)
{
    slong i, Ai, ai;
    ulong Aexp;
    slong Alen = A->length;

    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = n_poly_degree(A->coeffs + Ai);
        Aexp = pack_exp2(A->exps[Ai], ai);
    }

    for (i = 0; i < Z->length; i++)
    {
        if (Ai < Alen && Z->exps[i] == Aexp)
        {
            /* Z present, A present */
            Z->coeffs[i].coeffs[cur_length] = A->coeffs[Ai].coeffs[ai];
            Z->coeffs[i].length = cur_length + 1;

            do {
                ai--;
            } while (ai >= 0 && A->coeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = n_poly_degree(A->coeffs + Ai);
                    Aexp = pack_exp2(A->exps[Ai], ai);
                }
            }
            else
            {
                FLINT_ASSERT(Ai < A->length);
                Aexp = pack_exp2(A->exps[Ai], ai);
            }
        }
        else if (Ai > Alen || Z->exps[i] > Aexp)
        {
            /* Z present, A missing */
            Z->coeffs[i].coeffs[cur_length] = 0;
            Z->coeffs[i].length = cur_length + 1;
        }
        else
        {
            /* Z missing, A present */
            return 0;
        }
    }

    return 1;
}


/*
    If I was formed from evaluations at
        alpha^alphashift, alpha^(alphashift + 1), ...
    construct the corresponding mpoly if possible with coeffs in (-p/2, p/2]
    The substitution degrees and degree bounds in Ictx are used.
*/
static int _nmod_mpoly_bma_get_fmpz_mpoly2(
    fmpz * Acoeffs,
    ulong * Aexps,
    flint_bitcnt_t Abits,
    ulong e0, ulong e1,
    const mpoly_ctx_t mctx,
    slong * shifts, slong * offsets,
    ulong alphashift,
    nmod_berlekamp_massey_t I,
    const mpoly_bma_interpolate_ctx_t Ictx,
    nmod_t fpctx)
{
    int success;
    slong i, j, t;
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    ulong new_exp, this_exp;
    mp_limb_t * values, * roots;
    mp_limb_t T, S, V, V0, V1, V2, p0, p1, r;

    FLINT_ASSERT(mctx->ord == ORD_LEX);

    t = nmod_poly_degree(I->V1);
    FLINT_ASSERT(I->points->length >= t);

    if (t < 1)
        return 0;

    /* use the rt member of I as temp space for roots - slightly dirty */
    nmod_poly_fit_length(I->rt, t);
    I->rt->length = t;
    roots = I->rt->coeffs;
    values = I->points->coeffs;

    success = nmod_poly_find_distinct_nonzero_roots(roots, I->V1);
    if (!success)
        return 0;

    for (i = 0; i < t; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = roots[i];
        for (j = t; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, fpctx), I->V1->coeffs[j], fpctx);
            S = nmod_add(nmod_mul(r, S, fpctx), T, fpctx);
            umul_ppmm(p1, p0, values[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, fpctx), I->V1->coeffs[0], fpctx) == 0);
        NMOD_RED3(V, V2, V1, V0, fpctx);
        S = nmod_mul(S, nmod_pow_ui(r, alphashift, fpctx), fpctx);
        V0 = nmod_mul(V, nmod_inv(S, fpctx), fpctx);
        if (V0 == 0)
            return 0;

        if (fpctx.n - V0 < V0)
            fmpz_neg_ui(Acoeffs + i, fpctx.n - V0);
        else
            fmpz_set_ui(Acoeffs + i, V0);

        mpoly_monomial_zero(Aexps + N*i, N);
        (Aexps + N*i)[offsets[0]] |= e0 << shifts[0];
        (Aexps + N*i)[offsets[1]] |= e1 << shifts[1];

        new_exp = nmod_discrete_log_pohlig_hellman_run(Ictx->dlogenv_sp, roots[i]);
        for (j = mctx->nvars - 1; j >= 2; j--)
        {
            this_exp = new_exp % Ictx->subdegs[j];
            new_exp = new_exp / Ictx->subdegs[j];
            if (this_exp > Ictx->degbounds[j])
                return 0;
            (Aexps + N*i)[offsets[j]] |= this_exp << shifts[j];
        }
        if (new_exp != 0)
            return 0;
    }

    return 1;
}

static int _fmpz_mod_bma_get_fmpz_mpoly2(
    fmpz * Acoeffs,
    ulong * Aexps,
    flint_bitcnt_t Abits,
    ulong e0, ulong e1,
    const mpoly_ctx_t mctx,
    slong * shifts, slong * offsets,
    const fmpz_t alphashift,
    fmpz_mod_berlekamp_massey_t I,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mod_ctx_t fpctx)
{
    int success;
    slong i, j, t;
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    ulong this_exp;
    fmpz_t new_exp;
    fmpz * values, * roots;
    fmpz_t T, S, V, temp, halfp;

    fmpz_init(halfp);
    fmpz_init(T);
    fmpz_init(S);
    fmpz_init(V);
    fmpz_init(temp);
    fmpz_init(new_exp);

    fmpz_tdiv_q_2exp(halfp, fmpz_mod_ctx_modulus(fpctx), 1);

    t = fmpz_mod_poly_degree(I->V1, fpctx);
    FLINT_ASSERT(I->points->length >= t);

    if (t < 1)
    {
        success = 0;
        goto cleanup;
    }

    fmpz_mod_poly_fit_length(I->rt, t, fpctx);
    I->rt->length = t;
    roots = I->rt->coeffs;
    values = I->points->coeffs;

    success = fmpz_mod_poly_find_distinct_nonzero_roots(roots, I->V1, fpctx);
    if (!success)
        goto cleanup;

    for (i = 0; i < t; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        fmpz_zero(V);
        fmpz_zero(T);
        fmpz_zero(S);
        for (j = t; j > 0; j--)
        {
            fmpz_mod_mul(temp, roots + i, T, fpctx);
            fmpz_mod_add(T, temp, I->V1->coeffs + j, fpctx);
            fmpz_mod_mul(temp, roots + i, S, fpctx);
            fmpz_mod_add(S, temp, T, fpctx);
            fmpz_mod_mul(temp, values + j - 1, T, fpctx);
            fmpz_mod_add(V, V, temp, fpctx);
        }
        /* roots[i] should be a root of master */
#if FLINT_WANT_ASSERT
        fmpz_mod_mul(temp, roots + i, T, fpctx);
        fmpz_mod_add(temp, temp, I->V1->coeffs + 0, fpctx);
        FLINT_ASSERT(fmpz_is_zero(temp));
#endif
        fmpz_mod_pow_fmpz(temp, roots + i, alphashift, fpctx);
        fmpz_mod_mul(S, S, temp, fpctx);
        fmpz_mod_inv(temp, S, fpctx);
        fmpz_mod_mul(Acoeffs + i, V, temp, fpctx);
        if (fmpz_is_zero(Acoeffs + i))
        {
            success = 0;
            goto cleanup;
        }

        if (fmpz_cmp(Acoeffs + i, halfp) > 0)
            fmpz_sub(Acoeffs + i, Acoeffs + i, fmpz_mod_ctx_modulus(fpctx));

        mpoly_monomial_zero(Aexps + N*i, N);
        (Aexps + N*i)[offsets[0]] |= e0 << shifts[0];
        (Aexps + N*i)[offsets[1]] |= e1 << shifts[1];

        fmpz_mod_discrete_log_pohlig_hellman_run(new_exp, Ictx->dlogenv, roots + i);
        for (j = mctx->nvars - 1; j >= 2; j--)
        {
            this_exp = fmpz_fdiv_ui(new_exp, Ictx->subdegs[j]);
            fmpz_fdiv_q_ui(new_exp, new_exp, Ictx->subdegs[j]);
            if (this_exp > Ictx->degbounds[j])
            {
                success = 0;
                goto cleanup;
            }
            (Aexps + N*i)[offsets[j]] |= this_exp << shifts[j];
        }
        if (!fmpz_is_zero(new_exp))
        {
            success = 0;
            goto cleanup;
        }
    }

    success = 1;

cleanup:

    fmpz_clear(T);
    fmpz_clear(S);
    fmpz_clear(V);
    fmpz_clear(temp);
    fmpz_clear(halfp);

    return success;
}




int nmod_bma_mpoly_get_fmpz_mpoly2(
    fmpz_mpoly_t A,
    n_poly_t Amarks,
    const fmpz_mpoly_ctx_t ctx,
    ulong alphashift,
    const nmod_bma_mpoly_t L,
    const mpoly_bma_interpolate_ctx_t Ictx,
    nmod_t fpctx)
{
    int success;
    slong i, j, k;
    slong * shifts, * offsets;
    ulong * marks;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    TMP_INIT;

    if (L->length < 1)
        return 0;

    n_poly_fit_length(Amarks, L->length + 1);
    Amarks->length = L->length;
    marks = Amarks->coeffs;

    TMP_START;
    shifts = TMP_ARRAY_ALLOC(2*ctx->minfo->nvars, slong);
    offsets = shifts + ctx->minfo->nvars;

    for (j = 0; j < ctx->minfo->nvars; j++)
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);

    k = 0;
    for (i = 0; i < L->length; i++)
    {
        nmod_berlekamp_massey_reduce(L->coeffs + i);
        marks[i] = k;
        k += nmod_poly_degree(L->coeffs[i].V1);
    }
    marks[L->length] = k;

    fmpz_mpoly_fit_length(A, k, ctx);
    A->length = k;

    for (i = 0; i < L->length; i++)
    {
        ulong e0 = extract_exp(L->exps[i], 1, 2);
        ulong e1 = extract_exp(L->exps[i], 0, 2);

        success = _nmod_mpoly_bma_get_fmpz_mpoly2(A->coeffs + marks[i],
                      A->exps + N*marks[i], A->bits, e0, e1, ctx->minfo,
                      shifts, offsets, alphashift, L->coeffs + i, Ictx, fpctx);
        if (!success)
            goto cleanup;
    }

    fmpz_mpoly_sort_terms(A, ctx);

    success = 1;

cleanup:

    TMP_END;

    return success;
}


int fmpz_mod_bma_mpoly_get_fmpz_mpoly2(
    fmpz_mpoly_t A,
    n_poly_t Amarks,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_t alphashift,
    const fmpz_mod_bma_mpoly_t L,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mod_ctx_t fpctx)
{
    int success;
    slong i, j, k;
    slong * shifts, * offsets;
    ulong * marks;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    TMP_INIT;

    if (L->length < 1)
        return 0;

    n_poly_fit_length(Amarks, L->length + 1);
    Amarks->length = L->length;
    marks = Amarks->coeffs;

    TMP_START;
    shifts = TMP_ARRAY_ALLOC(2*ctx->minfo->nvars, slong);
    offsets = shifts + ctx->minfo->nvars;

    for (j = 0; j < ctx->minfo->nvars; j++)
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);

    k = 0;
    for (i = 0; i < L->length; i++)
    {
        fmpz_mod_berlekamp_massey_reduce(L->coeffs + i, fpctx);
        marks[i] = k;
        k += fmpz_mod_poly_degree(L->coeffs[i].V1, fpctx);
    }
    marks[L->length] = k;

    fmpz_mpoly_fit_length(A, k, ctx);
    A->length = k;

    for (i = 0; i < L->length; i++)
    {
        ulong e0 = extract_exp(L->exps[i], 1, 2);
        ulong e1 = extract_exp(L->exps[i], 0, 2);

        success = _fmpz_mod_bma_get_fmpz_mpoly2(A->coeffs + marks[i],
                      A->exps + N*marks[i], A->bits, e0, e1, ctx->minfo,
                      shifts, offsets, alphashift, L->coeffs + i, Ictx, fpctx);
        if (!success)
            goto cleanup;
    }

    fmpz_mpoly_sort_terms(A, ctx);

    success = 1;

cleanup:

    TMP_END;

    return success;
}




void fmpz_mod_bma_mpoly_init(fmpz_mod_bma_mpoly_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->exps = NULL;
    A->coeffs = NULL;
    A->pointcount = 0;
}


void fmpz_mod_bma_mpoly_clear(
    fmpz_mod_bma_mpoly_t A,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;

    for (i = 0; i < A->alloc; i++)
        fmpz_mod_berlekamp_massey_clear(A->coeffs + i, fpctx);

    flint_free(A->exps);
    flint_free(A->coeffs);
}


void fmpz_mod_bma_mpoly_fit_length(
    fmpz_mod_bma_mpoly_t A,
    slong length,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        A->exps = FLINT_ARRAY_REALLOC(A->exps, new_alloc, ulong);
        A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, new_alloc,
                                             fmpz_mod_berlekamp_massey_struct);

        for (i = old_alloc; i < new_alloc; i++)
            fmpz_mod_berlekamp_massey_init(A->coeffs + i, fpctx);

        A->alloc = new_alloc;
    }
}

void fmpz_mod_bma_mpoly_zero(fmpz_mod_bma_mpoly_t L)
{
    L->length = 0;
    L->pointcount = 0;
}

int fmpz_mod_bma_mpoly_reduce(fmpz_mod_bma_mpoly_t L, const fmpz_mod_ctx_t fpctx)
{
    slong i;
    int changed;

    changed = 0;

    for (i = 0; i < L->length; i++)
    {
        FLINT_ASSERT(L->pointcount == fmpz_mod_berlekamp_massey_point_count(L->coeffs + i));
        changed |= fmpz_mod_berlekamp_massey_reduce(L->coeffs + i, fpctx);
    }

    return changed;
}

static void fmpz_mod_bma_mpoly_add_point(
    fmpz_mod_bma_mpoly_t L,
    const fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx_mp)
{
    slong j;
    slong Alen = A->length;
    fmpz_mod_poly_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    fmpz_mod_berlekamp_massey_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
            for (ai = Acoeff[Ai].length - 1; ai >= 0; ai--)
                tot += !fmpz_is_zero(Acoeff[Ai].coeffs + ai);

        fmpz_mod_bma_mpoly_fit_length(L, tot, ctx_mp);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = fmpz_mod_poly_degree(A->coeffs + Ai, ctx_mp);
        Aexp = pack_exp2(A->exps[Ai], ai);
    }

    while (Li < Llen || Ai < Alen)
    {
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L term present, A term present */
add_same_exp:
            fmpz_mod_berlekamp_massey_add_point(Lcoeff + Li,
                                               Acoeff[Ai].coeffs + ai, ctx_mp);
            Li++;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Acoeff[Ai].coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = fmpz_mod_poly_degree(A->coeffs + Ai, ctx_mp);
                    Aexp = pack_exp2(A->exps[Ai], ai);
                }
            }
            else
            {
                FLINT_ASSERT(Ai < Alen);
                Aexp = pack_exp2(A->exps[Ai], ai);
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            /* L term present, A term missing */
            fmpz_mod_berlekamp_massey_add_zeros(Lcoeff + Li, 1, ctx_mp);
            Li++;
        }
        else
        {
            /* L term missing, A term present */
            FLINT_ASSERT(Ai < Alen && (Li >= Llen || Lexp[Li] < Aexp));
            {
                ulong texp;
                fmpz_mod_berlekamp_massey_struct tcoeff;

                fmpz_mod_bma_mpoly_fit_length(L, Llen + 1, ctx_mp);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Llen - 1; j >= Li; j--)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            fmpz_mod_berlekamp_massey_start_over(Lcoeff + Li, ctx_mp);
            fmpz_mod_berlekamp_massey_add_zeros(Lcoeff + Li, L->pointcount, ctx_mp);
            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;
            goto add_same_exp;
        }
    }

    L->pointcount++;
}


/*
    A is in ZZ[x_0, ..., x_(n-1)]

    After the substitutions
        x_0     = x ^ (sub[1] * sub[2] * ... * sub[n-1])

        x_(n-2) = x ^ (sub[n-1])
        x_(n-1) = x ^ (1)
    a univariate in ZZ[x] remains. Return the content of this poly.
*/
void _fmpz_mpoly_ksub_content(
    fmpz_t content,
    const fmpz * Acoeffs,
    const ulong * Aexps,
    slong Alength,
    flint_bitcnt_t Abits,
    const ulong * subdegs,
    const mpoly_ctx_t mctx)
{
    slong i, j;
    slong nvars = mctx->nvars;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * offsets, * shifts;
    fmpz_mpoly_t T;
    fmpz_mpoly_ctx_t Tctx;
    fmpz_t e;
    TMP_INIT;

    TMP_START;
    fmpz_init(e);

    fmpz_mpoly_ctx_init(Tctx, 1, ORD_LEX);
    fmpz_mpoly_init(T, Tctx);

    offsets = TMP_ALLOC(2*nvars*sizeof(slong));
    shifts = offsets + nvars;
    for (j = 2; j < nvars; j++)
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, Abits, mctx);

    for (i = 0; i < Alength; i++)
    {
        fmpz_zero(e);
        for (j = 2; j < mctx->nvars; j++)
        {
            fmpz_mul_ui(e, e, subdegs[j]);
            fmpz_add_ui(e, e, ((Aexps + N*i)[offsets[j]]>>shifts[j])&mask);
        }
        _fmpz_mpoly_push_exp_ffmpz(T, e, Tctx);
        fmpz_set(T->coeffs + T->length - 1, Acoeffs + i);
    }

    fmpz_mpoly_sort_terms(T, Tctx);
    fmpz_mpoly_combine_like_terms(T, Tctx);

    _fmpz_vec_content(content, T->coeffs, T->length);

    fmpz_mpoly_clear(T, Tctx);
    fmpz_mpoly_ctx_clear(Tctx);

    fmpz_clear(e);
    TMP_END;
}




/*
    return 1: good (or at least not bad)
           0: no match
          -1: better degree bound is in GevaldegXY
*/

int static _random_check_sp(
    ulong * GevaldegXY,
    ulong GdegboundXY,
    int which_check,
    n_polyun_t Aeval_sp, n_polyun_t Acur_sp, n_polyun_t Ainc_sp, n_polyun_t Acoeff_sp,
    n_polyun_t Beval_sp, n_polyun_t Bcur_sp, n_polyun_t Binc_sp, n_polyun_t Bcoeff_sp,
    n_polyun_t Geval_sp,
    n_polyun_t Abareval_sp,
    n_polyun_t Bbareval_sp,
    mp_limb_t * alphas_sp,
    n_poly_struct * alpha_caches_sp,
    const fmpz_mpoly_t H, n_poly_t Hmarks,
    const fmpz_mpoly_t A, n_poly_t Amarks, ulong Abidegree,
    const fmpz_mpoly_t B, n_poly_t Bmarks, ulong Bbidegree,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx,
    nmod_t ctx_sp,
    flint_rand_t randstate,
    n_poly_polyun_stack_t St_sp)
{
    mp_limb_t Gammaeval_sp;
    int success;
    int point_try_count;
    slong i;

    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        alphas_sp[0] = alphas_sp[1] = 0;
        for (i = 2; i < ctx->minfo->nvars; i++)
        {
            alphas_sp[i] = n_urandint(randstate, ctx_sp.n - 1) + 1;
            nmod_pow_cache_start(alphas_sp[i], alpha_caches_sp + 3*i,
                         alpha_caches_sp + 3*i + 1, alpha_caches_sp + 3*i + 2);
        }

        fmpz_mpoly2_eval_nmod(Aeval_sp, Acur_sp, Ainc_sp, Acoeff_sp, A,
                 Amarks->coeffs, Amarks->length, alpha_caches_sp, ctx, ctx_sp);
        fmpz_mpoly2_eval_nmod(Beval_sp, Bcur_sp, Binc_sp, Bcoeff_sp, B,
                 Bmarks->coeffs, Bmarks->length, alpha_caches_sp, ctx, ctx_sp);

        if (Aeval_sp->length < 1 || Beval_sp->length < 1 ||
            n_polyu1n_bidegree(Aeval_sp) != Abidegree ||
            n_polyu1n_bidegree(Beval_sp) != Bbidegree)
        {
            continue;
        }

        Gammaeval_sp = fmpz_mpoly_evaluate_all_nmod(Gamma, alphas_sp, ctx, ctx_sp);
        FLINT_ASSERT(Gammaeval_sp != 0);

        success = n_polyu1n_mod_gcd_brown_smprime(Geval_sp, Abareval_sp,
                               Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp, St_sp);
        if (success)
            continue;

        _n_poly_vec_mul_nmod_intertible(Geval_sp->coeffs, Geval_sp->length,
                                                         Gammaeval_sp, ctx_sp);
        FLINT_ASSERT(Geval_sp->length > 0);
        *GevaldegXY = n_polyu1n_bidegree(Geval_sp);

        if (GdegboundXY < *GevaldegXY)
            continue;
        else if (GdegboundXY > *GevaldegXY)
            return -1;

        if (which_check == 0)
        {
            fmpz_mpoly2_eval_nmod(Bbareval_sp, Bcur_sp, Binc_sp, Bcoeff_sp, H,
                 Hmarks->coeffs, Hmarks->length, alpha_caches_sp, ctx, ctx_sp);

            return n_polyun_equal(Bbareval_sp, Geval_sp);
        }
        else
        {
            fmpz_mpoly2_eval_nmod(Geval_sp, Bcur_sp, Binc_sp, Bcoeff_sp, H,
                 Hmarks->coeffs, Hmarks->length, alpha_caches_sp, ctx, ctx_sp);
    
            return n_polyun_equal(Geval_sp, which_check == 1 ? Abareval_sp : Bbareval_sp);
        }
    }

    return 1; /* Hmm */
}

int static _random_check_mp(
    ulong * GevaldegXY,
    ulong GdegboundXY,
    int which_check,
    fmpz_mod_polyun_t Aeval_mp, fmpz_mod_polyun_t Acur_mp, fmpz_mod_polyun_t Ainc_mp, fmpz_mod_polyun_t Acoeff_mp,
    fmpz_mod_polyun_t Beval_mp, fmpz_mod_polyun_t Bcur_mp, fmpz_mod_polyun_t Binc_mp, fmpz_mod_polyun_t Bcoeff_mp,
    fmpz_mod_polyun_t Geval_mp,
    fmpz_mod_polyun_t Abareval_mp,
    fmpz_mod_polyun_t Bbareval_mp,
    fmpz_t Gammaeval_mp,
    fmpz * alphas_mp,
    fmpz_mod_poly_struct * alpha_caches_mp,
    const fmpz_mpoly_t H, n_poly_t Hmarks,
    const fmpz_mpoly_t A, n_poly_t Amarks, ulong Abidegree,
    const fmpz_mpoly_t B, n_poly_t Bmarks, ulong Bbidegree,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t ctx_mp,
    flint_rand_t randstate,
    fmpz_mod_poly_polyun_stack_t St_mp)
{
    int success;
    int point_try_count;
    slong i;

    /* try to test H at a random evaluation point */
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        for (i = 2; i < ctx->minfo->nvars; i++)
        {
            fmpz_randm(alphas_mp + i, randstate, fmpz_mod_ctx_modulus(ctx_mp));
            fmpz_mod_pow_cache_start(alphas_mp + i, alpha_caches_mp + i, ctx_mp);
        }

        fmpz_mpoly2_eval_fmpz_mod(Aeval_mp, Acur_mp, Ainc_mp, Acoeff_mp, A,
                 Amarks->coeffs, Amarks->length, alpha_caches_mp, ctx, ctx_mp);
        fmpz_mpoly2_eval_fmpz_mod(Beval_mp, Bcur_mp, Binc_mp, Bcoeff_mp, B,
                 Bmarks->coeffs, Bmarks->length, alpha_caches_mp, ctx, ctx_mp);

        if (Aeval_mp->length < 1 || Beval_mp->length < 1 ||
            fmpz_mod_polyu1n_bidegree(Aeval_mp) != Abidegree ||
            fmpz_mod_polyu1n_bidegree(Beval_mp) != Bbidegree)
        {
            continue;
        }

        fmpz_mpoly_evaluate_all_fmpz_mod(Gammaeval_mp, Gamma, alphas_mp, ctx, ctx_mp);
        FLINT_ASSERT(!fmpz_is_zero(Gammaeval_mp));

        success = fmpz_mod_polyu1n_gcd_brown_smprime(Geval_mp, Abareval_mp, Bbareval_mp,
                 Aeval_mp, Beval_mp, ctx_mp, St_mp->poly_stack, St_mp->polyun_stack);
        if (!success)
            continue;

        _fmpz_mod_poly_vec_mul_fmpz_mod(Geval_mp->coeffs, Geval_mp->length,
                                                         Gammaeval_mp, ctx_mp);
        FLINT_ASSERT(Geval_mp->length > 0);
        *GevaldegXY = fmpz_mod_polyu1n_bidegree(Geval_mp);

        if (GdegboundXY < *GevaldegXY)
        {
            continue;
        }
        else if (GdegboundXY > *GevaldegXY)
        {
            return -1;
        }

        if (which_check == 0)
        {
            fmpz_mpoly2_eval_fmpz_mod(Bbareval_mp, Bcur_mp, Binc_mp, Bcoeff_mp, H,
                 Hmarks->coeffs, Hmarks->length, alpha_caches_mp, ctx, ctx_mp);

            return fmpz_mod_polyun_equal(Bbareval_mp, Geval_mp, ctx_mp);
        }
        else
        {
            fmpz_mpoly2_eval_fmpz_mod(Geval_mp, Bcur_mp, Binc_mp, Bcoeff_mp, H,
                 Hmarks->coeffs, Hmarks->length, alpha_caches_mp, ctx, ctx_mp);
    
            return fmpz_mod_polyun_equal(Geval_mp, which_check == 1 ?
                                            Abareval_mp : Bbareval_mp, ctx_mp);
        }
    }

    return 1; /* Hmm */
}


/*
    return 
        -1: singular
        0:  inconsistent
        1:  success
*/
static int zip_solve(
    mp_limb_t * Acoeffs,
    n_polyun_t Z,
    n_polyun_t H,
    n_polyun_t M,
    const nmod_t fpctx)
{
    int success;
    slong Ai, i, n;
    n_poly_t t;

    n_poly_init(t);

    FLINT_ASSERT(Z->length == H->length);
    FLINT_ASSERT(Z->length == M->length);

    Ai = 0;
    for (i = 0; i < H->length; i++)
    {
        n = H->coeffs[i].length;
        FLINT_ASSERT(M->coeffs[i].length == n + 1);
        FLINT_ASSERT(Z->coeffs[i].length >= n);

        n_poly_fit_length(t, n);

        success = _nmod_zip_vand_solve(Acoeffs + Ai,
                                 H->coeffs[i].coeffs, n,
                                 Z->coeffs[i].coeffs, Z->coeffs[i].length,
                                 M->coeffs[i].coeffs, t->coeffs, fpctx);
        if (success < 1)
        {
            n_poly_clear(t);
            return success;
        }

        Ai += n;
    }

    n_poly_clear(t);
    return 1;
}


int _fmpz_vec_crt_nmod(
    flint_bitcnt_t * maxbits_,
    fmpz * a,
    fmpz_t am,
    mp_limb_t * b,
    slong len,
    nmod_t mod)
{
    int changed = 0;
    flint_bitcnt_t bits, maxbits = 0;
    slong i;
    fmpz_t t;

    fmpz_init(t);

    for (i = 0; i < len; i++)
    {
        fmpz_CRT_ui(t, a + i, am, b[i], mod.n, 1);
        changed |= !fmpz_equal(t, a + i);
        bits = fmpz_bits(t);
        maxbits = FLINT_MAX(maxbits, bits);
        fmpz_swap(a + i, t);
    }

    fmpz_clear(t);

    *maxbits_ = maxbits;
    return changed;
}

int fmpz_mpolyl_gcd_zippel2(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx)
{
    int which_check, changed, success, point_try_count;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    mpoly_bma_interpolate_ctx_t Ictx;
    n_poly_t Amarks, Bmarks, Hmarks;
    flint_bitcnt_t Hbitbound = UWORD_MAX;
    flint_bitcnt_t Hbits; /* coeff size */
    fmpz_mpoly_t H;
    fmpz_mpoly_t Hcontent;
    fmpz_t p, pm1;
    /* multi precision workspace */
    fmpz_mod_ctx_t ctx_mp;
    fmpz_mod_poly_polyun_stack_t St_mp;
    fmpz_mod_bma_mpoly_t GLambda_mp, AbarLambda_mp, BbarLambda_mp;
    fmpz_mod_polyun_t Aeval_mp, Beval_mp, Geval_mp, Abareval_mp, Bbareval_mp;
    fmpz_mod_poly_t Gammacur_mp, Gammainc_mp, Gammacoeff_mp;
    fmpz_mod_polyun_t Acur_mp, Ainc_mp, Acoeff_mp;
    fmpz_mod_polyun_t Bcur_mp, Binc_mp, Bcoeff_mp;
    fmpz_t sshift_mp, last_unlucky_sshift_plus_1_mp, image_count_mp;
    fmpz_t Gammaeval_mp;
    fmpz * alphas_mp;
    fmpz_mod_poly_struct * alpha_caches_mp;
    /* single precision workspace */
    nmod_t ctx_sp;
    n_poly_polyun_stack_t St_sp;
    nmod_bma_mpoly_t GLambda_sp, AbarLambda_sp, BbarLambda_sp;
    n_polyun_t Aeval_sp, Beval_sp, Geval_sp, Abareval_sp, Bbareval_sp;
    n_poly_t Gammacur_sp, Gammainc_sp, Gammacoeff_sp;
    n_polyun_t Acur_sp, Ainc_sp, Acoeff_sp;
    n_polyun_t Bcur_sp, Binc_sp, Bcoeff_sp;
    mp_limb_t p_sp, sshift_sp, last_unlucky_sshift_plus_1_sp, image_count_sp;
    mp_limb_t Gammaeval_sp;
    mp_limb_t * alphas_sp;
    n_poly_struct * alpha_caches_sp;
    /* misc */
    n_polyun_t HH, MH, ZH;
    n_poly_t Hn;
    slong i, j;
    ulong GdegboundXY, GevaldegXY;
    slong * Adegs, * Bdegs;
    flint_rand_t randstate;
    fmpz_t subprod, cAksub, cBksub;
    int unlucky_count;
    fmpz_t Hmodulus;
    slong cur_zip_image, req_zip_images;
    ulong ABtotal_length;
    ulong Abidegree, Bbidegree;

    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);
    FLINT_ASSERT(bits == Gamma->bits);

    /* let's initialize everything at once to avoid complicated cleanup */
    Abidegree = _mpoly_bidegree(A->exps, bits, ctx->minfo);
    Bbidegree = _mpoly_bidegree(B->exps, bits, ctx->minfo);
    ABtotal_length = FLINT_MAX(A->length + B->length, 100);
    n_poly_init(Amarks);
    n_poly_init(Bmarks);
    n_poly_init(Hmarks);

    mpoly2_fill_marks(&Amarks->coeffs, &Amarks->length, &Amarks->alloc,
                                         A->exps, A->length, bits, ctx->minfo);

    mpoly2_fill_marks(&Bmarks->coeffs, &Bmarks->length, &Bmarks->alloc,
                                         B->exps, B->length, bits, ctx->minfo);


    flint_randinit(randstate);
    fmpz_init(p);
    fmpz_init(pm1); /* p - 1 */
    fmpz_init(subprod);
    fmpz_init(cAksub);
    fmpz_init(cBksub);
    fmpz_init(Hmodulus);
    fmpz_mpoly_init3(H, 0, bits, ctx);
    fmpz_mpoly_init3(Hcontent, 0, bits, ctx);

    mpoly_bma_interpolate_ctx_init(Ictx, ctx->minfo->nvars);

    /* multiprecision workspace */
    fmpz_init(image_count_mp);
    fmpz_init(sshift_mp);
    fmpz_init(last_unlucky_sshift_plus_1_mp);

    fmpz_set_ui(p, 2);    /* something positive */
    fmpz_mod_ctx_init_ui(ctx_mp, 2); /* modulus no care */

    fmpz_mod_poly_stack_init(St_mp->poly_stack);
    fmpz_mod_polyun_stack_init(St_mp->polyun_stack);

    fmpz_mod_bma_mpoly_init(GLambda_mp);
    fmpz_mod_bma_mpoly_init(AbarLambda_mp);
    fmpz_mod_bma_mpoly_init(BbarLambda_mp);

    fmpz_init(Gammaeval_mp);
    fmpz_mod_polyun_init(Aeval_mp, ctx_mp);
    fmpz_mod_polyun_init(Beval_mp, ctx_mp);
    fmpz_mod_polyun_init(Geval_mp, ctx_mp);
    fmpz_mod_polyun_init(Abareval_mp, ctx_mp);
    fmpz_mod_polyun_init(Bbareval_mp, ctx_mp);

    fmpz_mod_poly_init(Gammacur_mp, ctx_mp);
    fmpz_mod_poly_init(Gammainc_mp, ctx_mp);
    fmpz_mod_poly_init(Gammacoeff_mp, ctx_mp);
    fmpz_mod_polyun_init(Acur_mp, ctx_mp);
    fmpz_mod_polyun_init(Ainc_mp, ctx_mp);
    fmpz_mod_polyun_init(Acoeff_mp, ctx_mp);
    fmpz_mod_polyun_init(Bcur_mp, ctx_mp);
    fmpz_mod_polyun_init(Binc_mp, ctx_mp);
    fmpz_mod_polyun_init(Bcoeff_mp, ctx_mp);

    alphas_mp = _fmpz_vec_init(nvars);
    alpha_caches_mp = FLINT_ARRAY_ALLOC(nvars, fmpz_mod_poly_struct);
    for (i = 0; i < nvars; i++)
        fmpz_mod_poly_init(alpha_caches_mp + i, ctx_mp);

    /* machine precision workspace "sp" */
    nmod_init(&ctx_sp, 2);  /* modulus no care */

    n_poly_stack_init(St_sp->poly_stack);
    n_polyun_stack_init(St_sp->polyun_stack);

    nmod_bma_mpoly_init(GLambda_sp);
    nmod_bma_mpoly_init(AbarLambda_sp);
    nmod_bma_mpoly_init(BbarLambda_sp);

    n_polyun_init(Aeval_sp);
    n_polyun_init(Beval_sp);
    n_polyun_init(Geval_sp);
    n_polyun_init(Abareval_sp);
    n_polyun_init(Bbareval_sp);

    n_poly_init(Gammacur_sp);
    n_poly_init(Gammainc_sp);
    n_poly_init(Gammacoeff_sp);
    n_polyun_init(Acur_sp);
    n_polyun_init(Ainc_sp);
    n_polyun_init(Acoeff_sp);
    n_polyun_init(Bcur_sp);
    n_polyun_init(Binc_sp);
    n_polyun_init(Bcoeff_sp);

    alphas_sp = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    alpha_caches_sp = FLINT_ARRAY_ALLOC(3*nvars, n_poly_struct);
    for (i = 0; i < 3*nvars; i++)
        n_poly_init(alpha_caches_sp + i);

    n_polyun_init(HH);
    n_polyun_init(MH);
    n_polyun_init(ZH);
    n_poly_init(Hn);

    Adegs = FLINT_ARRAY_ALLOC(2*nvars, slong);
    Bdegs = Adegs + nvars;
    fmpz_mpoly_degrees_si(Adegs, A, ctx);
    fmpz_mpoly_degrees_si(Bdegs, B, ctx);

    GdegboundXY = FLINT_MIN(A->exps[0], B->exps[0]);
    p_sp = UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX);
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        p_sp = n_nextprime(p_sp, 1);
        nmod_init(&ctx_sp, p_sp);

        for (i = 2; i < ctx->minfo->nvars; i++)
        {
            alphas_sp[i] = n_urandint(randstate, p_sp - 1) + 1;
            nmod_pow_cache_start(alphas_sp[i], alpha_caches_sp + 3*i,
                         alpha_caches_sp + 3*i + 1, alpha_caches_sp + 3*i + 2);
        }

        fmpz_mpoly2_eval_nmod(Aeval_sp, Acur_sp, Ainc_sp, Acoeff_sp, A,
                 Amarks->coeffs, Amarks->length, alpha_caches_sp, ctx, ctx_sp);

        fmpz_mpoly2_eval_nmod(Beval_sp, Bcur_sp, Binc_sp, Bcoeff_sp, B,
                 Bmarks->coeffs, Bmarks->length, alpha_caches_sp, ctx, ctx_sp);

        if (Aeval_sp->length < 1 || Beval_sp->length < 1 ||
            n_polyu1n_bidegree(Aeval_sp) != Abidegree ||
            n_polyu1n_bidegree(Beval_sp) != Bbidegree)
        {
            continue;
        }
        success = n_polyu1n_mod_gcd_brown_smprime(Geval_sp, Abareval_sp,
                               Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp, St_sp);
        if (success)
        {
            FLINT_ASSERT(Geval_sp->length > 0);
            GdegboundXY = n_polyu1n_bidegree(Geval_sp);
            break;
        }
    }

    Ictx->degbounds[0] = 0;
    Ictx->degbounds[1] = 0;
    for (i = 2; i < nvars; i++)
    {
        Ictx->degbounds[i] = FLINT_MAX(Adegs[i], Bdegs[i]);
        Ictx->subdegs[i] = Ictx->degbounds[i]; /* will increment */
    }

    /* initialization done! */

    if (GdegboundXY == 0)
        goto gcd_is_trivial;

pick_ksub:

    fmpz_one(subprod);
    for (i = 2; i < nvars; i++)
    {
        if (n_add_checked(&Ictx->subdegs[i], Ictx->subdegs[i], 1))
        {
            success = 0;
            goto cleanup;
        }
        fmpz_mul_ui(subprod, subprod, Ictx->subdegs[i]);
    }

    _fmpz_mpoly_ksub_content(cAksub, A->coeffs, A->exps, Amarks->coeffs[1],
                                              bits, Ictx->subdegs, ctx->minfo);

    _fmpz_mpoly_ksub_content(cBksub, B->coeffs, B->exps, Bmarks->coeffs[1],
                                              bits, Ictx->subdegs, ctx->minfo);

    if (fmpz_is_zero(cAksub) || fmpz_is_zero(cBksub))
    {
        goto pick_ksub;
    }

pick_bma_prime:

    if (fmpz_cmp_ui(p, ABtotal_length) < 0)
        fmpz_set_ui(p, ABtotal_length);

    if (fmpz_cmp(p, subprod) < 0)
        fmpz_set(p, subprod);

    success = fmpz_next_smooth_prime(p, p);
    fmpz_sub_ui(pm1, p, 1);
    if (!success)
    {
        success = 0;
        goto cleanup;
    }

    if (fmpz_divisible(cAksub, p) || fmpz_divisible(cBksub, p))
    {
        goto pick_bma_prime;
    }

    for (i = 0; i < Gamma->length; i++)
        if (fmpz_divisible(Gamma->coeffs + i, p))
            goto pick_bma_prime;

    if (fmpz_abs_fits_ui(p))
    {
        p_sp = fmpz_get_ui(p);
        sshift_sp = 1;

        unlucky_count = 0;
        last_unlucky_sshift_plus_1_sp = 0;

        nmod_init(&ctx_sp, p_sp);
        nmod_discrete_log_pohlig_hellman_precompute_prime(Ictx->dlogenv_sp, p_sp);

        nmod_bma_mpoly_reset_prime(GLambda_sp, ctx_sp);
        nmod_bma_mpoly_reset_prime(AbarLambda_sp, ctx_sp);
        nmod_bma_mpoly_reset_prime(BbarLambda_sp, ctx_sp);
        nmod_bma_mpoly_zero(GLambda_sp);
        nmod_bma_mpoly_zero(AbarLambda_sp);
        nmod_bma_mpoly_zero(BbarLambda_sp);

        FLINT_ASSERT(sshift_sp == 1);
        nmod_mpoly_bma_interpolate_alpha_powers(alphas_sp, sshift_sp, 2,
                                                            Ictx, ctx, ctx_sp);
        for (j = 2; j < nvars; j++)
            nmod_pow_cache_start(alphas_sp[j], alpha_caches_sp + 3*j + 0,
                         alpha_caches_sp + 3*j + 1, alpha_caches_sp + 3*j + 2);

        mpoly2_nmod_monomial_evals(Ainc_sp, A->exps, bits, Amarks->coeffs,
                          Amarks->length, alpha_caches_sp, ctx->minfo, ctx_sp);

        mpoly2_nmod_monomial_evals(Binc_sp, B->exps, bits, Bmarks->coeffs,
                          Bmarks->length, alpha_caches_sp, ctx->minfo, ctx_sp);

        mpoly_nmod_monomial_evals(Gammainc_sp, Gamma->exps, Gamma->length, bits,
                          alpha_caches_sp + 3*2, 2, nvars, ctx->minfo, ctx_sp);

        fmpz_mpoly2_nmod_coeffs(Acoeff_sp, A->coeffs, Amarks->coeffs,
                                                       Amarks->length, ctx_sp);

        fmpz_mpoly2_nmod_coeffs(Bcoeff_sp, B->coeffs, Bmarks->coeffs,
                                                       Bmarks->length, ctx_sp);

        fmpz_mpoly_nmod_coeffs(Gammacoeff_sp, Gamma->coeffs, Gamma->length, ctx_sp);

        n_polyun_set(Acur_sp, Ainc_sp);
        n_polyun_set(Bcur_sp, Binc_sp);
        n_poly_set(Gammacur_sp, Gammainc_sp);

        image_count_sp = 0;

    next_bma_image_sp:

        /* image count is also the current power of alpha we are evaluating */
        image_count_sp++;

        FLINT_ASSERT(GLambda_sp->pointcount == AbarLambda_sp->pointcount);
        FLINT_ASSERT(GLambda_sp->pointcount == BbarLambda_sp->pointcount);
        FLINT_ASSERT(sshift_sp + GLambda_sp->pointcount == image_count_sp);

        if (image_count_sp >= p_sp - 1)
            goto pick_bma_prime;

        n_polyun_mod_zip_eval_cur_inc_coeff(Aeval_sp, Acur_sp, Ainc_sp, Acoeff_sp, ctx_sp);
        n_polyun_mod_zip_eval_cur_inc_coeff(Beval_sp, Bcur_sp, Binc_sp, Bcoeff_sp, ctx_sp);
        Gammaeval_sp = n_poly_mod_zip_eval_cur_inc_coeff(Gammacur_sp, Gammainc_sp, Gammacoeff_sp, ctx_sp);

        if (Aeval_sp->length < 1 || Beval_sp->length < 1 ||
            n_polyu1n_bidegree(Aeval_sp) != Abidegree ||
            n_polyu1n_bidegree(Beval_sp) != Bbidegree)
        {
            sshift_sp += GLambda_sp->pointcount + 1;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;
        }

        /* the evaluation killed neither lc(A) nor lc(B) */
        FLINT_ASSERT(Gammaeval_sp != 0);

        success = n_polyu1n_mod_gcd_brown_smprime(Geval_sp, Abareval_sp,
                             Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp, St_sp);
        if (!success)
        {
            sshift_sp += GLambda_sp->pointcount + 1;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;
        }

        FLINT_ASSERT(Geval_sp->length > 0);
        GevaldegXY = n_polyu1n_bidegree(Geval_sp);
        _n_poly_vec_mul_nmod_intertible(Geval_sp->coeffs, Geval_sp->length,
                                                         Gammaeval_sp, ctx_sp);

        if (GdegboundXY < GevaldegXY)
        {
            if (sshift_sp == last_unlucky_sshift_plus_1_sp)
                goto pick_ksub;

            if (++unlucky_count > 2)
                goto pick_bma_prime;

            last_unlucky_sshift_plus_1_sp = sshift_sp + 1;
            sshift_sp += GLambda_sp->pointcount + 1;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;        
        }
        else if (GdegboundXY > GevaldegXY)
        {
            /* new bound on deg_XY(G) */
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            sshift_sp += GLambda_sp->pointcount;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            nmod_bma_mpoly_add_point(GLambda_sp, Geval_sp, ctx_sp);
            nmod_bma_mpoly_add_point(AbarLambda_sp, Abareval_sp, ctx_sp);
            nmod_bma_mpoly_add_point(BbarLambda_sp, Bbareval_sp, ctx_sp);
            goto next_bma_image_sp;
        }

        nmod_bma_mpoly_add_point(GLambda_sp, Geval_sp, ctx_sp);
        nmod_bma_mpoly_add_point(AbarLambda_sp, Abareval_sp, ctx_sp);
        nmod_bma_mpoly_add_point(BbarLambda_sp, Bbareval_sp, ctx_sp);

        if ((GLambda_sp->pointcount & 7) != 0)
        {
            goto next_bma_image_sp;
        }

        if (GLambda_sp->pointcount/2 >= Gamma->length &&
            !nmod_bma_mpoly_reduce(GLambda_sp) &&
            nmod_bma_mpoly_get_fmpz_mpoly2(H, Hmarks, ctx, sshift_sp, GLambda_sp, Ictx, ctx_sp) &&
            Hmarks->coeffs[1] == Gamma->length)
        {
            which_check = 0;
            goto check_sp;
        }

        if (AbarLambda_sp->pointcount/2 >= Amarks->coeffs[1] &&
            !nmod_bma_mpoly_reduce(AbarLambda_sp) &&
            nmod_bma_mpoly_get_fmpz_mpoly2(H, Hmarks, ctx, sshift_sp, AbarLambda_sp, Ictx, ctx_sp) &&
            Hmarks->coeffs[1] == Amarks->coeffs[1])
        {
            which_check = 1;
            goto check_sp;
        }

        if (BbarLambda_sp->pointcount/2 >= Bmarks->coeffs[1] &&
            !nmod_bma_mpoly_reduce(BbarLambda_sp) &&
            nmod_bma_mpoly_get_fmpz_mpoly2(H, Hmarks, ctx, sshift_sp, BbarLambda_sp, Ictx, ctx_sp) &&
            Hmarks->coeffs[1] == Bmarks->coeffs[1])
        {
            which_check = 2;
            goto check_sp;
        }

        if (GLambda_sp->pointcount/2 > ABtotal_length)
        {
            success = 0;
            goto cleanup;
        }

        goto next_bma_image_sp;

    check_sp:

        FLINT_ASSERT(0 <= which_check && which_check <= 2);
        FLINT_ASSERT(which_check != 0 || GdegboundXY == _mpoly_bidegree(H->exps, bits, ctx->minfo));

        success =  _random_check_sp(
                    &GevaldegXY,
                    GdegboundXY,
                    which_check,
                    Aeval_sp, Acur_sp, Ainc_sp, Acoeff_sp,
                    Beval_sp, Bcur_sp, Binc_sp, Bcoeff_sp,
                    Geval_sp,
                    Abareval_sp,
                    Bbareval_sp,
                    alphas_sp,
                    alpha_caches_sp,
                    H, Hmarks,
                    A, Amarks, Abidegree,
                    B, Bmarks, Bbidegree,
                    Gamma,
                    ctx,
                    ctx_sp,
                    randstate,
                    St_sp);

        if (success < 1)
        {
            if (success == 0)
                goto next_bma_image_sp;

            FLINT_ASSERT(GdegboundXY > GevaldegXY);
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            sshift_sp += GLambda_sp->pointcount;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;
        }
    }
    else
    {
        fmpz_one(sshift_mp);

        unlucky_count = 0;
        fmpz_zero(last_unlucky_sshift_plus_1_mp);

        fmpz_mod_ctx_set_modulus(ctx_mp, p);
        fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(Ictx->dlogenv, p);
        fmpz_mod_bma_mpoly_zero(GLambda_mp);
        fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
        fmpz_mod_bma_mpoly_zero(BbarLambda_mp);

        FLINT_ASSERT(fmpz_is_one(sshift_mp));
        fmpz_mod_mpoly_bma_interpolate_alpha_powers(alphas_mp, sshift_mp, 2,
                                                            Ictx, ctx, ctx_mp); 
        for (j = 2; j < nvars; j++)
            fmpz_mod_pow_cache_start(alphas_mp + j, alpha_caches_mp + j, ctx_mp);

        mpoly2_monomial_evals_fmpz_mod(Ainc_mp, A->exps, bits, Amarks->coeffs,
               Amarks->length, alpha_caches_mp + 2, nvars, ctx->minfo, ctx_mp);

        mpoly2_monomial_evals_fmpz_mod(Binc_mp, B->exps, bits, Bmarks->coeffs,
               Bmarks->length, alpha_caches_mp + 2, nvars, ctx->minfo, ctx_mp);

        mpoly_monomial_evals_fmpz_mod(Gammainc_mp, Gamma->exps, Gamma->length,
                      bits, alpha_caches_mp + 2, 2, nvars, ctx->minfo, ctx_mp);

        fmpz_mpoly2_fmpz_mod_coeffs(Acoeff_mp, A->coeffs, Amarks->coeffs,
                                                       Amarks->length, ctx_mp);

        fmpz_mpoly2_fmpz_mod_coeffs(Bcoeff_mp, B->coeffs, Bmarks->coeffs,
                                                       Bmarks->length, ctx_mp);

        fmpz_mpoly_fmpz_mod_coeffs(Gammacoeff_mp, Gamma->coeffs, Gamma->length,
                                                                       ctx_mp);

        fmpz_mod_polyun_set(Acur_mp, Ainc_mp, ctx_mp);
        fmpz_mod_polyun_set(Bcur_mp, Binc_mp, ctx_mp);
        fmpz_mod_poly_set(Gammacur_mp, Gammainc_mp, ctx_mp);

        fmpz_zero(image_count_mp);

    next_bma_image_mp:

        /* image count is also the current power of alpha we are evaluating */
        fmpz_add_ui(image_count_mp, image_count_mp, 1);

        FLINT_ASSERT(GLambda_sp->pointcount == AbarLambda_sp->pointcount);
        FLINT_ASSERT(GLambda_sp->pointcount == BbarLambda_sp->pointcount);
    #if FLINT_WANT_ASSERT
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_add_ui(t, sshift_mp, GLambda_mp->pointcount);
            FLINT_ASSERT(fmpz_equal(t, image_count_mp));
            fmpz_clear(t);
        }
    #endif

        if (fmpz_cmp(image_count_mp, pm1) >= 0)
            goto pick_bma_prime;

        fmpz_mod_polyu2n_zip_eval_cur_inc_coeff(Aeval_mp, Acur_mp, Ainc_mp, Acoeff_mp, ctx_mp);
        fmpz_mod_polyu2n_zip_eval_cur_inc_coeff(Beval_mp, Bcur_mp, Binc_mp, Bcoeff_mp, ctx_mp);
        fmpz_mod_poly_zip_eval_cur_inc_coeff(Gammaeval_mp, Gammacur_mp, Gammainc_mp, Gammacoeff_mp, ctx_mp);

        if (Aeval_mp->length < 1 || Beval_mp->length < 1 ||
            fmpz_mod_polyu1n_bidegree(Aeval_mp) != Abidegree ||
            fmpz_mod_polyu1n_bidegree(Beval_mp) != Bbidegree)
        {
            fmpz_add_ui(sshift_mp, sshift_mp, GLambda_mp->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(GLambda_mp);
            fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(BbarLambda_mp);
            goto next_bma_image_mp;
        }

        /* the evaluation killed neither lc(A) nor lc(B) */
        FLINT_ASSERT(!fmpz_is_zero(Gammaeval_mp));

        success = fmpz_mod_polyu1n_gcd_brown_smprime(Geval_mp,
                         Abareval_mp, Bbareval_mp, Aeval_mp, Beval_mp, ctx_mp,
                                       St_mp->poly_stack, St_mp->polyun_stack);
        if (!success)
        {
            fmpz_add_ui(sshift_mp, sshift_mp, GLambda_mp->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(GLambda_mp);
            fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(BbarLambda_mp);
            goto next_bma_image_mp;
        }

        FLINT_ASSERT(Geval_mp->length > 0);
        GevaldegXY = fmpz_mod_polyu1n_bidegree(Geval_mp);
        _fmpz_mod_poly_vec_mul_fmpz_mod(Geval_mp->coeffs, Geval_mp->length,
                                                         Gammaeval_mp, ctx_mp);

        FLINT_ASSERT(fmpz_equal(Gammaeval_mp, fmpz_mod_polyun_leadcoeff(Geval_mp)));

        if (GdegboundXY < GevaldegXY)
        {
            if (fmpz_equal(sshift_mp, last_unlucky_sshift_plus_1_mp))
                goto pick_ksub;

            if (++unlucky_count > 2)
                goto pick_bma_prime;

            fmpz_add_ui(last_unlucky_sshift_plus_1_mp, sshift_mp, 1);
            fmpz_add_ui(sshift_mp, sshift_mp, GLambda_mp->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(GLambda_mp);
            fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(BbarLambda_mp);
            goto next_bma_image_mp;
        }
        else if (GdegboundXY > GevaldegXY)
        {
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            fmpz_add_ui(sshift_mp, sshift_mp, GLambda_mp->pointcount);
            fmpz_mod_bma_mpoly_zero(GLambda_mp);
            fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(BbarLambda_mp);
            fmpz_mod_bma_mpoly_add_point(GLambda_mp, Geval_mp, ctx_mp);
            fmpz_mod_bma_mpoly_add_point(AbarLambda_mp, Abareval_mp, ctx_mp);
            fmpz_mod_bma_mpoly_add_point(BbarLambda_mp, Bbareval_mp, ctx_mp);
            goto next_bma_image_mp;
        }

        fmpz_mod_bma_mpoly_add_point(GLambda_mp, Geval_mp, ctx_mp);
        fmpz_mod_bma_mpoly_add_point(AbarLambda_mp, Abareval_mp, ctx_mp);
        fmpz_mod_bma_mpoly_add_point(BbarLambda_mp, Bbareval_mp, ctx_mp);

        if ((GLambda_mp->pointcount & 7) != 0)
        {
            goto next_bma_image_mp;
        }

        if (GLambda_mp->pointcount/2 >= Gamma->length &&
            !fmpz_mod_bma_mpoly_reduce(GLambda_mp, ctx_mp) &&
            fmpz_mod_bma_mpoly_get_fmpz_mpoly2(H, Hmarks, ctx, sshift_mp, GLambda_mp, Ictx, ctx_mp) &&
            Hmarks->coeffs[1] == Gamma->length)
        {
            which_check = 0;
            goto check_mp;
        }

        if (AbarLambda_mp->pointcount/2 >= Bmarks->coeffs[1] &&
            !fmpz_mod_bma_mpoly_reduce(AbarLambda_mp, ctx_mp) &&
            fmpz_mod_bma_mpoly_get_fmpz_mpoly2(H, Hmarks, ctx, sshift_mp, AbarLambda_mp, Ictx, ctx_mp) &&
            Hmarks->coeffs[1] == Amarks->coeffs[1])
        {
            which_check = 1;
            goto check_mp;
        }

        if (BbarLambda_mp->pointcount/2 >= Bmarks->coeffs[1] &&
            !fmpz_mod_bma_mpoly_reduce(BbarLambda_mp, ctx_mp) &&
            fmpz_mod_bma_mpoly_get_fmpz_mpoly2(H, Hmarks, ctx, sshift_mp, BbarLambda_mp, Ictx, ctx_mp) &&
            Hmarks->coeffs[1] == Bmarks->coeffs[1])
        {
            which_check = 2;
            goto check_mp;
        }

        if (GLambda_mp->pointcount/2 > ABtotal_length)
        {
            success = 0;
            goto cleanup;
        }

        goto next_bma_image_mp;

    check_mp:

        FLINT_ASSERT(0 <= which_check && which_check <= 2);
        FLINT_ASSERT(which_check != 0 || GdegboundXY == _mpoly_bidegree(H->exps, bits, ctx->minfo));

        success = _random_check_mp(
                    &GevaldegXY,
                    GdegboundXY,
                    which_check,
                    Aeval_mp, Acur_mp, Ainc_mp, Acoeff_mp,
                    Beval_mp, Bcur_mp, Binc_mp, Bcoeff_mp,
                    Geval_mp,
                    Abareval_mp,
                    Bbareval_mp,
                    Gammaeval_mp,
                    alphas_mp,
                    alpha_caches_mp,
                    H, Hmarks,
                    A, Amarks, Abidegree,
                    B, Bmarks, Bbidegree,
                    Gamma,
                    ctx,
                    ctx_mp,
                    randstate,
                    St_mp);

        if (success < 1)
        {
            if (success == 0)
                goto next_bma_image_mp;

            FLINT_ASSERT(GdegboundXY > GevaldegXY);
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            fmpz_add_ui(sshift_mp, sshift_mp, GLambda_mp->pointcount);
            fmpz_mod_bma_mpoly_zero(GLambda_mp);
            fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(BbarLambda_mp);
            goto next_bma_image_mp;
        }
    }

    /* assume that H is correct modulo Hmodulus = p */
    FLINT_ASSERT(0 <= which_check && which_check <= 2);
    fmpz_set(Hmodulus, p);

    p_sp = UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX);

pick_zip_prime:

    if (p_sp >= UWORD_MAX_PRIME)
    {
        success = 0;
        goto cleanup;
    }
    p_sp = n_nextprime(p_sp, 1);

    if (0 == fmpz_fdiv_ui(Hmodulus, p_sp))
        goto pick_zip_prime;

    nmod_init(&ctx_sp, p_sp);

    FLINT_ASSERT(p_sp > 3);
    for (i = 2; i < ctx->minfo->nvars; i++)
    {
        alphas_sp[i] = n_urandint(randstate, p_sp - 3) + 2;
        nmod_pow_cache_start(alphas_sp[i], alpha_caches_sp + 3*i + 0,
                         alpha_caches_sp + 3*i + 1, alpha_caches_sp + 3*i + 2);
    }

    mpoly2_nmod_monomial_evals(HH, H->exps, bits, Hmarks->coeffs,
                    Hmarks->length, alpha_caches_sp, ctx->minfo, ctx_sp);

    req_zip_images = n_polyun_product_roots(MH, HH, ctx_sp);
    req_zip_images++;

    mpoly2_nmod_monomial_evals(Ainc_sp, A->exps, bits, Amarks->coeffs,
                      Amarks->length, alpha_caches_sp, ctx->minfo, ctx_sp);

    mpoly2_nmod_monomial_evals(Binc_sp, B->exps, bits, Bmarks->coeffs,
                      Bmarks->length, alpha_caches_sp, ctx->minfo, ctx_sp);

    mpoly_nmod_monomial_evals(Gammainc_sp, Gamma->exps, Gamma->length, bits,
                      alpha_caches_sp + 3*2, 2, nvars, ctx->minfo, ctx_sp);

    fmpz_mpoly2_nmod_coeffs(Acoeff_sp, A->coeffs, Amarks->coeffs,
                                                   Amarks->length, ctx_sp);

    fmpz_mpoly2_nmod_coeffs(Bcoeff_sp, B->coeffs, Bmarks->coeffs,
                                                   Bmarks->length, ctx_sp);

    fmpz_mpoly_nmod_coeffs(Gammacoeff_sp, Gamma->coeffs, Gamma->length, ctx_sp);

    n_polyun_set(Acur_sp, Ainc_sp);
    n_polyun_set(Bcur_sp, Binc_sp);
    n_poly_set(Gammacur_sp, Gammainc_sp);

    n_polyun_zip_start(ZH, HH, req_zip_images);

    for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
    {
        n_polyun_mod_zip_eval_cur_inc_coeff(Aeval_sp, Acur_sp, Ainc_sp, Acoeff_sp, ctx_sp);
        n_polyun_mod_zip_eval_cur_inc_coeff(Beval_sp, Bcur_sp, Binc_sp, Bcoeff_sp, ctx_sp);
        Gammaeval_sp = n_poly_mod_zip_eval_cur_inc_coeff(Gammacur_sp, Gammainc_sp, Gammacoeff_sp, ctx_sp);

        if (Aeval_sp->length < 1 || Beval_sp->length < 1 ||
            n_polyu1n_bidegree(Aeval_sp) != Abidegree ||
            n_polyu1n_bidegree(Beval_sp) != Bbidegree)
        {
            goto pick_zip_prime;
        }

        FLINT_ASSERT(Gammaeval_sp != 0);

        success = n_polyu1n_mod_gcd_brown_smprime(Geval_sp, Abareval_sp,
                               Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp, St_sp);
        if (!success)
            goto pick_zip_prime;

        FLINT_ASSERT(Geval_sp->length > 0);
        GevaldegXY = n_polyu1n_bidegree(Geval_sp);

        if (GevaldegXY > GdegboundXY)
            goto pick_zip_prime;

        if (GevaldegXY < GdegboundXY)
        {
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            goto pick_bma_prime;
        }

        _n_poly_vec_mul_nmod_intertible(Geval_sp->coeffs, Geval_sp->length,
                                                         Gammaeval_sp, ctx_sp);

        success = n_polyu2n_add_zipun_must_match(ZH,
                                which_check == 1 ? Abareval_sp :
                                which_check == 2 ? Bbareval_sp : Geval_sp,
                                                                cur_zip_image);
        if (!success)
            goto pick_bma_prime;
    }

    FLINT_ASSERT(H->length == Hmarks->coeffs[Hmarks->length]);
    n_poly_fit_length(Hn, H->length);

    success = zip_solve(Hn->coeffs, ZH, HH, MH, ctx_sp);
    if (success < 0)
        goto pick_zip_prime; /* singular */
    if (success == 0)
        goto pick_bma_prime; /* no match */

    changed = _fmpz_vec_crt_nmod(&Hbits, H->coeffs, Hmodulus,
                                                Hn->coeffs, H->length, ctx_sp);
    fmpz_mul_ui(Hmodulus, Hmodulus, ctx_sp.n);
    if (changed)
    {
        if (Hbits > Hbitbound)
            goto pick_bma_prime;
        goto pick_zip_prime;
    }

    success = fmpz_mpolyl_content(Hcontent, H, 2, ctx);
    if (!success)
        goto cleanup;

    if (which_check == 1)
    {
        success = fmpz_mpoly_divides(Abar, H, Hcontent, ctx);
        FLINT_ASSERT(success);
        if (!fmpz_mpoly_divides(G, A, Abar, ctx) ||
            !fmpz_mpoly_divides(Bbar, B, G, ctx))
        {
            goto pick_zip_prime;
        }
    }
    else if (which_check == 2)
    {
        success = fmpz_mpoly_divides(Bbar, H, Hcontent, ctx);
        FLINT_ASSERT(success);
        if (!fmpz_mpoly_divides(G, B, Bbar, ctx) ||
            !fmpz_mpoly_divides(Abar, A, G, ctx))
        {
            goto pick_zip_prime;
        }
    }
    else
    {
        FLINT_ASSERT(which_check == 0);
        success = fmpz_mpoly_divides(G, H, Hcontent, ctx);
        FLINT_ASSERT(success);
        if (!fmpz_mpoly_divides(Abar, A, G, ctx) ||
            !fmpz_mpoly_divides(Bbar, B, G, ctx))
        {
            goto pick_zip_prime;
        }
    }

    success = 1;

cleanup:

    n_poly_clear(Amarks);
    n_poly_clear(Bmarks);
    n_poly_clear(Hmarks);

    n_polyun_clear(HH);
    n_polyun_clear(MH);
    n_polyun_clear(ZH);
    n_poly_clear(Hn);

    /* machine precision workspace */
    flint_free(alphas_sp);

    for (i = 0; i < 3*nvars; i++)
        n_poly_clear(alpha_caches_sp + i);
    flint_free(alpha_caches_sp);

    n_polyun_clear(Aeval_sp);
    n_polyun_clear(Beval_sp);
    n_polyun_clear(Geval_sp);
    n_polyun_clear(Abareval_sp);
    n_polyun_clear(Bbareval_sp);

    n_poly_clear(Gammacur_sp);
    n_polyun_clear(Acur_sp);
    n_polyun_clear(Bcur_sp);
    n_poly_clear(Gammainc_sp);
    n_polyun_clear(Ainc_sp);
    n_polyun_clear(Binc_sp);
    n_poly_clear(Gammacoeff_sp);
    n_polyun_clear(Acoeff_sp);
    n_polyun_clear(Bcoeff_sp);

    nmod_bma_mpoly_clear(GLambda_sp);
    nmod_bma_mpoly_clear(AbarLambda_sp);
    nmod_bma_mpoly_clear(BbarLambda_sp);
    n_poly_stack_clear(St_sp->poly_stack);
    n_polyun_stack_clear(St_sp->polyun_stack);

    /* multiprecision workspace */
    _fmpz_vec_clear(alphas_mp, ctx->minfo->nvars);

    for (i = 0; i < nvars; i++)
        fmpz_mod_poly_clear(alpha_caches_mp + i, ctx_mp);
    flint_free(alpha_caches_mp);

    fmpz_clear(Gammaeval_mp);
    fmpz_mod_polyun_clear(Aeval_mp, ctx_mp);
    fmpz_mod_polyun_clear(Beval_mp, ctx_mp);
    fmpz_mod_polyun_clear(Geval_mp, ctx_mp);
    fmpz_mod_polyun_clear(Abareval_mp, ctx_mp);
    fmpz_mod_polyun_clear(Bbareval_mp, ctx_mp);

    fmpz_mod_poly_clear(Gammacur_mp, ctx_mp);
    fmpz_mod_polyun_clear(Acur_mp, ctx_mp);
    fmpz_mod_polyun_clear(Bcur_mp, ctx_mp);
    fmpz_mod_poly_clear(Gammainc_mp, ctx_mp);
    fmpz_mod_polyun_clear(Ainc_mp, ctx_mp);
    fmpz_mod_polyun_clear(Binc_mp, ctx_mp);
    fmpz_mod_poly_clear(Gammacoeff_mp, ctx_mp);
    fmpz_mod_polyun_clear(Acoeff_mp, ctx_mp);
    fmpz_mod_polyun_clear(Bcoeff_mp, ctx_mp);
    fmpz_mod_poly_stack_clear(St_mp->poly_stack);
    fmpz_mod_polyun_stack_clear(St_mp->polyun_stack);

    fmpz_mod_bma_mpoly_clear(GLambda_mp, ctx_mp);
    fmpz_mod_bma_mpoly_clear(AbarLambda_mp, ctx_mp);
    fmpz_mod_bma_mpoly_clear(BbarLambda_mp, ctx_mp);
    fmpz_mod_ctx_clear(ctx_mp);
    fmpz_clear(last_unlucky_sshift_plus_1_mp);
    fmpz_clear(sshift_mp);
    fmpz_clear(image_count_mp);

    fmpz_mpoly_clear(Hcontent, ctx);
    fmpz_mpoly_clear(H, ctx);
    flint_free(Adegs);

    mpoly_bma_interpolate_ctx_clear(Ictx);

    fmpz_clear(Hmodulus);
    fmpz_clear(cBksub);
    fmpz_clear(cAksub);
    fmpz_clear(subprod);
    fmpz_clear(pm1);
    fmpz_clear(p);
    flint_randclear(randstate);

    FLINT_ASSERT(G->bits == bits);
    FLINT_ASSERT(Abar->bits == bits);
    FLINT_ASSERT(Bbar->bits == bits);

    return success;

gcd_is_trivial:

    fmpz_mpoly_one(G, ctx);
    fmpz_mpoly_set(Abar, A, ctx);
    fmpz_mpoly_set(Bbar, B, ctx);
    success = 1;
    goto cleanup;
}

