/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

/*
    evaluation at
        gen(start)   -> alpha[0]
        gen(start+1) -> alpha[1]
        ...
        gen(stop-1)  -> alpha[stop-start-1]

    the other gen are assumed to not appear in A
    alpha[i] is represented by alpha_caches[3*i], ..., alpha_caches[3*i+2]
*/
void mpoly_monomial_evals_nmod(
    n_poly_t EH,
    const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
    n_poly_struct * alpha_caches,
    slong start,
    slong stop,
    const mpoly_ctx_t mctx,
    const nmod_t fpctx)
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
            p[i] = nmod_pow_cache_mulpow_ui(p[i], ei, alpha_caches + 3*k + 0, 
                                                alpha_caches + 3*k + 1,
                                                alpha_caches + 3*k + 2, fpctx);
        }
    }

    TMP_END;
}

/*
    evaluation at

    gen(0) -> x
    gen(1) -> alpha[0]
    gen(2) -> alpha[1]
    gen(3) -> alpha[2]
    ...
    gen(m-1) -> alpha[m-2]

    the univariate marks should be filled in by mpoly1_fill_marks
    alpha[i] is represented by alpha_caches[3*i], ..., alpha_caches[3*i+2]
*/
void mpoly1_monomial_evals_nmod(
    n_polyun_t EH,
    const ulong * Aexps, flint_bitcnt_t Abits, const ulong * Amarks, slong Amarkslen,
    n_poly_struct * alpha_caches,
    slong m,
    const mpoly_ctx_t mctx,
    const nmod_t fpctx)
{
    slong start, stop, i, j, k, n;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    mp_limb_t * p;
    TMP_INIT;

    FLINT_ASSERT(1 < m && m <= mctx->nvars);
    FLINT_ASSERT(Amarkslen > 0);

    TMP_START;

    off = TMP_ARRAY_ALLOC(2*m, slong);
    shift = off + m;
    for (k = 0; k < m; k++)
        mpoly_gen_offset_shift_sp(&off[k], &shift[k], k, Abits, mctx);

    n_polyun_fit_length(EH, Amarkslen);

    for (i = 0; i < Amarkslen; i++)
    {
        start = Amarks[i];
        stop = Amarks[i + 1];
        FLINT_ASSERT(start < stop);
        n = stop - start;

        EH->exps[i] = (Aexps[N*start + off[0]] >> shift[0]) & mask;
        n_poly_fit_length(EH->coeffs + i, n);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            p[j] = 1;
            for (k = 1; k < m; k++)
            {
                ulong ei = (Aexps[N*(start + j) + off[k]] >> shift[k]) & mask;
                p[j] = nmod_pow_cache_mulpow_ui(p[j], ei,
                                          alpha_caches + 3*(k - 1) + 0,
                                          alpha_caches + 3*(k - 1) + 1,
                                          alpha_caches + 3*(k - 1) + 2, fpctx);
            }
        }
    }

    EH->length = Amarkslen;

    TMP_END;
}

/*
    evaluation at

    gen(0) -> x
    gen(1) -> y
    gen(2) -> alpha[0]
    gen(3) -> alpha[1]
    ...
    gen(m-1) -> alpha[m-3]

    the bivariate marks should be filled in by mpoly2_fill_marks
    alpha[i] is represented by alpha_caches[3*i], ..., alpha_caches[3*i+2]
*/
void mpoly2_monomial_evals_nmod(
    n_polyun_t EH,
    const ulong * Aexps, flint_bitcnt_t Abits, ulong * Amarks, slong Amarkslen,
    n_poly_struct * alpha_caches,
    slong m,
    const mpoly_ctx_t mctx,
    const nmod_t fpctx)
{
    slong start, stop, i, j, k, n;
    ulong e0, e1;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    mp_limb_t * p;
    TMP_INIT;

    FLINT_ASSERT(2 < m && m <= mctx->nvars);
    FLINT_ASSERT(Amarkslen > 0);

    TMP_START;

    off = TMP_ARRAY_ALLOC(2*m, slong);
    shift = off + m;
    for (k = 0; k < m; k++)
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
            for (k = 2; k < m; k++)
            {
                ulong ei = (Aexps[N*(start + j) + off[k]] >> shift[k]) & mask;
                p[j] = nmod_pow_cache_mulpow_ui(p[j], ei,
                                          alpha_caches + 3*(k - 2) + 0,
                                          alpha_caches + 3*(k - 2) + 1,
                                          alpha_caches + 3*(k - 2) + 2, fpctx);
            }
        }
    }

    EH->length = Amarkslen;

    TMP_END;
}

void n_polyun_zip_start(n_polyun_t Z, n_polyun_t H, slong req_images)
{
    slong j;
    n_polyun_fit_length(Z, H->length);
    Z->length = H->length;
    for (j = 0; j < H->length; j++)
    {
        Z->exps[j] = H->exps[j];
        n_poly_fit_length(Z->coeffs + j, req_images);
        Z->coeffs[j].length = 0;
    }
}


int n_polyu2n_add_zip_must_match(
    n_polyun_t Z,
    const n_bpoly_t A,
    slong cur_length)
{
    slong i, Ai, ai;
    const n_poly_struct * Acoeffs = A->coeffs;

    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(A->coeffs + Ai);

    for (i = 0; i < Z->length; i++)
    {
        if (Ai >= 0 && Z->exps[i] == pack_exp2(Ai, ai))
        {
            /* Z present, A present */
            Z->coeffs[i].coeffs[cur_length] = Acoeffs[Ai].coeffs[ai];
            Z->coeffs[i].length = cur_length + 1;
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
        else if (Ai < 0 || Z->exps[i] > pack_exp2(Ai, ai))
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

int n_polyun_zip_solve(
    nmod_mpoly_t A,
    n_polyun_t Z,
    n_polyun_t H,
    n_polyun_t M,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong Ai, i, n;
    mp_limb_t * Acoeffs = A->coeffs;
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
        FLINT_ASSERT(Ai + n <= A->length);

        n_poly_fit_length(t, n);

        success = _nmod_zip_vand_solve(Acoeffs + Ai,
                                 H->coeffs[i].coeffs, n,
                                 Z->coeffs[i].coeffs, Z->coeffs[i].length,
                                 M->coeffs[i].coeffs, t->coeffs, ctx->mod);
        if (success < 1)
        {
            n_poly_clear(t);
            return success;
        }

        Ai += n;
        FLINT_ASSERT(Ai <= A->length);

    }

    FLINT_ASSERT(Ai == A->length);

    n_poly_clear(t);
    return 1;
}

