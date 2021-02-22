/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"
#include "fmpz_mod_mpoly_factor.h"

/*
    evaluation at
        gen(start) -> caches[0]
        gen(start+1) -> caches[1]
        ...
        gen(stop-1) -> caches[stop-start-1]

    the other gen are assumed to not appear in A
*/
void mpoly_monomial_evals_fmpz_mod(
    fmpz_mod_poly_t EH,
    const ulong * Aexps, slong Alen, flint_bitcnt_t Abits,
    fmpz_mod_poly_struct * alpha_caches,
    slong start,
    slong stop,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, k;
    fmpz * p;
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

    fmpz_mod_poly_fit_length(EH, Alen, fpctx);
    EH->length = Alen;
    p = EH->coeffs;

    for (i = 0; i < Alen; i++)
    {
        fmpz_one(p + i);
        for (k = 0; k < num; k++)
        {
            ulong ei = (Aexps[N*i + off[k]] >> shift[k]) & mask;
            fmpz_mod_pow_cache_mulpow_ui(p + i, p + i, ei, alpha_caches + k, fpctx);
        }
    }

    TMP_END;
}

/*
    evaluation at

    gen(0) -> x
    gen(1) -> alphas[0]
    gen(2) -> alphas[1]
    gen(3) -> alphas[2]
    ...
    gen(m-1) -> alphas[m-2]

    the uniivariate marks should be filled in by mpoly1_fill_marks
*/
void mpoly1_monomial_evals_fmpz_mod(
    fmpz_mod_polyun_t EH,
    const ulong * Aexps, flint_bitcnt_t Abits, const ulong * Amarks, slong Amarkslen,
    fmpz_mod_poly_struct * alpha_caches,
    slong m,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong start, stop, i, j, k, n;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    fmpz * p;
    TMP_INIT;

    FLINT_ASSERT(1 < m && m <= mctx->nvars);
    FLINT_ASSERT(Amarkslen > 0);

    TMP_START;

    off = TMP_ARRAY_ALLOC(2*m, slong);
    shift = off + m;
    for (k = 0; k < m; k++)
        mpoly_gen_offset_shift_sp(&off[k], &shift[k], k, Abits, mctx);

    fmpz_mod_polyun_fit_length(EH, Amarkslen, fpctx);

    for (i = 0; i < Amarkslen; i++)
    {
        start = Amarks[i];
        stop = Amarks[i + 1];
        FLINT_ASSERT(start < stop);
        n = stop - start;

        EH->exps[i] = (Aexps[N*start + off[0]] >> shift[0]) & mask;
        fmpz_mod_poly_fit_length(EH->coeffs + i, n, fpctx);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            fmpz_one(p + j);
            for (k = 1; k < m; k++)
            {
                ulong ei = (Aexps[N*(start + j) + off[k]] >> shift[k]) & mask;
                fmpz_mod_pow_cache_mulpow_ui(p + j, p + j, ei,
                                                  alpha_caches + k - 1, fpctx);
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
    gen(2) -> alphas[0]
    gen(3) -> alphas[1]
    ...
    gen(m-1) -> alphas[m-3]

    the bivariate marks should be filled in by mpoly2_fill_marks
*/
void mpoly2_monomial_evals_fmpz_mod(
    fmpz_mod_polyun_t EH,
    const ulong * Aexps, flint_bitcnt_t Abits, ulong * Amarks, slong Amarkslen,
    fmpz_mod_poly_struct * alpha_caches,
    slong m,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong start, stop, i, j, k, n;
    ulong e0, e1;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    fmpz * p;
    TMP_INIT;

    FLINT_ASSERT(2 < m && m <= mctx->nvars);
    FLINT_ASSERT(Amarkslen > 0);

    TMP_START;

    off = TMP_ARRAY_ALLOC(2*m, slong);
    shift = off + m;
    for (k = 0; k < m; k++)
        mpoly_gen_offset_shift_sp(&off[k], &shift[k], k, Abits, mctx);

    fmpz_mod_polyun_fit_length(EH, Amarkslen, fpctx);

    for (i = 0; i < Amarkslen; i++)
    {
        start = Amarks[i];
        stop = Amarks[i + 1];
        FLINT_ASSERT(start < stop);
        n = stop - start;

        e0 = (Aexps[N*start + off[0]] >> shift[0]) & mask;
        e1 = (Aexps[N*start + off[1]] >> shift[1]) & mask;

        EH->exps[i] = pack_exp2(e0, e1);
        fmpz_mod_poly_fit_length(EH->coeffs + i, n, fpctx);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            fmpz_one(p + j);
            for (k = 2; k < m; k++)
            {
                ulong ei = (Aexps[N*(start + j) + off[k]] >> shift[k]) & mask;
                fmpz_mod_pow_cache_mulpow_ui(p + j, p + j, ei,
                                                  alpha_caches + k - 2, fpctx);
            }
        }
    }

    EH->length = Amarkslen;

    TMP_END;
}


/* return the largest degree */
slong fmpz_mod_polyun_product_roots(
    fmpz_mod_polyun_t M,
    const fmpz_mod_polyun_t H,
    const fmpz_mod_ctx_t ctx)
{
    slong i, max_length = 0;

    fmpz_mod_polyun_fit_length(M, H->length, ctx);
    M->length = H->length;
    for (i = 0; i < H->length; i++)
    {
        slong len = H->coeffs[i].length;
        M->exps[i] = H->exps[i];
        max_length = FLINT_MAX(max_length, len);
        fmpz_mod_poly_product_roots_fmpz_vec(M->coeffs + i,
                                             H->coeffs[i].coeffs, len, ctx);
    }

    return max_length;
}


void fmpz_mod_polyun_zip_start(
    fmpz_mod_polyun_t Z,
    fmpz_mod_polyun_t H,
    slong req_images,
    const fmpz_mod_ctx_t fctx)
{
    slong j;
    fmpz_mod_polyun_fit_length(Z, H->length, fctx);
    Z->length = H->length;
    for (j = 0; j < H->length; j++)
    {
        Z->exps[j] = H->exps[j];
        fmpz_mod_poly_fit_length(Z->coeffs + j, req_images, fctx);
        Z->coeffs[j].length = 0;
    }
}


/*
    return 
        -1: singular
        0:  inconsistent
        1:  success
*/
int _fmpz_mod_zip_vand_solve(
    fmpz * coeffs,             /* in Fp: size mlength */
    const fmpz * monomials,    /* in Fp: size mlength */
    slong mlength,
    const fmpz * evals,        /* in Fp: size elength */
    slong elength,
    const fmpz * master,       /* in Fp: size mlength + 1 */
    fmpz * scratch,            /* in Fp: size mlength */
    const fmpz_mod_ctx_t ctx)
{
    int success;
    slong i, j;
    fmpz_t V, T, S, r;

    FLINT_ASSERT(elength >= mlength);

    fmpz_init(V);
    fmpz_init(T);
    fmpz_init(S);
    fmpz_init(r);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        fmpz_zero(V);
        fmpz_zero(T);
        fmpz_zero(S);
        fmpz_set(r, monomials + i);
        for (j = mlength; j > 0; j--)
        {
            fmpz_mod_addmul(T, master + j, r, T, ctx);
            fmpz_mod_addmul(S, T, r, S, ctx);
            fmpz_mod_addmul(V, V, T, evals + j - 1, ctx);
        }
        /* roots[i] should be a root of master */
    #if FLINT_WANT_ASSERT
        fmpz_mod_addmul(T, master + 0, r, T, ctx);
        FLINT_ASSERT(fmpz_is_zero(T));
    #endif
        fmpz_mod_mul(S, S, r, ctx);
        if (fmpz_is_zero(S))
        {
            success = -1;
            goto cleanup;
        }
        fmpz_mod_inv(S, S, ctx);
        fmpz_mod_mul(coeffs + i, V, S, ctx);
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
        fmpz_mod_pow_ui(scratch + j, monomials + j, mlength, ctx);

    for (i = mlength; i < elength; i++)
    {
        fmpz_zero(V);
        for (j = 0; j < mlength; j++)
        {
            fmpz_mod_mul(scratch + j, scratch + j, monomials + j, ctx);
            fmpz_mod_addmul(V, V, coeffs + j, scratch + j, ctx);
        }

        if (!fmpz_equal(V, evals + i))
        {
            success = 0;
            goto cleanup;
        }
    }

    success = 1;

cleanup:

    fmpz_clear(V);
    fmpz_clear(T);
    fmpz_clear(S);
    fmpz_clear(r);

    return success;
}



void _fmpz_mod_zip_eval_step(
    fmpz_t ev,
    fmpz * cur,            /* in Fp */
    const fmpz * inc,      /* in Fp */
    const fmpz * coeffs,   /* in Fp */
    slong length,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_zero(ev);
    for (i = 0; i < length; i++)
    {
        fmpz_mod_addmul(ev, ev, cur + i, coeffs + i, ctx);
        fmpz_mod_mul(cur + i, cur + i, inc + i, ctx);
    }
}
