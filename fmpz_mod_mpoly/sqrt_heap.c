/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

/* try to prove that A is not a square */
static int _is_proved_not_square(
    int count,
    flint_rand_t state,
    const fmpz * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fctx)
{
    int tries_left, success = 0;
    slong i, N = mpoly_words_per_exp(Abits, mctx);
    fmpz eval[1], * alphas;
    ulong * t;
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);

    TMP_START;
    t = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    if (count == 1)
    {
        success = mpoly_is_proved_not_square(Aexps, Alen, Abits, N, t);
        if (success)
            goto cleanup;
    }

    tries_left = 3*count;

    fmpz_init(eval);
    alphas = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (i = 0; i < mctx->nvars; i++)
        fmpz_init(alphas + i);

next_p:

    for (i = 0; i < mctx->nvars; i++)
        fmpz_randm(alphas + i, state, fmpz_mod_ctx_modulus(fctx));

    _fmpz_mod_mpoly_eval_all_fmpz_mod(eval, Acoeffs, Aexps, Alen, Abits,
                                                           alphas, mctx, fctx);

    success = fmpz_jacobi(eval, fmpz_mod_ctx_modulus(fctx)) < 0;

    if (!success && --tries_left >= 0)
        goto next_p;

    fmpz_clear(eval);
    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(alphas + i);

cleanup:

    TMP_END;

    return success;
}


static int _fmpz_mod_mpoly_sqrt_heap1(
    fmpz_mod_mpoly_t Q,
    const fmpz * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fctx)
{
    slong i, j, Qlen, Ai;
    slong next_loc, heap_len = 1, heap_alloc;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain_nodes[64];
    mpoly_heap_t ** chain;
    slong exp_alloc;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    ulong mask, exp, exp3 = 0;
    ulong cmpmask;
    mpz_t acc, acc2, modulus;
    fmpz zero = 0;
    const fmpz * s;
    fmpz_t lc_inv;
    int lt_divides;
    flint_rand_t heuristic_state;
    int heuristic_count = 0;

    fmpz_init(lc_inv);
    mpz_init(modulus);
    mpz_init(acc);
    mpz_init(acc2);
    fmpz_get_mpz(modulus, fmpz_mod_ctx_modulus(fctx));

    FLINT_ASSERT(mpoly_words_per_exp(bits, mctx) == 1);
    mpoly_get_cmpmask(&cmpmask, 1, bits, mctx);

    flint_randinit(heuristic_state);

    /* alloc array of heap nodes which can be chained together */
    next_loc = 2*n_sqrt(Alen) + 4;   /* something bigger than heap can ever be */
    heap_alloc = next_loc - 3;
    heap = (mpoly_heap1_s *) flint_malloc((heap_alloc + 1)*sizeof(mpoly_heap1_s));
    chain_nodes[0] = (mpoly_heap_t *) flint_malloc(heap_alloc*sizeof(mpoly_heap_t));
    chain = (mpoly_heap_t **) flint_malloc(heap_alloc*sizeof(mpoly_heap_t*));
    store = store_base = (slong *) flint_malloc(2*heap_alloc*sizeof(mpoly_heap_t *));

    for (i = 0; i < heap_alloc; i++)
       chain[i] = chain_nodes[0] + i;

    exp_alloc = 1;

    mask = mpoly_overflow_mask_sp(bits);

    /* "insert" (-1, 1, Aexps[1]) into "heap" */
    Ai = 1;

    /* compute first term */
    Qlen = 0;
    _fmpz_mod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                               &Qexps, &Q->exps_alloc, 1, Qlen + 1);

    if (!fmpz_sqrtmod(Qcoeffs + 0, Acoeffs + 0, fmpz_mod_ctx_modulus(fctx)))
        goto not_sqrt;

    Qlen = 1;

    /* precompute leading cofficient info */
    fmpz_mod_add(lc_inv, Qcoeffs + 0, Qcoeffs + 0, fctx);
    fmpz_mod_inv(lc_inv, lc_inv, fctx);

    if (!mpoly_monomial_halves1(Qexps + 0, Aexps[0], mask))
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final exponent */
    {
        if (fmpz_jacobi(Acoeffs + Alen - 1, fmpz_mod_ctx_modulus(fctx)) < 0)
            goto not_sqrt;

        if (!mpoly_monomial_halves1(&exp3, Aexps[Alen - 1], mask))
            goto not_sqrt; /* exponent is not square */

        exp3 += Qexps[0]; /* overflow not possible */
    }

    while (heap_len > 1 || Ai < Alen)
    {
        _fmpz_mod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                                   &Qexps, &Q->exps_alloc, 1, Qlen + 1);

        if (heap_len > 1 && Ai < Alen && Aexps[Ai] == heap[1].exp)
        {
            /* take from both A and heap */
            exp = Aexps[Ai];
            s = Acoeffs + Ai;
            Ai++;
        }
        else if (heap_len > 1 && (Ai >= Alen ||
                          mpoly_monomial_gt1(heap[1].exp, Aexps[Ai], cmpmask)))
        {
            /* take only from heap */
            exp = heap[1].exp;
            s = &zero;
            if (mpoly_monomial_overflows1(exp, mask))
                goto not_sqrt;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen);

            /* take only from A */
            exp = Aexps[Ai];
            s = Acoeffs + Ai;
            Ai++;

            goto skip_heap;
        }

        /* total is always acc + 2*acc2 */
        mpz_set_ui(acc, 0);
        mpz_set_ui(acc2, 0);
        do {
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);
            do {
                mpz_ptr t;
                fmpz Qi, Qj;

                *store++ = x->i;
                *store++ = x->j;

                Qi = Qcoeffs[x->i];
                Qj = Qcoeffs[x->j];
                t = (x->i != x->j) ? acc2 : acc;

                if (COEFF_IS_MPZ(Qi) && COEFF_IS_MPZ(Qj))
                {
                    mpz_addmul(t, COEFF_TO_PTR(Qi), COEFF_TO_PTR(Qj));
                }
                else if (COEFF_IS_MPZ(Qi) && !COEFF_IS_MPZ(Qj))
                {
                    flint_mpz_addmul_ui(t, COEFF_TO_PTR(Qi), Qj);
                }
                else if (!COEFF_IS_MPZ(Qi) && COEFF_IS_MPZ(Qj))
                {
                    flint_mpz_addmul_ui(t, COEFF_TO_PTR(Qj), Qi);
                }
                else
                {
                    ulong pp1, pp0;
                    umul_ppmm(pp1, pp0, Qcoeffs[x->i], Qcoeffs[x->j]);
                    flint_mpz_add_uiui(t, t, pp1, pp0);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && heap[1].exp == exp);

        mpz_addmul_ui(acc, acc2, 2);
        mpz_tdiv_qr(acc2, _fmpz_promote(Qcoeffs + Qlen), acc, modulus);
        _fmpz_demote_val(Qcoeffs + Qlen);

        fmpz_mod_sub(Qcoeffs + Qlen, s, Qcoeffs + Qlen, fctx);
        s = Qcoeffs + Qlen;

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right */
            if (j < i)
            {
                x = chain[i];
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                _mpoly_heap_insert1(heap, Qexps[x->i] + Qexps[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
            }
        }

    skip_heap:

        fmpz_mod_mul(Qcoeffs + Qlen, s, lc_inv, fctx);
        if (fmpz_is_zero(Qcoeffs + Qlen))
            continue;

        lt_divides = mpoly_monomial_divides1(Qexps + Qlen, exp, Qexps[0], mask);
        if (!lt_divides)
            goto not_sqrt;

        if (Qlen >= heap_alloc)
        {
            /* run some tests if the square root is getting long */
            if (Qlen > Alen && _is_proved_not_square(
                                        ++heuristic_count, heuristic_state,
                                       Acoeffs, Aexps, Alen, bits, mctx, fctx))
            {
                goto not_sqrt;
            }

            heap_alloc *= 2;
            heap = (mpoly_heap1_s *) flint_realloc(heap, (heap_alloc + 1)*sizeof(mpoly_heap1_s));
            chain_nodes[exp_alloc] = (mpoly_heap_t *) flint_malloc((heap_alloc/2)*sizeof(mpoly_heap_t));
            chain = (mpoly_heap_t **) flint_realloc(chain, heap_alloc*sizeof(mpoly_heap_t*));
            store = store_base = (slong *) flint_realloc(store_base, 2*heap_alloc*sizeof(mpoly_heap_t *));
            for (i = 0; i < heap_alloc/2; i++)
                chain[i + heap_alloc/2] = chain_nodes[exp_alloc] + i;
            exp_alloc++;
        }

        /* put (Qlen, 1) in heap */
        i = Qlen;
        x = chain[i];
        x->i = i;
        x->j = 1;
        x->next = NULL;

        _mpoly_heap_insert1(heap, Qexps[i] + Qexps[1], x,
                                                &next_loc, &heap_len, cmpmask);

        Qlen++;
    }

cleanup:

    flint_randclear(heuristic_state);

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    fmpz_clear(lc_inv);
    mpz_clear(modulus);
    mpz_clear(acc);
    mpz_clear(acc2);

    flint_free(heap);
    flint_free(chain);
    flint_free(store_base);
    for (i = 0; i < exp_alloc; i++)
        flint_free(chain_nodes[i]);

    return Qlen > 0;

not_sqrt:

    Qlen = 0;
    goto cleanup;
}

static int _fmpz_mod_mpoly_sqrt_heap(
    fmpz_mod_mpoly_t Q,
    const fmpz * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx,
    const fmpz_mod_ctx_t fctx)
{
    slong N = mpoly_words_per_exp(bits, mctx);
    ulong * cmpmask;
    slong i, j, Qlen, Ai;
    slong next_loc;
    slong heap_len = 1, heap_alloc;
    int exp_alloc;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain_nodes[64];
    mpoly_heap_t ** chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    ulong * exp, * exp3;
    ulong * exps[64];
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    mpz_t acc, acc2, modulus;
    fmpz zero = 0;
    const fmpz * s;
    fmpz_t lc_inv;
    int halves, lt_divides;
    flint_rand_t heuristic_state;
    int heuristic_count = 0;
    TMP_INIT;

    if (N == 1)
        return _fmpz_mod_mpoly_sqrt_heap1(Q, Acoeffs, Aexps, Alen, bits,
                                                                   mctx, fctx);

    fmpz_init(lc_inv);
    mpz_init(modulus);
    mpz_init(acc);
    mpz_init(acc2);
    fmpz_get_mpz(modulus, fmpz_mod_ctx_modulus(fctx));

    TMP_START;

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, mctx);

    flint_randinit(heuristic_state);

    /* alloc array of heap nodes which can be chained together */
    next_loc = 2*sqrt(Alen) + 4;   /* something bigger than heap can ever be */
    heap_alloc = next_loc - 3;
    heap = (mpoly_heap_s *) flint_malloc((heap_alloc + 1)*sizeof(mpoly_heap_s));
    chain_nodes[0] = (mpoly_heap_t *) flint_malloc(heap_alloc*sizeof(mpoly_heap_t));
    chain = (mpoly_heap_t **) flint_malloc(heap_alloc*sizeof(mpoly_heap_t*));
    store = store_base = (slong *) flint_malloc(2*heap_alloc*sizeof(mpoly_heap_t *));

    for (i = 0; i < heap_alloc; i++)
       chain[i] = chain_nodes[0] + i;

    /* array of exponent vectors, each of "N" words */
    exps[0] = (ulong *) flint_malloc(heap_alloc*N*sizeof(ulong));
    exp_alloc = 1;
    /* list of pointers to available exponent vectors */
    exp_list = (ulong **) flint_malloc(heap_alloc*sizeof(ulong *));
    /* space to save copy of current exponent vector */
    exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    /* final exponent */
    exp3 = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    /* set up list of available exponent vectors */
    exp_next = 0;
    for (i = 0; i < heap_alloc; i++)
        exp_list[i] = exps[0] + i*N;

    mask = (bits <= FLINT_BITS) ? mpoly_overflow_mask_sp(bits) : 0;

    /* "insert" (-1, 1, Aexps[0]) into "heap" */
    Ai = 1;

    /* compute first term */
    Qlen = 0;
    _fmpz_mod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                               &Qexps, &Q->exps_alloc, 1, Qlen + 1);

FLINT_ASSERT(Alen > 0);
FLINT_ASSERT(!fmpz_is_zero(Acoeffs + 0));
FLINT_ASSERT(fmpz_mod_is_canonical(Acoeffs + 0, fctx));

    if (!fmpz_sqrtmod(Qcoeffs + 0, Acoeffs + 0, fmpz_mod_ctx_modulus(fctx)))
        goto not_sqrt;

    Qlen = 1;

    /* precompute leading cofficient info */
    fmpz_mod_add(lc_inv, Qcoeffs + 0, Qcoeffs + 0, fctx);
    fmpz_mod_inv(lc_inv, lc_inv, fctx);

    if (bits <= FLINT_BITS)
        halves = mpoly_monomial_halves(Qexps + 0, Aexps + 0, N, mask);
    else
        halves = mpoly_monomial_halves_mp(Qexps + 0, Aexps + 0, N, bits);

    if (!halves)
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final term */
    {
        if (fmpz_jacobi(Acoeffs + Alen - 1, fmpz_mod_ctx_modulus(fctx)) < 0)
            goto not_sqrt;

        if (bits <= FLINT_BITS)
            halves = mpoly_monomial_halves(exp3, Aexps + (Alen - 1)*N, N, mask);
        else
            halves = mpoly_monomial_halves_mp(exp3, Aexps + (Alen - 1)*N, N, bits);

        if (!halves)
            goto not_sqrt; /* exponent is not square */

        if (bits <= FLINT_BITS)
            mpoly_monomial_add(exp3, exp3, Qexps + 0, N);
        else
            mpoly_monomial_add_mp(exp3, exp3, Qexps + 0, N);
    }

    while (heap_len > 1 || Ai < Alen)
    {
        _fmpz_mod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc,
                                   &Qexps, &Q->exps_alloc, N, Qlen + 1);

        if (heap_len > 1 && Ai < Alen &&
            mpoly_monomial_equal(Aexps + N*Ai, heap[1].exp, N))
        {
            /* take from both A and heap */
            mpoly_monomial_set(exp, Aexps + N*Ai, N);
            s = Acoeffs + Ai;
            Ai++;
        }
        else if (heap_len > 1 && (Ai >= Alen || mpoly_monomial_lt(
                                       Aexps + N*Ai, heap[1].exp, N, cmpmask)))
        {
            /* take only from heap */
            mpoly_monomial_set(exp, heap[1].exp, N);
            s = &zero;
            if (bits <= FLINT_BITS ? mpoly_monomial_overflows(exp, N, mask)
                                   : mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_sqrt;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen);

            /* take only from A */
            mpoly_monomial_set(exp, Aexps + N*Ai, N);
            s = Acoeffs + Ai;
            Ai++;

            goto skip_heap;
        }

        /* total is always acc + 2*acc2 */
        mpz_set_ui(acc, 0);
        mpz_set_ui(acc2, 0);
        do {
            exp_list[--exp_next] = heap[1].exp;
            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
            do {
                mpz_ptr t;
                fmpz Qi, Qj;

                *store++ = x->i;
                *store++ = x->j;

                Qi = Qcoeffs[x->i];
                Qj = Qcoeffs[x->j];
                t = (x->i != x->j) ? acc2 : acc;

                if (COEFF_IS_MPZ(Qi) && COEFF_IS_MPZ(Qj))
                {
                    mpz_addmul(t, COEFF_TO_PTR(Qi), COEFF_TO_PTR(Qj));
                }
                else if (COEFF_IS_MPZ(Qi) && !COEFF_IS_MPZ(Qj))
                {
                    flint_mpz_addmul_ui(t, COEFF_TO_PTR(Qi), Qj);
                }
                else if (!COEFF_IS_MPZ(Qi) && COEFF_IS_MPZ(Qj))
                {
                    flint_mpz_addmul_ui(t, COEFF_TO_PTR(Qj), Qi);
                }
                else
                {
                    ulong pp1, pp0;
                    umul_ppmm(pp1, pp0, Qcoeffs[x->i], Qcoeffs[x->j]);
                    flint_mpz_add_uiui(t, t, pp1, pp0);
                }
            } while ((x = x->next) != NULL);
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        mpz_addmul_ui(acc, acc2, 2);
        mpz_tdiv_qr(acc2, _fmpz_promote(Qcoeffs + Qlen), acc, modulus);
        _fmpz_demote_val(Qcoeffs + Qlen);

        fmpz_mod_sub(Qcoeffs + Qlen, s, Qcoeffs + Qlen, fctx);
        s = Qcoeffs + Qlen;

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right */
            if (j < i)
            {
                x = chain[i];
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Qexps + N*x->i,
                                                            Qexps + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Qexps + N*x->i,
                                                            Qexps + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }

    skip_heap:

        fmpz_mod_mul(Qcoeffs + Qlen, s, lc_inv, fctx);
        if (fmpz_is_zero(Qcoeffs + Qlen))
            continue;

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(Qexps + N*Qlen,
                                                exp, Qexps + N*0, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(Qexps + N*Qlen,
                                                exp, Qexps + N*0, N, bits);
        if (!lt_divides)
            goto not_sqrt;

        if (Qlen >= heap_alloc)
        {
            /* run some tests if the square root is getting long */
            if (Qlen > Alen && _is_proved_not_square(
                                          ++heuristic_count, heuristic_state,
                                       Acoeffs, Aexps, Alen, bits, mctx, fctx))
            {
                goto not_sqrt;
            }

            heap_alloc *= 2;
            heap = (mpoly_heap_s *) flint_realloc(heap, (heap_alloc + 1)*sizeof(mpoly_heap_s));
            chain_nodes[exp_alloc] = (mpoly_heap_t *) flint_malloc((heap_alloc/2)*sizeof(mpoly_heap_t));
            chain = (mpoly_heap_t **) flint_realloc(chain, heap_alloc*sizeof(mpoly_heap_t*));
            store = store_base = (slong *) flint_realloc(store_base, 2*heap_alloc*sizeof(mpoly_heap_t *));
            exps[exp_alloc] = (ulong *) flint_malloc((heap_alloc/2)*N*sizeof(ulong));
            exp_list = (ulong **) flint_realloc(exp_list, heap_alloc*sizeof(ulong *));
            for (i = 0; i < heap_alloc/2; i++)
            {
               chain[i + heap_alloc/2] = chain_nodes[exp_alloc] + i;
               exp_list[i + heap_alloc/2] = exps[exp_alloc] + i*N;
            }
            exp_alloc++;
        }

        /* put (Qlen, 1) in heap */
        i = Qlen;
        x = chain[i];
        x->i = i;
        x->j = 1;
        x->next = NULL;

        if (bits <= FLINT_BITS)
            mpoly_monomial_add(exp_list[exp_next], Qexps + x->i*N,
                                                      Qexps + x->j*N, N);
        else
            mpoly_monomial_add_mp(exp_list[exp_next], Qexps + x->i*N,
                                                         Qexps + x->j*N, N);

        exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                         &next_loc, &heap_len, N, cmpmask);

        Qlen++;
    }

cleanup:

    flint_randclear(heuristic_state);

    Q->coeffs = Qcoeffs;
    Q->exps = Qexps;
    Q->length = Qlen;

    fmpz_clear(lc_inv);
    mpz_clear(modulus);
    mpz_clear(acc);
    mpz_clear(acc2);

    flint_free(heap);
    flint_free(chain);
    flint_free(store_base);
    flint_free(exp_list);
    for (i = 0; i < exp_alloc; i++)
    {
        flint_free(exps[i]);
        flint_free(chain_nodes[i]);
    }

    TMP_END;

    return Qlen > 0;

not_sqrt:
    Qlen = 0;
    goto cleanup;
}


int fmpz_mod_mpoly_sqrt_heap(
    fmpz_mod_mpoly_t Q,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong lenq_est;

    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        fmpz_mod_mpoly_zero(Q, ctx);
        return 1;
    }

    if (fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)))
    {
        nmod_mpoly_ctx_t nctx;
        nmod_mpoly_t nQ, nA;

        nctx->minfo[0] = ctx->minfo[0];
        nmod_init(&nctx->mod, fmpz_get_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)));
        nmod_mpoly_init(nQ, nctx);
        nmod_mpoly_init(nA, nctx);

        _fmpz_mod_mpoly_get_nmod_mpoly(nA, nctx, A, ctx);
        success = nmod_mpoly_sqrt_heap(nQ, nA, nctx);
        _fmpz_mod_mpoly_set_nmod_mpoly(Q, ctx, nQ, nctx);

        nmod_mpoly_clear(nA, nctx);
        nmod_mpoly_clear(nQ, nctx);

        return success;
    }

    lenq_est = n_sqrt(A->length);

    if (Q == A)
    {
        fmpz_mod_mpoly_t T;
        fmpz_mod_mpoly_init3(T, lenq_est, A->bits, ctx);
        success = _fmpz_mod_mpoly_sqrt_heap(T, A->coeffs, A->exps, A->length,
                                             A->bits, ctx->minfo, ctx->ffinfo);
        fmpz_mod_mpoly_swap(Q, T, ctx);
        fmpz_mod_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mod_mpoly_fit_length_reset_bits(Q, lenq_est, A->bits, ctx);
        success = _fmpz_mod_mpoly_sqrt_heap(Q, A->coeffs, A->exps, A->length,
                                             A->bits, ctx->minfo, ctx->ffinfo);
    }

    return success;
}
