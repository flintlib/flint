/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/* try to prove that A is not a square */
static int _is_proved_not_square(
    int count,
    flint_rand_t state,
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t Abits,
    const mpoly_ctx_t mctx,
    const fq_nmod_ctx_t fqctx)
{
    int tries_left, success = 0;
    slong i, N = mpoly_words_per_exp(Abits, mctx);
    fq_nmod_struct eval[1], * alphas, ** alpha_ptrs;
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

    fq_nmod_init(eval, fqctx);

    alphas = (fq_nmod_struct *) TMP_ALLOC(mctx->nvars*sizeof(fq_nmod_struct));
    alpha_ptrs = (fq_nmod_struct **) TMP_ALLOC(mctx->nvars*sizeof(fq_nmod_struct *));
    for (i = 0; i < mctx->nvars; i++)
    {
        alpha_ptrs[i] = alphas + i;
        fq_nmod_init(alphas + i, fqctx);
    }

next_p:

    for (i = 0; i < mctx->nvars; i++)
        fq_nmod_rand(alphas + i, state, fqctx);

    _fq_nmod_mpoly_eval_all_fq_nmod(eval, Acoeffs, Aexps, Alen, Abits,
                                                      alpha_ptrs, mctx, fqctx);

    success = !fq_nmod_is_square(eval, fqctx);

    if (!success && --tries_left >= 0)
        goto next_p;

    fq_nmod_clear(eval, fqctx);
    for (i = 0; i < mctx->nvars; i++)
        fq_nmod_clear(alphas + i, fqctx);

cleanup:

    TMP_END;

    return success;
}

static int n_fq_sqrt(mp_limb_t * q, const mp_limb_t * a, const fq_nmod_ctx_t ctx)
{
    int res;
    fq_nmod_t t;
    fq_nmod_init(t, ctx);
    n_fq_get_fq_nmod(t, a, ctx);
    res = fq_nmod_sqrt(t, t, ctx);
    n_fq_set_fq_nmod(q, t, ctx);
    fq_nmod_clear(t, ctx);
    return res;
}


static int _fq_nmod_mpoly_sqrt_heap(
    fq_nmod_mpoly_t Q,
    const mp_limb_t * Acoeffs,
    const ulong * Aexps,
    slong Alen,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx,
    const fq_nmod_ctx_t fqctx)
{
    slong d = fq_nmod_ctx_degree(fqctx);
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
    mp_limb_t * Qcoeffs = Q->coeffs;
    ulong * Qexps = Q->exps;
    ulong * exp, * exp3;
    ulong * exps[64];
    ulong ** exp_list;
    slong exp_next;
    ulong mask;
    mp_limb_t * t, * t2, * lc_inv;
    int lt_divides, halves;
    flint_rand_t heuristic_state;
    int heuristic_count = 0;
    TMP_INIT;

    TMP_START;

    t = (mp_limb_t *) TMP_ALLOC(13*d*sizeof(mp_limb_t));
    t2 = t + 6*d;
    lc_inv = t2 + 6*d;

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
    _fq_nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc, d,
                              &Qexps, &Q->exps_alloc, N, Qlen + 1);

    if (!n_fq_sqrt(Qcoeffs + d*0, Acoeffs + d*0, fqctx))
        goto not_sqrt;

    Qlen = 1;

    /* precompute leading cofficient info */
    _n_fq_add(t2, Qcoeffs + d*0, Qcoeffs + d*0, d, fqctx->mod);
    _n_fq_inv(lc_inv, t2, fqctx, t);

    if (bits <= FLINT_BITS)
        halves = mpoly_monomial_halves(Qexps + 0, Aexps + 0, N, mask);
    else
        halves = mpoly_monomial_halves_mp(Qexps + 0, Aexps + 0, N, bits);

    if (!halves)
        goto not_sqrt; /* exponent is not square */

    /* optimisation, compute final term */
    {
        if (!n_fq_sqrt(t, Acoeffs + d*(Alen - 1), fqctx))
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
        _fq_nmod_mpoly_fit_length(&Qcoeffs, &Q->coeffs_alloc, d,
                                  &Qexps, &Q->exps_alloc, N, Qlen + 1);

        if (heap_len > 1 && Ai < Alen &&
            mpoly_monomial_equal(Aexps + N*Ai, heap[1].exp, N))
        {
            /* take from both A and heap */
            mpoly_monomial_set(exp, Aexps + N*Ai, N);
            _n_fq_set(Qcoeffs + d*Qlen, Acoeffs + d*Ai, d);
            Ai++;
        }
        else if (heap_len > 1 && (Ai >= Alen || mpoly_monomial_lt(
                                       Aexps + N*Ai, heap[1].exp, N, cmpmask)))
        {
            /* take only from heap */
            mpoly_monomial_set(exp, heap[1].exp, N);
            _n_fq_zero(Qcoeffs + d*Qlen, d);

            if (bits <= FLINT_BITS ? mpoly_monomial_overflows(exp, N, mask)
                                   : mpoly_monomial_overflows_mp(exp, N, bits))
                goto not_sqrt;
        }
        else
        {
            FLINT_ASSERT(Ai < Alen);

            /* take only from A */
            mpoly_monomial_set(exp, Aexps + N*Ai, N);
            _n_fq_set(Qcoeffs + d*Qlen, Acoeffs + d*Ai, d);
            Ai++;

            goto skip_heap;
        }

        _nmod_vec_zero(t, 6*d);
        _nmod_vec_zero(t2, 6*d);
        /* TODO lazy_size = _n_fq_dot_lazy_size(heap_alloc + 1, fqctx) */
        {
            do {
                exp_list[--exp_next] = heap[1].exp;
                x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);
                do {
                    mp_limb_t * dest;
                    *store++ = x->i;
                    *store++ = x->j;
                    dest = (x->i != x->j) ? t2 : t;
                    _n_fq_madd2(dest, Qcoeffs + d*x->i,
                                      Qcoeffs + d*x->j, fqctx, dest + 2*d);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));
            _nmod_vec_add(t, t, t2, 2*d, fqctx->mod);
            _nmod_vec_add(t, t, t2, 2*d, fqctx->mod);
        }

        _n_fq_reduce2(t2, t, fqctx, t + 2*d);
        _nmod_vec_sub(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, t2, d, fqctx->mod);

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

        if (_n_fq_is_zero(Qcoeffs + d*Qlen, d))
            continue;

    skip_heap:

        if (bits <= FLINT_BITS)
            lt_divides = mpoly_monomial_divides(Qexps + N*Qlen,
                                                exp, Qexps + N*0, N, mask);
        else
            lt_divides = mpoly_monomial_divides_mp(Qexps + N*Qlen,
                                                exp, Qexps + N*0, N, bits);
        if (!lt_divides)
            goto not_sqrt;

        _n_fq_mul(Qcoeffs + d*Qlen, Qcoeffs + d*Qlen, lc_inv, fqctx, t);

        if (Qlen >= heap_alloc)
        {
            /* run some tests if the square root is getting long */
            if (Qlen > Alen && _is_proved_not_square(
                                          ++heuristic_count, heuristic_state,
                                      Acoeffs, Aexps, Alen, bits, mctx, fqctx))
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

int fq_nmod_mpoly_sqrt_heap(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong lenq_est;

    if ((ctx->fqctx->mod.n % 2) == 0)
    {
        slong d = fq_nmod_ctx_degree(ctx->fqctx);
        flint_bitcnt_t bits = A->bits;
        mp_limb_t * Aexps = A->exps;
        slong Alen = A->length;
        slong i, j, N = mpoly_words_per_exp(bits, ctx->minfo);
        ulong mask = (bits <= FLINT_BITS) ? mpoly_overflow_mask_sp(bits) : 0;
        mp_limb_t * t;

        if (Q != A)
            fq_nmod_mpoly_fit_length_reset_bits(Q, Alen, bits, ctx);

        for (i = 0; i < Alen; i++)
        {
            if (bits <= FLINT_BITS ?
                !mpoly_monomial_halves(Q->exps + N*i, Aexps + N*i, N, mask) :
                !mpoly_monomial_halves_mp(Q->exps + N*i, Aexps + N*i, N, bits))
            {
                Q->length = 0;
                return 0;
            }
        }

        t = FLINT_ARRAY_ALLOC(N_FQ_MUL_ITCH*d, mp_limb_t);

        for (i = 0; i < Alen; i++)
        {
            _n_fq_set(Q->coeffs + d*i, A->coeffs + d*i, d);
            for (j = 1; j < d; j++)
                _n_fq_mul(Q->coeffs + d*i, Q->coeffs + d*i, Q->coeffs + d*i,
                                                                ctx->fqctx, t);
        }

        flint_free(t);

        Q->length = Alen;
        return 1;
    }

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        fq_nmod_mpoly_zero(Q, ctx);
        return 1;
    }

    lenq_est = n_sqrt(A->length);

    if (Q == A)
    {
        fq_nmod_mpoly_t T;
        fq_nmod_mpoly_init3(T, lenq_est, A->bits, ctx);
        success = _fq_nmod_mpoly_sqrt_heap(T, A->coeffs, A->exps, A->length,
                                              A->bits, ctx->minfo, ctx->fqctx);
        fq_nmod_mpoly_swap(Q, T, ctx);
        fq_nmod_mpoly_clear(T, ctx);
    }
    else
    {
        fq_nmod_mpoly_fit_length_reset_bits(Q, lenq_est, A->bits, ctx);
        success = _fq_nmod_mpoly_sqrt_heap(Q, A->coeffs, A->exps, A->length,
                                              A->bits, ctx->minfo, ctx->fqctx);
    }

    return success;
}

